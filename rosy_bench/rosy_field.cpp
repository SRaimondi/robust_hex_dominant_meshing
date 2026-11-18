// rosy_field -- headless driver for the RoSy (3D frame field) stage of
//   Gao, Jakob, Tarini, Panozzo, "Robust Hex-Dominant Mesh Generation using
//   Field-Guided Polyhedral Agglomeration", SIGGRAPH 2017.
//
// Pipeline (identical calls to src/batch.cpp + Optimizer::run(), executed synchronously):
//   MultiResolutionHierarchy::load        read the surface OBJ
//   MultiResolutionHierarchy::tet_meshing TetGen ("pqm" with a constant-size background mesh)
//   MultiResolutionHierarchy::build       normals, adjacency, multi-resolution hierarchy, random init
//   smoothOrientationsTet / prolongOrientations, coarse -> fine, N iterations per level
//
// Outputs (for -o <prefix>):
//   <prefix>_field.txt         per-vertex quaternions + edge connectivity (+ tets), see write_field()
//   <prefix>_singularities.obj singularity graph as line segments
//   <prefix>_stats.json        timings, sizes, status (also written on failure)
//   <prefix>_tet.mesh          with --save-tet: Medit tet mesh, written right after tet meshing

#include "hierarchy.h"
#include "orientations.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <thread>

#if defined(TBB_VERSION_MAJOR) && TBB_VERSION_MAJOR >= 2021
#include <tbb/global_control.h>
#define RHDM_ONETBB 1
#endif

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

std::string json_escape(const std::string &s) {
    std::string out;
    for (char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\t': out += "\\t"; break;
            default:
                if ((unsigned char) c < 0x20) { char buf[8]; snprintf(buf, sizeof(buf), "\\u%04x", c); out += buf; }
                else out += c;
        }
    }
    return out;
}

// Minimal ordered JSON object writer (flat key/value pairs).
struct Json {
    std::vector<std::pair<std::string, std::string>> kv;
    void num(const std::string &k, double v) { std::ostringstream o; o.precision(10); o << v; kv.emplace_back(k, o.str()); }
    void integer(const std::string &k, long long v) { kv.emplace_back(k, std::to_string(v)); }
    void str(const std::string &k, const std::string &v) { kv.emplace_back(k, "\"" + json_escape(v) + "\""); }
    void raw(const std::string &k, const std::string &v) { kv.emplace_back(k, v); }
    bool write(const std::string &path) const {
        std::ofstream f(path);
        if (!f) return false;
        f << "{\n";
        for (size_t i = 0; i < kv.size(); ++i)
            f << "  \"" << kv[i].first << "\": " << kv[i].second << (i + 1 < kv.size() ? ",\n" : "\n");
        f << "}\n";
        return (bool) f;
    }
};

void usage(const char *argv0) {
    fprintf(stderr,
        "Usage: %s -i <input.obj> -o <output_prefix> [options]\n"
        "\n"
        "Options:\n"
        "  --tet-len-factor <k>  target tet edge length = k * (average input surface edge length).\n"
        "                        Default 1.05, which is what the original batch/GUI default resolves to.\n"
        "  --tet-len <L>         absolute target tet edge length (overrides --tet-len-factor).\n"
        "  --iters <n>           RoSy iterations per hierarchy level (default 200, as hard-coded\n"
        "                        in Optimizer::setOptimizeOrientations for tet meshes).\n"
        "  --threads <n>         number of TBB worker threads (default: all cores).\n"
        "  --no-hierarchy        smooth only on the finest level (Optimizer::setHierarchy(false)).\n"
        "  --save-tet            write the TetGen tet mesh to <prefix>_tet.mesh (Medit: vertices,\n"
        "                        boundary triangles, tets) right after tet meshing succeeds.\n",
        argv0);
}

enum ExitCode { OK = 0, BAD_ARGS = 2, LOAD_FAILED = 3, TET_FAILED = 4, FIELD_FAILED = 5, WRITE_FAILED = 6 };

// ------------------------------------------------------------------ field file
//
//   # comments ...
//   FRAMEFIELD 1
//   QUATERNION xyzw            unit quaternion, Hamilton convention; the frame axes are the
//                              columns of the rotation matrix R(q) (same as Quaternion::toMatrix)
//   VERTICES <nV>
//   x y z qx qy qz qw b        b = 1 for boundary (surface) vertices, else 0
//   EDGES <nE>
//   i j                        0-based, i < j; the neighbourhood graph used by the smoother
//   TETS <nT>                  (optional for other producers; 0 allowed)
//   a b c d                    0-based
bool write_field(const std::string &path, const MultiResolutionHierarchy &mRes) {
    const MatrixXf &V = mRes.mV[0];
    const MatrixXf &Q = mRes.mQ[0];
    const SMatrix &L = mRes.mL[0];
    const MatrixXu &T = mRes.mT;
    const std::vector<bool> &B = mRes.nV_boundary_flag[0];

    size_t nE = 0;
    for (int k = 0; k < L.outerSize(); ++k)
        for (SMatrix::InnerIterator it(L, k); it; ++it)
            if (it.col() > it.row()) ++nE;

    FILE *f = fopen(path.c_str(), "w");
    if (!f) return false;
    fprintf(f, "# Frame field exported by rosy_field (robust_hex_dominant_meshing, Gao et al. 2017)\n");
    fprintf(f, "# VERTICES: x y z qx qy qz qw boundary | EDGES: i j (0-based) | TETS: a b c d (0-based)\n");
    fprintf(f, "# Frame axes = columns of the rotation matrix of q (Hamilton, unit, stored x y z w).\n");
    fprintf(f, "FRAMEFIELD 1\nQUATERNION xyzw\n");
    fprintf(f, "VERTICES %lld\n", (long long) V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        Quaternion q = Quaternion(Q.col(i)).normalized();
        fprintf(f, "%.9g %.9g %.9g %.9g %.9g %.9g %.9g %d\n",
                V(0, i), V(1, i), V(2, i), q.x(), q.y(), q.z(), q.w(), B[i] ? 1 : 0);
    }
    fprintf(f, "EDGES %zu\n", nE);
    for (int k = 0; k < L.outerSize(); ++k)
        for (SMatrix::InnerIterator it(L, k); it; ++it)
            if (it.col() > it.row())
                fprintf(f, "%d %d\n", (int) it.row(), (int) it.col());
    fprintf(f, "TETS %lld\n", (long long) T.cols());
    for (int t = 0; t < T.cols(); ++t)
        fprintf(f, "%u %u %u %u\n", T(0, t), T(1, t), T(2, t), T(3, t));
    bool ok = ferror(f) == 0;
    ok &= fclose(f) == 0;
    return ok;
}

// ------------------------------------------------------------------ singularities
//
// Same criterion as MultiResolutionHierarchy::detectOrientationSingularitiesTet():
// a tet face is singular if composing the matched octahedral symmetries
// (AxisToAxisRot) around its three edges is not the identity. For every
// (singular face, incident tet) pair a segment face-centroid -> tet-centroid is
// emitted; the segments of neighbouring tets share endpoints, so the output is a
// connected graph tracing the singular curves. Vertices are de-duplicated.
//
// OBJ groups, by the holonomy's conjugacy class:
//   rot90_pos / rot90_neg  +-90 deg about a frame axis (regular singular lines; the sign is the
//                          original code's red/blue rule: rotation axis . face normal > 0)
//   rot180                 180 deg (about a frame axis or a face diagonal)
//   rot120                 120 deg about a body diagonal
struct SingularityStats { size_t faces = 0, segments = 0, rot90_pos = 0, rot90_neg = 0, rot180 = 0, rot120 = 0; };

bool write_singularities(const std::string &path, const MultiResolutionHierarchy &mRes, SingularityStats &st) {
    const MatrixXf &V = mRes.mV[0];
    const MatrixXf &Q = mRes.mQ[0];
    const MatrixXu &T = mRes.mT;

    std::map<std::array<uint32_t, 3>, uint32_t> faceVertex;  // sorted face -> obj vertex
    std::map<uint32_t, uint32_t> tetVertex;                  // tet -> obj vertex
    std::vector<Vector3f> verts;
    std::vector<std::pair<uint32_t, uint32_t>> groups[4];    // rot90_pos, rot90_neg, rot180, rot120
    std::set<std::array<uint32_t, 3>> singularFaces;

    for (uint32_t t = 0; t < (uint32_t) T.cols(); ++t) {
        for (auto face : tet_faces) {
            AxisToAxisRot rot1, rot2, rot3;
            for (int i = 0; i < 3; ++i) {
                int i0 = T(face[i], t), i1 = T(face[(i + 1) % 3], t);
                rot1 *= AxisToAxisRot(Quaternion(Q.col(i0)), Quaternion(Q.col(i1)));
            }
            if (rot1.isIdentity())
                continue;

            std::array<uint32_t, 3> key = { T(face[0], t), T(face[1], t), T(face[2], t) };
            std::sort(key.begin(), key.end());
            singularFaces.insert(key);

            auto fit = faceVertex.find(key);
            if (fit == faceVertex.end()) {
                Vector3f c = (V.col(key[0]) + V.col(key[1]) + V.col(key[2])) / 3.f;
                fit = faceVertex.emplace(key, (uint32_t) verts.size()).first;
                verts.push_back(c);
            }
            auto tit = tetVertex.find(t);
            if (tit == tetVertex.end()) {
                Vector3f c = 0.25f * (V.col(T(0, t)) + V.col(T(1, t)) + V.col(T(2, t)) + V.col(T(3, t)));
                tit = tetVertex.emplace(t, (uint32_t) verts.size()).first;
                verts.push_back(c);
            }

            int g;
            if (rot1.index >= 4 && rot1.index <= 9) {
                // Original red/blue classification (detectOrientationSingularitiesTet).
                for (int i = 0; i < 3; ++i) {
                    int i0 = T(face[(i + 1) % 3], t), i1 = T(face[(i + 2) % 3], t);
                    rot2 *= AxisToAxisRot(Quaternion(Q.col(i0)), Quaternion(Q.col(i1)));
                }
                for (int i = 0; i < 3; ++i) {
                    int i0 = T(face[(i + 2) % 3], t), i1 = T(face[(i + 3) % 3], t);
                    rot3 *= AxisToAxisRot(Quaternion(Q.col(i0)), Quaternion(Q.col(i1)));
                }
                Quaternion q0 = Q.col(T(face[0], t)), q1 = Q.col(T(face[1], t)), q2 = Q.col(T(face[2], t));
                Vector3f normal = Vector3f(V.col(T(face[1], t)) - V.col(T(face[0], t))).cross(
                                  Vector3f(V.col(T(face[2], t)) - V.col(T(face[0], t)))).normalized();
                Vector3f a1 = q0.toMatrix().col(rot1.preservedAxis()) * (rot1.index > 6 ? 1 : -1);
                Vector3f a2 = q1.toMatrix().col(rot2.preservedAxis()) * (rot2.index > 6 ? 1 : -1);
                Vector3f a3 = q2.toMatrix().col(rot3.preservedAxis()) * (rot3.index > 6 ? 1 : -1);
                g = (a1 + a2 + a3).normalized().dot(normal) > 0 ? 0 : 1;
            } else if (rot1.index <= 15) {  // 1-3: about a frame axis, 10-15: about a face diagonal
                g = 2;
            } else {
                g = 3;
            }
            groups[g].emplace_back(fit->second, tit->second);
        }
    }

    st.faces = singularFaces.size();
    st.rot90_pos = groups[0].size();
    st.rot90_neg = groups[1].size();
    st.rot180 = groups[2].size();
    st.rot120 = groups[3].size();
    st.segments = st.rot90_pos + st.rot90_neg + st.rot180 + st.rot120;

    FILE *f = fopen(path.c_str(), "w");
    if (!f) return false;
    fprintf(f, "# Singularity graph of the RoSy frame field (segments face centroid -> tet centroid)\n");
    fprintf(f, "# singular faces: %zu, segments: %zu\n", st.faces, st.segments);
    for (const Vector3f &v : verts)
        fprintf(f, "v %.9g %.9g %.9g\n", v.x(), v.y(), v.z());
    const char *names[4] = { "rot90_pos", "rot90_neg", "rot180", "rot120" };
    for (int g = 0; g < 4; ++g) {
        if (groups[g].empty()) continue;
        fprintf(f, "g %s\n", names[g]);
        for (auto &s : groups[g])
            fprintf(f, "l %u %u\n", s.first + 1, s.second + 1);
    }
    bool ok = ferror(f) == 0;
    ok &= fclose(f) == 0;
    return ok;
}

// Medit .mesh of the TetGen output as stored by tet_meshing(): vertices (same indices and
// float coordinates as in _field.txt), boundary triangles (oriented, mF) and tetrahedra
// (TetGen's vertex order, identical to what TetGen's own -g output writes). 1-based.
bool write_tet_medit(const std::string &path, const MultiResolutionHierarchy &mRes) {
    const MatrixXf &V = mRes.mV[0];
    const MatrixXu &F = mRes.mF;
    const MatrixXu &T = mRes.mT;
    FILE *f = fopen(path.c_str(), "w");
    if (!f) return false;
    fprintf(f, "MeshVersionFormatted 1\nDimension 3\nVertices\n%lld\n", (long long) V.cols());
    for (int i = 0; i < V.cols(); ++i)
        fprintf(f, "%.9g %.9g %.9g 0\n", V(0, i), V(1, i), V(2, i));
    fprintf(f, "Triangles\n%lld\n", (long long) F.cols());
    for (int i = 0; i < F.cols(); ++i)
        fprintf(f, "%u %u %u 1\n", F(0, i) + 1, F(1, i) + 1, F(2, i) + 1);
    fprintf(f, "Tetrahedra\n%lld\n", (long long) T.cols());
    for (int t = 0; t < T.cols(); ++t)
        fprintf(f, "%u %u %u %u 0\n", T(0, t) + 1, T(1, t) + 1, T(2, t) + 1, T(3, t) + 1);
    fprintf(f, "End\n");
    bool ok = ferror(f) == 0;
    ok &= fclose(f) == 0;
    return ok;
}

} // namespace

int main(int argc, char **argv) {
    std::string input, prefix;
    double tetLenFactor = -1, tetLen = -1;
    int iters = 200, threads = -1;
    bool hierarchy = true, saveTet = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> const char * {
            if (i + 1 >= argc) { fprintf(stderr, "Missing value for %s\n", a.c_str()); exit(BAD_ARGS); }
            return argv[++i];
        };
        if (a == "-i" || a == "--input") input = next();
        else if (a == "-o" || a == "--output") prefix = next();
        else if (a == "--tet-len-factor") tetLenFactor = atof(next());
        else if (a == "--tet-len") tetLen = atof(next());
        else if (a == "--iters") iters = atoi(next());
        else if (a == "--threads" || a == "-t") threads = atoi(next());
        else if (a == "--no-hierarchy") hierarchy = false;
        else if (a == "--save-tet") saveTet = true;
        else if (a == "-h" || a == "--help") { usage(argv[0]); return OK; }
        else { fprintf(stderr, "Unknown argument: %s\n", a.c_str()); usage(argv[0]); return BAD_ARGS; }
    }
    if (input.empty() || prefix.empty() || iters < 1) { usage(argv[0]); return BAD_ARGS; }

    int nThreads = threads > 0 ? threads : (int) std::thread::hardware_concurrency();
#ifdef RHDM_ONETBB
    tbb::global_control tbbControl(tbb::global_control::max_allowed_parallelism, nThreads);
#else
    tbb::task_scheduler_init tbbInit(nThreads);
#endif

    Json js;
    js.str("input", input);
    js.integer("threads", nThreads);
    js.integer("iters_per_level", iters);
    js.raw("hierarchy", hierarchy ? "true" : "false");
    const std::string statsPath = prefix + "_stats.json";

    auto fail = [&](int code, const std::string &stage, const std::string &msg) {
        fprintf(stderr, "ERROR [%s]: %s\n", stage.c_str(), msg.c_str());
        js.str("status", "failed");
        js.str("failed_stage", stage);
        js.str("error", msg);
        js.write(statsPath);
        return code;
    };

    MultiResolutionHierarchy mRes;

    // ---- load
    {   // MultiResolutionHierarchy::load() segfaults on an OBJ without vertices/faces: pre-check.
        std::ifstream is(input);
        if (!is) return fail(LOAD_FAILED, "load", "cannot open " + input);
        size_t nv = 0, nf = 0;
        std::string line;
        while (std::getline(is, line)) {
            if (line.compare(0, 2, "v ") == 0) ++nv;
            else if (line.compare(0, 2, "f ") == 0) ++nf;
        }
        if (nv < 4 || nf < 4)
            return fail(LOAD_FAILED, "load", "OBJ has " + std::to_string(nv) + " vertices and " +
                        std::to_string(nf) + " faces; a closed surface needs at least 4 of each");
    }
    auto t0 = Clock::now();
    bool loaded = false;
    try { loaded = mRes.load(input); }
    catch (const std::exception &e) { return fail(LOAD_FAILED, "load", e.what()); }
    if (!loaded || mRes.mF.cols() == 0 || mRes.mV[0].cols() < 4)
        return fail(LOAD_FAILED, "load", "could not read a triangle mesh from the OBJ file");
    double tLoad = seconds_since(t0);
    js.integer("input_vertices", mRes.mV[0].cols());
    js.integer("input_triangles", mRes.mF.cols());
    js.num("input_avg_edge_length", mRes.ms.mAverageEdgeLength);
    js.num("bbox_diagonal", (mRes.mAABB.max - mRes.mAABB.min).norm());

    // ---- tet meshing
    // load() sets tet_elen = 0.3 * 3.5 * avg_edge = 1.05 * avg_edge (tElen_ratio = 1).
    // The original batch mode never changes it afterwards, so that is the default here.
    if (tetLen > 0) mRes.tet_elen = (Float) tetLen;
    else if (tetLenFactor > 0) mRes.tet_elen = (Float) (tetLenFactor * mRes.ms.mAverageEdgeLength);
    js.num("tet_target_edge_length", mRes.tet_elen);

    t0 = Clock::now();
    try {
        if (!mRes.tet_meshing())
            return fail(TET_FAILED, "tet_meshing", "tet_meshing() returned false");
    } catch (int code) {
        return fail(TET_FAILED, "tet_meshing", "TetGen terminated with error code " + std::to_string(code) +
                    " (1: out of memory, 2: internal error, 3: self-intersecting input, "
                    "4: very small input feature, 5: very close input facets, 10: input error)");
    } catch (const std::exception &e) {
        return fail(TET_FAILED, "tet_meshing", e.what());
    }
    double tTet = seconds_since(t0);
    if (mRes.mT.cols() == 0)
        return fail(TET_FAILED, "tet_meshing", "TetGen produced no tetrahedra");
    js.integer("tet_vertices", mRes.mV[0].cols());
    js.integer("tets", mRes.mT.cols());
    js.integer("boundary_triangles", mRes.mF.cols());

    // Written as soon as the tet mesh exists (so it survives a later RoSy failure); not timed.
    if (saveTet) {
        const std::string tetPath = prefix + "_tet.mesh";
        auto tw = Clock::now();
        if (!write_tet_medit(tetPath, mRes))
            return fail(WRITE_FAILED, "write", "could not write " + tetPath);
        js.str("tet_mesh_file", tetPath);
        js.num("time_tet_export_s", seconds_since(tw));
    }

    // ---- hierarchy + RoSy
    double tBuild = 0, tRosy = 0;
    try {
        t0 = Clock::now();
        mRes.build();
        tBuild = seconds_since(t0);
        js.integer("hierarchy_levels", mRes.levels());

        // Optimizer::run() with mOptimizeOrientations = true, alignment = randomization = true.
        t0 = Clock::now();
        int level = hierarchy ? std::max(0, mRes.levels() - 2) : 0;
        while (true) {
            for (int it = 0; it < iters; ++it)
                mRes.smoothOrientationsTet((uint32_t) level, true, true);
            if (level == 0)
                break;
            --level;
            mRes.prolongOrientations(level);
        }
        tRosy = seconds_since(t0);
    } catch (const std::exception &e) {
        return fail(FIELD_FAILED, "rosy", e.what());
    }

    js.num("time_load_s", tLoad);
    js.num("time_tet_meshing_s", tTet);
    js.num("time_hierarchy_build_s", tBuild);
    js.num("time_rosy_solve_s", tRosy);
    js.num("time_rosy_total_s", tBuild + tRosy);
    js.num("time_field_total_s", tTet + tBuild + tRosy);

    // ---- outputs (not timed as part of the field)
    SingularityStats ss;
    t0 = Clock::now();
    if (!write_field(prefix + "_field.txt", mRes))
        return fail(WRITE_FAILED, "write", "could not write " + prefix + "_field.txt");
    if (!write_singularities(prefix + "_singularities.obj", mRes, ss))
        return fail(WRITE_FAILED, "write", "could not write " + prefix + "_singularities.obj");
    js.num("time_export_s", seconds_since(t0));

    js.integer("singular_faces", ss.faces);
    js.integer("singular_segments", ss.segments);
    js.integer("singular_segments_rot90_pos", ss.rot90_pos);
    js.integer("singular_segments_rot90_neg", ss.rot90_neg);
    js.integer("singular_segments_rot180", ss.rot180);
    js.integer("singular_segments_rot120", ss.rot120);
    js.str("status", "ok");
    if (!js.write(statsPath)) {
        fprintf(stderr, "ERROR [write]: could not write %s\n", statsPath.c_str());
        return WRITE_FAILED;
    }

    printf("\n[rosy_field] %s\n", input.c_str());
    printf("  tet mesh: %lld vertices, %lld tets (target edge %.4g)\n",
           (long long) mRes.mV[0].cols(), (long long) mRes.mT.cols(), mRes.tet_elen);
    printf("  time: tet %.3fs | hierarchy %.3fs | rosy %.3fs | field total %.3fs\n",
           tTet, tBuild, tRosy, tTet + tBuild + tRosy);
    printf("  singular faces: %zu (segments %zu)\n", ss.faces, ss.segments);
    return OK;
}
