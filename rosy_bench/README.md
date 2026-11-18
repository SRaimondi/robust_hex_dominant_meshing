# rosy_bench — headless RoSy frame-field extraction

Runs only the orientation-field part of Gao et al. 2017 (load OBJ → TetGen tet mesh →
multi-resolution hierarchy → quaternion RoSy smoothing) and exports the result.
The GUI, position field and hex-dominant extraction are not built.

## Build

Requires CMake ≥ 3.16, a C++14 compiler, oneTBB (`brew install tbb` / `apt install libtbb-dev`)
and the submodules (`git submodule update --init --recursive`; only Eigen and pcg32 are used).

```sh
cmake -S rosy_bench -B build/rosy-release -DCMAKE_BUILD_TYPE=Release
cmake --build build/rosy-release -j
```

## Run

Single mesh:

```sh
build/rosy-release/rosy_field -i mesh.obj -o out/mesh [--tet-len-factor 1.05] [--tet-len L]
                               [--iters 200] [--threads N] [--no-hierarchy] [--save-tet]
```

Folder of meshes (each mesh in its own process; failures, crashes and timeouts are recorded
and the batch continues):

```sh
python3 rosy_bench/run_batch.py <input_dir> <output_dir> [--timeout 3600] [--skip-existing]
                                [--recursive] [--tet-len-factor k] [--iters n] [--threads n]
```

Writes `<output_dir>/<name>/{<name>_field.txt, <name>_singularities.obj, <name>_stats.json, log.txt}`
and `<output_dir>/summary.csv` (one row per mesh: status, failed stage, error, sizes, timings).

## Defaults (match the original code)

* Target tet edge length = 1.05 × average input edge length: this is what `load()` sets and
  the original batch mode never changes (its `-tl`/`-s` flags have no effect on tet sizing).
* 200 smoothing iterations per hierarchy level (hard-coded for tet meshes in
  `Optimizer::setOptimizeOrientations`), boundary alignment and neighbour randomisation on.
* `Float` = `float` (`SINGLE_PRECISION`), as in the original build. Results are deterministic
  and independent of the thread count.

## Timings (`_stats.json`, `summary.csv`)

| key | what |
|---|---|
| `time_tet_meshing_s` | `tet_meshing()`: background mesh + TetGen `pqm` + boundary extraction |
| `time_hierarchy_build_s` | `build()`: normals, adjacency, multi-resolution hierarchy, random init |
| `time_rosy_solve_s` | all `smoothOrientationsTet` + `prolongOrientations` calls |
| `time_field_total_s` | sum of the three above |

Loading the OBJ and writing outputs are reported separately (`time_load_s`, `time_export_s`).
The BVH that the original `build()` constructs is used only by the GUI, so it is compiled out
(`RHDM_HEADLESS`) and is not counted.

## `_field.txt` format (text, 0-based indices)

```
# comment lines
FRAMEFIELD 1
QUATERNION xyzw
VERTICES <nV>
x y z qx qy qz qw b        # unit quaternion (Hamilton); frame axes = columns of R(q); b = 1 on the boundary
EDGES <nE>
i j                        # i < j, undirected neighbourhood graph (= tet edges)
TETS <nT>                  # may be 0 for producers without a tet mesh
a b c d
```

## `_singularities.obj`

Same test as `detectOrientationSingularitiesTet()`: a tet face is singular when composing the
matched cube symmetries around its edges is not the identity. Each (singular face, incident tet)
pair becomes a line segment from the face centroid to the tet centroid. Shared endpoints are
merged, so the segments form a connected graph. OBJ groups give the holonomy class:
`rot90_pos` / `rot90_neg` (±90° about a frame axis, signed with the original code's red/blue
rule), `rot180`, and `rot120`.
