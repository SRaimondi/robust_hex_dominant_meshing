#!/usr/bin/env python3
"""Run `rosy_field` on every OBJ file in a folder, isolating failures.

Each mesh runs in its own process, so a TetGen error, a crash (segfault) or a
hang (timeout) only affects that mesh; the batch continues and the failure is
reported in the summary.

Usage:
    python3 rosy_bench/run_batch.py <input_dir> <output_dir> [options]

Per mesh <name>.obj the outputs go to <output_dir>/<name>/:
    <name>_field.txt, <name>_singularities.obj, <name>_stats.json, log.txt
and one row is appended to <output_dir>/summary.csv.
"""

import argparse
import csv
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_BIN = Path(__file__).resolve().parent.parent / "build" / "rosy-release" / "rosy_field"

COLUMNS = [
    "mesh", "status", "failed_stage", "error", "return_code", "wall_time_s",
    "input_vertices", "input_triangles", "tet_target_edge_length",
    "tet_vertices", "tets", "hierarchy_levels", "threads", "iters_per_level",
    "time_tet_meshing_s", "time_hierarchy_build_s", "time_rosy_solve_s",
    "time_rosy_total_s", "time_field_total_s",
    "singular_faces", "singular_segments",
]


def tail(text, n=5):
    lines = [l for l in text.strip().splitlines() if l.strip()]
    return " | ".join(lines[-n:])


def run_one(binary, obj, out_dir, tool_args, timeout):
    name = obj.stem
    mesh_dir = out_dir / name
    mesh_dir.mkdir(parents=True, exist_ok=True)
    prefix = mesh_dir / name
    stats_path = Path(f"{prefix}_stats.json")
    if stats_path.exists():
        stats_path.unlink()  # never report a stale result

    cmd = [str(binary), "-i", str(obj), "-o", str(prefix)] + tool_args
    row = {"mesh": obj.name}
    t0 = time.monotonic()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, errors="replace", timeout=timeout or None)
        output, rc = proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as e:
        output = e.stdout.decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        rc = None
    row["wall_time_s"] = f"{time.monotonic() - t0:.3f}"
    (mesh_dir / "log.txt").write_text(" ".join(cmd) + "\n\n" + output)

    stats = {}
    if stats_path.exists():
        try:
            stats = json.loads(stats_path.read_text())
        except json.JSONDecodeError:
            pass
    for k in COLUMNS:
        if k in stats and k not in row:
            row[k] = stats[k]

    row["return_code"] = "" if rc is None else rc
    if rc is None:
        row["status"] = "timeout"
        row["error"] = f"killed after {timeout}s; last output: {tail(output)}"
    elif rc < 0:
        try:
            sig = signal.Signals(-rc).name
        except ValueError:
            sig = str(-rc)
        row["status"] = "crashed"
        row["error"] = f"terminated by {sig}; last output: {tail(output)}"
    elif rc != 0:
        row["status"] = "failed"
        row.setdefault("error", tail(output))
    else:
        row["status"] = stats.get("status", "failed")
        if row["status"] != "ok":
            row.setdefault("error", "exit code 0 but no valid stats file")
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir", type=Path)
    ap.add_argument("output_dir", type=Path)
    ap.add_argument("--bin", type=Path, default=DEFAULT_BIN, help=f"rosy_field executable (default: {DEFAULT_BIN})")
    ap.add_argument("--recursive", action="store_true", help="also search sub-folders for .obj files")
    ap.add_argument("--timeout", type=float, default=3600, help="seconds per mesh, 0 = no limit (default 3600)")
    ap.add_argument("--skip-existing", action="store_true",
                    help="skip meshes whose previous run in output_dir finished with status ok")
    # forwarded to rosy_field
    ap.add_argument("--tet-len-factor", type=float, help="tet edge length = k * avg input edge length (default 1.05)")
    ap.add_argument("--tet-len", type=float, help="absolute tet edge length")
    ap.add_argument("--iters", type=int, help="RoSy iterations per level (default 200)")
    ap.add_argument("--threads", type=int, help="TBB threads (default: all cores)")
    ap.add_argument("--no-hierarchy", action="store_true")
    ap.add_argument("--save-tet", action="store_true", help="also write <name>_tet.mesh")
    args = ap.parse_args()

    if not args.bin.is_file() or not os.access(args.bin, os.X_OK):
        sys.exit(f"rosy_field executable not found: {args.bin} (build it first, see rosy_bench/README.md)")
    pattern = "**/*" if args.recursive else "*"
    objs = sorted(p for p in args.input_dir.glob(pattern) if p.is_file() and p.suffix.lower() == ".obj")
    if not objs:
        sys.exit(f"no .obj files in {args.input_dir}")
    names = [p.stem for p in objs]
    dups = {n for n in names if names.count(n) > 1}
    if dups:
        sys.exit(f"duplicate mesh names would share an output folder: {sorted(dups)}")

    tool_args = []
    if args.tet_len_factor is not None: tool_args += ["--tet-len-factor", str(args.tet_len_factor)]
    if args.tet_len is not None: tool_args += ["--tet-len", str(args.tet_len)]
    if args.iters is not None: tool_args += ["--iters", str(args.iters)]
    if args.threads is not None: tool_args += ["--threads", str(args.threads)]
    if args.no_hierarchy: tool_args += ["--no-hierarchy"]
    if args.save_tet: tool_args += ["--save-tet"]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "summary.csv"

    # With --skip-existing, keep previous successful rows; otherwise start a fresh summary.
    previous = {}
    if args.skip_existing and summary_path.exists():
        with open(summary_path, newline="") as f:
            previous = {r["mesh"]: r for r in csv.DictReader(f) if r.get("status") == "ok"}

    rows = []
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for k, obj in enumerate(objs, 1):
            stats_file = args.output_dir / obj.stem / f"{obj.stem}_stats.json"
            if args.skip_existing and obj.name in previous and stats_file.exists():
                row = previous[obj.name]
                print(f"[{k}/{len(objs)}] {obj.name}: skipped (already ok)")
            else:
                print(f"[{k}/{len(objs)}] {obj.name} ...", end=" ", flush=True)
                row = run_one(args.bin, obj, args.output_dir, tool_args, args.timeout)
                if row["status"] == "ok":
                    print(f"ok  tets={row.get('tets')}  field={float(row['time_field_total_s']):.3f}s "
                          f"(tet {float(row['time_tet_meshing_s']):.3f} + hier {float(row['time_hierarchy_build_s']):.3f} "
                          f"+ rosy {float(row['time_rosy_solve_s']):.3f})  sing.faces={row.get('singular_faces')}")
                else:
                    print(f"{row['status'].upper()} [{row.get('failed_stage', '?')}] {row.get('error', '')}")
            rows.append(row)
            writer.writerow(row)
            f.flush()

    failed = [r for r in rows if r["status"] != "ok"]
    print(f"\n{len(rows) - len(failed)}/{len(rows)} succeeded. Summary: {summary_path}")
    if failed:
        print("Failures:")
        for r in failed:
            print(f"  {r['mesh']}: {r['status']} [{r.get('failed_stage') or '-'}] {r.get('error', '')}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
