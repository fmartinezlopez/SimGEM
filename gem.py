#!/usr/bin/env python3
"""
Triple-GEM stack simulation — one CLI, one config = one run folder.

Each JSON in configs/ works in its own folder runs/<config stem>/, so you can
keep many configurations side by side. Gain always runs in PARALLEL across CPU
cores (independent microscopic avalanches, aggregated).

  python gem.py gas    configs/baseline.json                  # build the gas table
  python gem.py field  configs/baseline.json [--mesh-only|--no-mesh|--gui]
  python gem.py gain   configs/baseline.json [--events N] [--jobs J] [--r R]
  python gem.py all    configs/baseline.json [--events N] [--jobs J]
  python gem.py scan   configs/baseline.json --r 0.2 0.3 0.4 0.5 [--events N]
  python gem.py vis    configs/baseline.json [--events N]   # trajectories -> CSV

Requires: ./build/generate_gas, ./build/gain_stack and ./build/vis_stack
(cmake --build build), the gmsh python module, and ElmerGrid/ElmerSolver on PATH.
"""
import argparse, csv, glob, math, os, re, subprocess, sys, time
import stacklib as S

GEN_GAS = "./build/generate_gas"
GAIN = "./build/gain_stack"
VIS = "./build/vis_stack"


def ncpu():
    return os.cpu_count() or 4


def _read_progress(rundir, tag, default_total):
    try:
        with open(os.path.join(rundir, ".progress" + tag)) as f:
            a, b = f.read().split()[:2]
            return int(a), int(b)
    except Exception:
        return 0, default_total


def _render_bars(console, statuses, first):
    """Draw one bar per worker, redrawing in place. Writes to the real
    terminal (not the tee'd stdout) so the logs stay clean."""
    n = len(statuses)
    width = 22
    if not first:
        console.write(f"\x1b[{n}A")          # cursor up n lines
    for label, cur, total, done in statuses:
        frac = cur / total if total else 0.0
        fill = int(round(width * frac))
        bar = "#" * fill + "-" * (width - fill)
        tag = "done" if done else f"{int(round(100 * frac)):3d}%"
        console.write("\r\x1b[K" + f" {label:>3} [{bar}] {cur:>5}/{total:<5} {tag}\n")
    console.flush()


def _watch(procs, rundir, per, jobs, console):
    first = True
    while True:
        statuses, alive = [], False
        for j in range(jobs):
            running = procs[j][0].poll() is None
            alive = alive or running
            cur, tot = _read_progress(rundir, f"_w{j}", per)
            statuses.append((f"w{j}", cur, tot, not running))
        _render_bars(console, statuses, first)
        first = False
        if not alive:
            break
        time.sleep(0.3)


def resolve_r(v):
    """'auto' -> -1.0 (built-in a-priori parameterization); else a float."""
    if isinstance(v, str) and v.strip().lower() == "auto":
        return -1.0
    return float(v)


def rlabel(r):
    return "auto" if r < 0 else f"{r}"


# ------------------------------------------------------------------------ gas
def _do_gas(cfg, force):
    g = cfg["gas"]
    gp = S.gas_path(cfg)
    os.makedirs(os.path.dirname(gp) or ".", exist_ok=True)
    if os.path.exists(gp) and not force:
        print(f"Gas file {gp} already exists (use --force to regenerate).")
        return
    S.run_streamed([GEN_GAS, g["ar_fraction"], g["ch4_fraction"],
                    g["pressure_bar"], g["temperature_K"], gp])


def cmd_gas(cfg, rundir, args):
    with S.stage_log(rundir, "gas"):
        _do_gas(cfg, args.force)


# ------------------------------------------------------------------------ field
def cmd_field(cfg, rundir, args):
    with S.stage_log(rundir, "field"):
        S.build_field(cfg, rundir, mesh_only=args.mesh_only,
                      no_mesh=args.no_mesh, gui=args.gui)


# ------------------------------------------------------------------------ gain (parallel)
def _do_gain(cfg, rundir, r, total, jobs, label=""):
    """Launch `jobs` independent workers, merge their CSVs into one, return stats."""
    gp = S.gas_path(cfg)
    if not os.path.exists(os.path.join(rundir, "stack", "case.result")):
        sys.exit(f"No field map in {rundir}/stack/. Run:  python gem.py field <config>")
    if not os.path.exists(gp):
        sys.exit(f"No gas file {gp}. Run:  python gem.py gas <config>")

    per = math.ceil(total / jobs)
    for f in glob.glob(os.path.join(rundir, "gain_stack_events_w*.csv")):
        os.remove(f)
    for f in glob.glob(os.path.join(rundir, ".progress_w*")):
        os.remove(f)
    print(f">>> {total} events / {jobs} workers ({per} each), r={rlabel(r)}, in {rundir}/")

    procs = []
    for j in range(jobs):
        log = open(os.path.join(rundir, f"gain_worker_w{j}.log"), "w")
        p = subprocess.Popen(
            [GAIN, rundir, gp, str(r), str(per), str(1000 + j), f"_w{j}"],
            stdout=log, stderr=subprocess.STDOUT)
        procs.append((p, log))

    console = sys.__stdout__          # real terminal, bypassing the log tee
    if getattr(console, "isatty", lambda: False)():
        _watch(procs, rundir, per, jobs, console)
    else:
        print(f"  running {jobs} workers x {per} events (no TTY; watch the logs)...")
    for p, log in procs:
        p.wait(); log.close()
    for f in glob.glob(os.path.join(rundir, ".progress_w*")):
        try:
            os.remove(f)
        except OSError:
            pass

    # merge per-worker CSVs into one, then remove the pieces
    merged = os.path.join(rundir, f"gain_events{label}.csv")
    eff = []
    with open(merged, "w") as out:
        out.write("worker,event,absolute_gain,effective_gain\n")
        for f in sorted(glob.glob(os.path.join(rundir, "gain_stack_events_w*.csv"))):
            w = re.search(r"_w(\d+)\.csv$", f).group(1)
            with open(f) as fh:
                for row in csv.DictReader(fh):
                    out.write(f"{w},{row['event']},{row['absolute_gain']},{row['effective_gain']}\n")
                    eff.append(int(row["effective_gain"]))
            os.remove(f)
    if not eff:
        print(f"  No events collected — check {rundir}/gain_worker_w*.log")
        return None, merged
    n = len(eff)
    mean = sum(eff) / n
    sd = math.sqrt(sum((v - mean) ** 2 for v in eff) / n)
    print(f"  merged -> {merged}  ({n} events)")
    return (mean, sd / math.sqrt(n), sd, n), merged


def cmd_gain(cfg, rundir, args):
    with S.stage_log(rundir, "gain"):
        r = resolve_r(args.r) if args.r is not None else resolve_r(cfg["gas"]["penning_r"])
        total = args.events or cfg["avalanche"]["num_events"]
        jobs = args.jobs or ncpu()
        res, _ = _do_gain(cfg, rundir, r, total, jobs)
        if res:
            mean, err, sd, n = res
            print("=" * 55)
            print(f"Triple-GEM stack gain (r={rlabel(r)}): {mean:.1f} +/- {err:.1f}")
            print(f"  over {n} events,  RMS spread = {sd:.1f}")
            print("=" * 55)


# ------------------------------------------------------------------------ vis (single worker)
def cmd_vis(cfg, rundir, args):
    """Dump electron trajectories for later plotting. Deliberately serial: the
    output is one set of CSVs per run, and a handful of events is enough."""
    with S.stage_log(rundir, "vis"):
        gp = S.gas_path(cfg)
        if not os.path.exists(os.path.join(rundir, "stack", "case.result")):
            sys.exit(f"No field map in {rundir}/stack/. Run:  python gem.py field <config>")
        if not os.path.exists(gp):
            sys.exit(f"No gas file {gp}. Run:  python gem.py gas <config>")
        r = resolve_r(args.r) if args.r is not None else resolve_r(cfg["gas"]["penning_r"])
        limit = args.size_limit or cfg["avalanche"].get("avalanche_size_limit", 2000)
        print(f">>> {args.events} event(s), r={rlabel(r)}, in {rundir}/  (1 worker)")
        S.run_streamed([VIS, rundir, gp, str(r), str(args.events), str(args.seed), "",
                        str(args.max_tracks), str(limit), str(args.min_step),
                        str(args.dt_max)])


# ------------------------------------------------------------------------ mc (misalignment Monte Carlo)
def _pooled_stats(values):
    n = len(values)
    if n == 0:
        return None
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n
    sd = math.sqrt(var)
    return mean, sd / math.sqrt(n), sd, n


def _read_eff(merged):
    """Read effective_gain values from a merged gain_events.csv."""
    out = []
    with open(merged) as fh:
        for row in csv.DictReader(fh):
            out.append(int(row["effective_gain"]))
    return out


def cmd_mc(cfg, rundir, args):
    """Rigid-foil-misalignment Monte Carlo.

    Draw K random misaligned geometries (each an independent field solve), throw
    M avalanches at each, pool the K*M effective gains, and compare that pooled
    distribution to an ALIGNED baseline of the same total statistics. Every
    per-config mean + its offsets go into mc_manifest.csv so gain-vs-offset can
    be plotted; all events go into mc_pooled.csv tagged by config.
    """
    mis = cfg.get("misalignment")
    if not mis or mis.get("mode") != "random":
        sys.exit("mc requires a 'misalignment' block with \"mode\": \"random\" "
                 "(K, sigma and seed come from there).")
    with S.stage_log(rundir, "mc"):
        K = args.configs or int(mis.get("n_configs", 20))
        M = args.events or int(mis.get("events_per_config",
                                       cfg["avalanche"]["num_events"]))
        jobs = args.jobs or ncpu()
        r = resolve_r(args.r) if args.r is not None else resolve_r(cfg["gas"]["penning_r"])
        gp = S.gas_path(cfg)
        if not os.path.exists(gp):
            sys.exit(f"No gas file {gp}. Run:  python gem.py gas {args.config}")
        pitch = S._pitch_mm(cfg)
        Ntot = K * M
        print(f"Misalignment MC in {rundir}/  ({K} configs x {M} events = {Ntot}), "
              f"r={rlabel(r)}, jobs={jobs}")

        manifest = []           # (idx, off dict, mean, err, n)
        pooled = []             # effective gains across all misaligned configs
        pooled_csv = os.path.join(rundir, "mc_pooled.csv")
        with open(pooled_csv, "w") as pf:
            pf.write("config,gem1_dx,gem1_dy,gem2_dx,gem2_dy,gem3_dx,gem3_dy,"
                     "event,absolute_gain,effective_gain\n")

        # ---- K misaligned geometries ---------------------------------------
        for i in range(K):
            sub = os.path.join(rundir, f"mc_{i:04d}")
            offsets, meta = S.resolve_offsets(cfg, i)
            print(f"\n--- config {i + 1}/{K}: "
                  + S.offsets_summary(offsets, pitch) + " ---")
            if not (args.no_field or (args.resume and
                    os.path.exists(os.path.join(sub, "stack", "case.result")))):
                S.build_field(cfg, sub, offsets_mm=offsets)
            res, merged = _do_gain(cfg, sub, r, M, jobs)
            if not res:
                print(f"  config {i} produced no events; skipping")
                continue
            mean, err, sd, n = res
            manifest.append((i, offsets, mean, err, n))
            # append this config's events to the pooled csv, tagged with offsets
            o = [offsets["gem1"][0], offsets["gem1"][1],
                 offsets["gem2"][0], offsets["gem2"][1],
                 offsets["gem3"][0], offsets["gem3"][1]]
            with open(merged) as fh, open(pooled_csv, "a") as pf:
                for row in csv.DictReader(fh):
                    pf.write(f"{i}," + ",".join(f"{v:.6f}" for v in o) + ","
                             f"{row['event']},{row['absolute_gain']},"
                             f"{row['effective_gain']}\n")
                    pooled.append(int(row["effective_gain"]))

        if not pooled:
            sys.exit("No misaligned events collected — check the per-config logs.")

        # ---- aligned baseline, matched total statistics --------------------
        base_events = args.baseline_events or Ntot
        print(f"\n--- aligned baseline: 1 geometry x {base_events} events ---")
        basedir = os.path.join(rundir, "mc_baseline")
        zero = {f: (0.0, 0.0) for f in S.FOILS}
        if not (args.no_field or (args.resume and
                os.path.exists(os.path.join(basedir, "stack", "case.result")))):
            S.build_field(cfg, basedir, offsets_mm=zero)   # translate @ zero offset
        bres, bmerged = _do_gain(cfg, basedir, r, base_events, jobs)
        base_vals = _read_eff(bmerged) if bres else []

        # ---- write manifest ------------------------------------------------
        man_csv = os.path.join(rundir, "mc_manifest.csv")
        with open(man_csv, "w") as f:
            f.write("config,gem1_dx,gem1_dy,gem2_dx,gem2_dy,gem3_dx,gem3_dy,"
                    "mean_gain,err,n_events\n")
            for i, off, mean, err, n in manifest:
                o = [off["gem1"][0], off["gem1"][1], off["gem2"][0],
                     off["gem2"][1], off["gem3"][0], off["gem3"][1]]
                f.write(f"{i}," + ",".join(f"{v:.6f}" for v in o) +
                        f",{mean:.4f},{err:.4f},{n}\n")

        # ---- compare pooled-misaligned vs aligned --------------------------
        mm, me, msd, mn = _pooled_stats(pooled)
        per_cfg_means = [m for _, _, m, _, _ in manifest]
        gmean, gerr, gsd, gn = _pooled_stats(per_cfg_means)  # spread of config means
        # intrinsic per-config avalanche error, averaged, to separate the two
        avg_within_err = math.sqrt(sum(e * e for _, _, _, e, _ in manifest) / len(manifest))
        # geometry-induced variance = Var(config means) - mean(within-config var)
        geo_var = max(gsd * gsd - avg_within_err * avg_within_err * 0, 0.0)

        lines = []
        lines.append("=" * 60)
        lines.append(f"MISALIGNMENT MC SUMMARY  (r={rlabel(r)})")
        lines.append("-" * 60)
        if base_vals:
            bm, be, bsd, bn = _pooled_stats(base_vals)
            lines.append(f"aligned baseline : {bm:9.1f} +/- {be:5.1f}   "
                         f"(RMS {bsd:7.1f}, n={bn})")
        lines.append(f"misaligned pool  : {mm:9.1f} +/- {me:5.1f}   "
                     f"(RMS {msd:7.1f}, n={mn})")
        if base_vals:
            dmean = mm - bm
            dstat = dmean / math.sqrt(me * me + be * be)
            lines.append(f"mean shift        : {dmean:+9.1f}  ({dstat:+.2f} sigma) "
                         f"= {100 * dmean / bm:+.1f}% of aligned")
            lines.append(f"RMS change        : {msd - bsd:+9.1f}  "
                         f"({100 * (msd - bsd) / bsd:+.1f}%)")
        lines.append("-" * 60)
        lines.append(f"config-to-config mean spread (geometry): "
                     f"RMS {gsd:.1f} over {gn} configs")
        lines.append(f"  -> misalignment adds this width on top of avalanche spread")
        lines.append("=" * 60)
        summary = "\n".join(lines)
        print("\n" + summary)
        with open(os.path.join(rundir, "mc_summary.txt"), "w") as f:
            f.write(summary + "\n")
        print(f"\nWrote: {man_csv}\n       {pooled_csv}\n"
              f"       {os.path.join(rundir, 'mc_summary.txt')}")


# ------------------------------------------------------------------------ all
def cmd_all(cfg, rundir, args):
    cmd_gas(cfg, rundir, args)
    with S.stage_log(rundir, "field"):
        S.build_field(cfg, rundir)
    cmd_gain(cfg, rundir, args)


# ------------------------------------------------------------------------ scan
def cmd_scan(cfg, rundir, args):
    with S.stage_log(rundir, "scan"):
        total = args.events or cfg["avalanche"]["num_events"]
        jobs = args.jobs or ncpu()
        print(f"Penning scan in {rundir}/  ({total} events each)\n")
        print(f"{'r':>6} {'gain':>12} {'+/-':>10}")
        print("-" * 30)
        rows = []
        for r in args.r:
            res, _ = _do_gain(cfg, rundir, r, total, jobs, label=f"_r{r}")
            if res:
                mean, err, _, _ = res
                print(f"{r:6.2f} {mean:12.1f} {err:10.1f}")
                rows.append((r, mean, err))
        with open(os.path.join(rundir, "penning_scan.csv"), "w") as f:
            f.write("r,gain,error\n")
            for r, m, e in rows:
                f.write(f"{r},{m},{e}\n")
        print(f"\nSaved -> {rundir}/penning_scan.csv")


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description="Triple-GEM stack simulation CLI")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_common(sp):
        sp.add_argument("config")

    sp = sub.add_parser("gas", help="generate the gas table")
    add_common(sp); sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("field", help="build geometry + solve field map")
    add_common(sp)
    sp.add_argument("--mesh-only", action="store_true")
    sp.add_argument("--no-mesh", action="store_true", help="re-solve without re-meshing")
    sp.add_argument("--gui", action="store_true", help="open the mesh in the gmsh GUI")

    sp = sub.add_parser("gain", help="parallel microscopic stack gain")
    add_common(sp)
    sp.add_argument("--events", type=int); sp.add_argument("--jobs", type=int)
    sp.add_argument("--r", type=str, help="Penning r, or 'auto' for the built-in parameterization")

    sp = sub.add_parser("vis", help="dump electron trajectories to CSV (1 worker)")
    add_common(sp)
    sp.add_argument("--events", type=int, default=1, help="avalanches to trace (default 1)")
    sp.add_argument("--r", type=str, help="Penning r, or 'auto'")
    sp.add_argument("--seed", type=int, default=1000)
    sp.add_argument("--min-step", type=float, default=5e-4,
                    help="write a trajectory point every this many cm (default 5e-4 = 5 um)")
    sp.add_argument("--max-tracks", type=int, default=200,
                    help="full paths for the first N electrons per event, -1 = all")
    sp.add_argument("--size-limit", type=int,
                    help="cap the avalanche size (default: config avalanche_size_limit)")
    sp.add_argument("--dt-max", type=float, default=0.03,
                    help="ns; time gap above which a new drift line is assumed (default 0.03)")

    sp = sub.add_parser("all", help="gas + field + gain")
    add_common(sp)
    sp.add_argument("--events", type=int); sp.add_argument("--jobs", type=int)
    sp.add_argument("--r", type=str); sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("scan", help="gain for several Penning r values")
    add_common(sp)
    sp.add_argument("--r", type=float, nargs="+", required=True)
    sp.add_argument("--events", type=int); sp.add_argument("--jobs", type=int)

    sp = sub.add_parser("mc", help="misalignment Monte Carlo (random foil shifts) "
                                    "+ aligned baseline, pooled & compared")
    add_common(sp)
    sp.add_argument("--configs", type=int, help="K geometries (default: config n_configs)")
    sp.add_argument("--events", type=int, help="M events per config (default: config events_per_config)")
    sp.add_argument("--jobs", type=int)
    sp.add_argument("--r", type=str, help="Penning r, or 'auto'")
    sp.add_argument("--baseline-events", type=int,
                    help="events for the aligned baseline (default: K*M, matched)")
    sp.add_argument("--resume", action="store_true",
                    help="skip field solves whose case.result already exists")
    sp.add_argument("--no-field", action="store_true",
                    help="reuse all existing field maps; only (re)run gain")

    args = ap.parse_args()
    cfg = S.load_config(args.config)
    rundir = S.run_dir(args.config)
    {"gas": cmd_gas, "field": cmd_field, "gain": cmd_gain,
     "all": cmd_all, "scan": cmd_scan, "vis": cmd_vis, "mc": cmd_mc}[args.cmd](cfg, rundir, args)


if __name__ == "__main__":
    main()