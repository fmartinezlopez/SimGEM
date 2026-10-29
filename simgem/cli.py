"""Command-line interface (gem.py)."""
import argparse
import json
import os
import re
import shutil
import sys

from . import config as C
from . import analysis, corrections, field, misalignment, penning, runner, sweep, transport
from .util import binary, run_streamed, stage_log

# Outputs live under the working directory (override with environment variables).
RUNS_DIR = os.path.abspath(os.environ.get("SIMGEM_RUNS", "runs"))
GAS_DIR = os.path.abspath(os.environ.get("SIMGEM_GAS", "gas"))


def rundir_for(args):
    if args.out:
        return args.out
    d = os.path.join(RUNS_DIR, C.run_name(args.config))
    if args.set:
        slug = "_".join(re.sub(r"[^A-Za-z0-9.=+-]+", "-", s) for s in args.set)
        d = os.path.join(d, "set_" + slug[:120])
    return d


def snapshot(cfg, rundir, args, extra=None):
    os.makedirs(rundir, exist_ok=True)
    d = {"config_file": os.path.abspath(args.config), "overrides": args.set or [],
         "resolved": cfg}
    d.update(extra or {})
    with open(os.path.join(rundir, "config.json"), "w") as f:
        json.dump(d, f, indent=2, default=str)


# ----------------------------------------------------------------- commands
def cmd_check(cfg, rundir, args):
    st = C.Stack(cfg)
    h = st.field_hash()
    print(f"Gas      : {C.gas_label(cfg)}   Penning r = {cfg['gas']['penning']['r']}"
          + (f" on {cfg['gas']['penning']['gas']}" if cfg["gas"]["penning"]["gas"] else ""))
    print(f"Cell     : {st.cell[0] * 1e3:.2f} x {st.cell[1] * 1e3:.2f} um (translational-periodic)")
    print("Foils (top -> bottom):")
    for f in st.foils:
        h_ = f["hole"]
        off = st.offsets[f["name"]]
        print(f"  {f['name']:<6} {f['pattern']:<6} pitch {f['pitch_um']:g} um, kapton {f['kapton_thickness_um']:g} um, "
              f"Cu {f['copper_thickness_um']:g} um, hole {h_['profile']} {h_['outer_diameter_um']:g}/"
              f"{h_['inner_diameter_um']:g} um, dV = {f['dV_V']:g} V, "
              f"offset ({off[0] * 1e3:+.1f}, {off[1] * 1e3:+.1f}) um")
    if cfg.get("hv"):
        hv = cfg["hv"]
        vd = hv.get("V_drift_gap", (hv.get("V_drift") or 0) - hv["V_gem"])
        print(f"HV       : V_gem = {hv['V_gem']:g} V (divider), V_drift = {hv['V_gem'] + vd:g} V "
              f"({vd:g} V over the {hv['drift_length_mm']:g} mm drift gap)")
        print("           divider: " + ", ".join(f"{k} {100 * v:.2f}%" for k, v in hv["shares"].items()))
    print("Gaps (top -> bottom):")
    for g in cfg["gaps"]:
        print(f"  {g['name']:<10} {g['thickness_mm']:6.3f} mm  {g['field_kV_cm']:7.3f} kV/cm")
    print("Electrodes:")
    for n in st.electrodes():
        print(f"  {n:<10} z = {st.z[n]:8.4f} mm   V = {st.V[n]:10.2f} V")
    print(f"Primaries start at z = {st.z_start:.4f} mm; readout counts z < {st.z_readout * 1e3:g} um")
    top = st.z[st.foils[0]["name"] + "_top"]
    if st.z_start - top < 2 * max(st.cell):
        print(f"  WARNING: primaries start only {(st.z_start - top) * 1e3:.0f} um above the top foil "
              "(inside its field distortion); increase the drift gap or reduce start_depth_um")
    print(f"Field map: {h} ({'cached' if field.is_ready(h) else 'not built yet'})")
    if cfg.get("misalignment_mc"):
        mc = cfg["misalignment_mc"]
        print(f"Misalignment MC: {mc['geometries']} geometries x {mc['events_per_geometry']} events, "
              f"{mc['distribution']}" + (f" sigma {mc['sigma_um']} um" if mc["distribution"] == "gaussian" else ""))


def cmd_field(cfg, rundir, args):
    st = C.Stack(cfg)
    with stage_log(rundir, "field"):
        if args.mesh_only or args.gui:
            d = os.path.join(rundir, "mesh_preview")
            field.preview_mesh(st, d, gui=args.gui)
            print(f"Mesh written to {d}/stack.msh (not solved, not cached)")
            return
        h = field.build(st, force=args.force, vtu=args.vtu, keep_msh=args.keep_msh)
        snapshot(cfg, rundir, args, {"field_hash": h})
        print(f"Field map: {field.field_dir(h)}")


def _gain(cfg, rundir, args):
    st = C.Stack(cfg)
    h = st.field_hash()
    if not field.is_ready(h):
        if args.cmd == "gain":
            sys.exit(f"No field map for this geometry yet ({h}). Run:  python gem.py field {args.config}"
                     f"{' ' + ' '.join('--set ' + s for s in args.set) if args.set else ''}\n"
                     "(or use `run`, which builds it)")
        h = field.build(st)
    snapshot(cfg, rundir, args, {"field_hash": h})
    n = args.events or cfg["simulation"]["events"]
    s = runner.run_gain(cfg, st, h, rundir, n, jobs=args.jobs, seed=args.seed)
    print(runner.format_summary(s, f"Stack gain ({C.gas_label(cfg)}, r = {cfg['gas']['penning']['r']})"))
    k = corrections.factor(cfg)
    if k != 1.0:
        c = cfg["corrections"]
        print(f"  corrected effective gain: {k * s['effective_gain']:.4g} +/- {k * s['effective_gain_err']:.2g}"
              f"   (x{k:.3g}: charging-up {c['charging_up_per_foil']:g} and scale "
              f"{c['gain_scale_per_foil']:g} per foil, {len(cfg['foils'])} foil(s))")


def cmd_gain(cfg, rundir, args):
    with stage_log(rundir, "gain"):
        _gain(cfg, rundir, args)


def cmd_run(cfg, rundir, args):
    with stage_log(rundir, "run"):
        _gain(cfg, rundir, args)


def cmd_scan(args):
    """Scan file (scans/*.json), or CONFIG --param P --values ... (one axis)."""
    try:
        if sweep.is_scan_file(args.spec):
            if args.param:
                sys.exit("--param/--values are for a config; a scan file defines its own axes")
            spec = sweep.load_spec(args.spec)
        else:
            if not (args.param and args.values):
                sys.exit(f"{args.spec} is not a scan file (needs 'base' and 'axes'/'points'); "
                         "for a config give --param and --values")
            base_out = os.path.join(RUNS_DIR, C.run_name(args.spec))
            spec = sweep.adhoc_spec(args.spec, args.param, args.values, base_out)
        outdir = args.out or spec.get("output") or os.path.join(RUNS_DIR, "scans", spec["_name"])
        if args.dry_run:
            sweep.run(spec, outdir, args.set or [], events=args.events, jobs=args.jobs,
                      only=args.only, dry_run=True)
            return
        with stage_log(outdir, "scan"):
            rows = sweep.run(spec, outdir, args.set or [], events=args.events, jobs=args.jobs,
                             only=args.only, force=args.force, analyze_only=args.analyze_only,
                             stop_on_error=args.stop_on_error, plot=not args.no_plot)
    except (sweep.ScanError, C.ConfigError) as e:
        sys.exit(f"Scan error: {e}")
    if any(r["status"] in ("failed", "field-failed") for r in rows or []):
        sys.exit(3)                     # some points failed (the scan itself completed)


def cmd_analyze(args):
    try:
        path = analysis.find_events(args.target)
    except FileNotFoundError as e:
        sys.exit(str(e))
    analysis.analyze_run(path, column=args.column, threshold=args.threshold,
                         max_value=args.max, bin_width=args.bin_width, plot=not args.no_plot)


def cmd_mc(cfg, rundir, args):
    with stage_log(rundir, "mc"):
        snapshot(cfg, rundir, args)
        misalignment.run_mc(cfg, rundir, geometries=args.geometries, events=args.events,
                                reference_events=args.reference_events, jobs=args.jobs,
                                skip_transport=args.skip_transport, force=args.force)


def cmd_transport(cfg, rundir, args):
    with stage_log(rundir, "transport"):
        st = C.Stack(cfg)
        _, rep = transport.run(cfg, st, rundir, electrons=args.electrons,
                               distance_um=args.distance_um)
        print(transport.format_report(rep))
        with open(os.path.join(rundir, "transport.json"), "w") as f:
            json.dump(rep, f, indent=2)


def cmd_vis(cfg, rundir, args):
    st = C.Stack(cfg)
    h = st.field_hash()
    if not field.is_ready(h):
        sys.exit(f"No field map for this geometry yet ({h}). Run:  python gem.py field {args.config}")
    with stage_log(rundir, "vis"):
        limit = args.size_limit or cfg["simulation"]["avalanche_size_limit"]
        run_streamed([binary("simgem_vis"), "--field", field.field_dir(h)] + runner.gas_args(penning.effective(cfg)) + [
            "--events", args.events, "--seed", args.seed, "--outdir", rundir,
            "--max-tracks", args.max_tracks, "--size-limit", limit,
            "--min-step", args.min_step, "--dt-max", args.dt_max,
            "--z-start", f"{st.z_start / 10:.8f}", "--z-readout", f"{st.z_readout / 10:.8f}",
            "--primary-energy", cfg["simulation"]["primary_energy_eV"]])
        # add the stack layout for plotting
        p = os.path.join(rundir, "vis_meta.json")
        with open(p) as f:
            meta = json.load(f)
        fm = field.load_meta(h)
        meta.update({"field_hash": h, "electrodes_cm": {k: v / 10 for k, v in st.z.items()},
                     "kapton_cm": {k: [a / 10, b / 10] for k, (a, b) in st.kapton.items()},
                     "holes_cm": {k: [[x / 10, y / 10] for x, y in v] for k, v in fm["holes_mm"].items()},
                     "penning_r": cfg["gas"]["penning"]["r"], "gas": C.gas_label(cfg)})
        with open(p, "w") as f:
            json.dump(meta, f, indent=2)


def cmd_gastable(cfg, rundir, args):
    import hashlib
    key = dict(C.gas_key(cfg), table=cfg["gas"]["table"])
    h = hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()[:8]
    name = "_".join(f"{k}{v:g}" for k, v in cfg["gas"]["components"].items())
    out = os.path.join(GAS_DIR, f"{name}_{cfg['gas']['pressure_bar']:g}bar_{h}.gas")
    if os.path.exists(out) and not args.force:
        print(f"{out} already exists (--force to regenerate)")
        return
    os.makedirs(os.path.dirname(out), exist_ok=True)
    t = cfg["gas"]["table"]
    with stage_log(rundir, "gastable"):
        run_streamed([binary("simgem_gastable")] + runner.gas_args(cfg, penning=False) + [
            "--out", out, "--field-min", t["field_min_V_cm"], "--field-max", t["field_max_V_cm"],
            "--field-points", t["field_points"], "--max-energy", t["max_electron_energy_eV"],
            "--collisions", t["collisions_1e7"]])


def cmd_fields(args):
    rows = field.list_cache()
    if args.remove:
        for h in args.remove:
            d = field.field_dir(h)
            if os.path.isdir(d):
                shutil.rmtree(d)
                print(f"removed {d}")
            else:
                print(f"no field map {h}")
        return
    tot = 0
    print(f"Field-map cache: {field.FIELDS_DIR}")
    for h, meta, size in rows:
        tot += size
        offs = meta.get("offsets_mm", {})
        moved = [f"{n}({v[0] * 1e3:+.0f},{v[1] * 1e3:+.0f})" for n, v in offs.items()
                 if abs(v[0]) > 1e-6 or abs(v[1]) > 1e-6]
        nf = len(meta["key"]["foils"])
        print(f"  {h}  {meta['created']}  {size / 1e6:6.0f} MB  {nf} foils, "
              f"V_drift {meta['potentials_V']['drift']:.0f} V"
              + (f", offsets {' '.join(moved)}" if moved else ""))
    print(f"  total {tot / 1e6:.0f} MB in {len(rows)} map(s)")


# ----------------------------------------------------------------- parser
def build_parser():
    ap = argparse.ArgumentParser(
        prog="gem.py", description="Microscopic simulation of a full GEM stack.",
        epilog="Any config value can be overridden with --set, e.g.  "
               "--set gas.penning.r=0.3  --set gaps.1.thickness_mm=1.0  --set 'foils[*].dV_V=560' "
               "(list indices: gaps.1 or 'gaps[1]'; quote brackets in zsh)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("config", help="JSON config file")
        sp.add_argument("--set", action="append", metavar="PATH=VALUE",
                        help="override a config value (repeatable)")
        sp.add_argument("--out", help="output folder (default runs/<config name>[/set_...])")
        return sp

    def runopts(sp):
        sp.add_argument("--events", type=int, help="events (default: simulation.events)")
        sp.add_argument("--jobs", type=int, help="parallel workers (default: all cores)")
        sp.add_argument("--seed", type=int, help="base random seed (default: simulation.seed)")

    common(sub.add_parser("check", help="validate a config and show the derived stack"))

    sp = common(sub.add_parser("field", help="build + solve + validate the field map (cached)"))
    sp.add_argument("--force", action="store_true", help="rebuild even if cached")
    sp.add_argument("--vtu", action="store_true", help="also write case.vtu (+ E field) for ParaView")
    sp.add_argument("--keep-msh", action="store_true", help="keep the gmsh .msh file")
    sp.add_argument("--mesh-only", action="store_true", help="only mesh (preview, not cached)")
    sp.add_argument("--gui", action="store_true", help="mesh and open the gmsh GUI")

    runopts(common(sub.add_parser("gain", help="parallel microscopic stack gain")))
    runopts(common(sub.add_parser("run", help="field (if needed) + gain")))

    sp = sub.add_parser("scan", help="parameter scan from a scan file (scans/*.json)",
                        description="Run a scan file, or a quick one-axis scan of a config "
                                    "with --param/--values. Resumable; see README.")
    sp.add_argument("spec", help="scan file, or a config together with --param/--values")
    sp.add_argument("--param", help="(config only) path to vary, e.g. gas.penning.r")
    sp.add_argument("--values", nargs="+", help="(config only) values of --param")
    sp.add_argument("--set", action="append", metavar="PATH=VALUE",
                    help="extra override of the base config for every point")
    sp.add_argument("--out", help="output folder (default runs/scans/<scan name>)")
    sp.add_argument("--events", type=int, help="events per point (overrides the scan file)")
    sp.add_argument("--jobs", type=int, help="workers per point")
    sp.add_argument("--only", help="subset of points, e.g. 0,3-5")
    sp.add_argument("--dry-run", action="store_true", help="list points, field maps and work; run nothing")
    sp.add_argument("--force", action="store_true", help="rerun points that have results")
    sp.add_argument("--analyze-only", action="store_true", help="refit/tabulate existing results only")
    sp.add_argument("--stop-on-error", action="store_true", help="abort on the first failing point")
    sp.add_argument("--no-plot", action="store_true")

    sp = sub.add_parser("analyze", help="Polya fit of a run's gain distribution")
    sp.add_argument("target", help="gain_events.csv, a run folder, or a config (-> runs/<name>/)")
    sp.add_argument("--column", default="effective_gain", choices=["effective_gain", "absolute_gain"])
    sp.add_argument("--threshold", type=int,
                    help="fit only gains >= this with a plain Polya (e.g. to compare with a "
                         "single-electron spectrum); default: zero-inflated Polya of all events")
    sp.add_argument("--max", type=int, help="highest gain included (default: 99.9%% quantile)")
    sp.add_argument("--bin-width", type=int, help="histogram bin width (default: auto)")
    sp.add_argument("--no-plot", action="store_true")

    sp = common(sub.add_parser("mc", help="foil-misalignment Monte Carlo"))
    sp.add_argument("--geometries", type=int)
    sp.add_argument("--events", type=int, help="events per geometry")
    sp.add_argument("--reference-events", type=int)
    sp.add_argument("--jobs", type=int)
    sp.add_argument("--skip-transport", action="store_true", help="skip the sensitivity check")
    sp.add_argument("--force", action="store_true", help="rerun geometries that have results")

    sp = common(sub.add_parser("transport", help="drift/diffusion in the gaps + misalignment sensitivity"))
    sp.add_argument("--electrons", type=int, default=300)
    sp.add_argument("--distance-um", type=float, help="drift distance per field (default min(500 um, gap))")

    sp = common(sub.add_parser("vis", help="dump electron trajectories to CSV (1 worker)"))
    sp.add_argument("--events", type=int, default=1)
    sp.add_argument("--seed", type=int, default=1000)
    sp.add_argument("--min-step", type=float, default=5e-4, help="cm between written points")
    sp.add_argument("--max-tracks", type=int, default=200, help="full paths per event, -1 = all")
    sp.add_argument("--size-limit", type=int, help="avalanche size cap")
    sp.add_argument("--dt-max", type=float, default=0.03, help="ns; drift-line split threshold")

    sp = common(sub.add_parser("gastable", help="OPTIONAL Magboltz gas table (not used for the gain)"))
    sp.add_argument("--force", action="store_true")

    sp = sub.add_parser("fields", help="list (or --remove) cached field maps")
    sp.add_argument("--remove", nargs="+", metavar="HASH")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.cmd == "fields":
        return cmd_fields(args)
    if args.cmd == "scan":
        return cmd_scan(args)
    if args.cmd == "analyze":
        return cmd_analyze(args)
    try:
        cfg = C.load(args.config, args.set or [])
    except C.ConfigError as e:
        sys.exit(f"Config error in {args.config}: {e}")
    rundir = rundir_for(args)
    cmds = {"check": cmd_check, "field": cmd_field, "gain": cmd_gain, "run": cmd_run,
            "mc": cmd_mc, "transport": cmd_transport, "vis": cmd_vis,
            "gastable": cmd_gastable}
    try:
        cmds[args.cmd](cfg, rundir, args)
    except (runner.WorkerError, RuntimeError, C.ConfigError) as e:
        sys.exit(f"ERROR: {e}")
