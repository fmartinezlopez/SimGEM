#!/usr/bin/env python3
"""
Export SimGEM scan results as GORG reference gain curves.

GORG (gorgana) overlays every CSV in GORG/metadata/reference_curves/ on its
gain-curve figures; the columns it reads are label, gas, gas_ratio,
pressure_psi, V_gem, gain. Extra columns (errors, Penning r, ...) are kept for
reference and ignored by the plots.

    python scripts/export_gorg_curves.py runs/scans/gorg/gain_curve_*psi \\
        --out ../GORG/metadata/reference_curves/simgem.csv

Every finished point of the given scan folders becomes one row. The gain is the
fitted gain (zero-inflated Polya mean, see simgem/analysis.py) when the fit
succeeded, otherwise the mean effective gain over all events; both are the
quantity a Fe-55 photopeak measures. Points must come from configs with an
"hv" block (V_gem as in the lab).
"""
import argparse
import csv
import datetime
import glob
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from simgem.config import BAR_PER_PSI as PSI, GAS_DISPLAY  # noqa: E402


def gas_strings(components):
    names = [GAS_DISPLAY.get(k.lower(), k) for k in components]
    tot = sum(components.values())
    ratio = "/".join(f"{100 * v / tot:g}" for v in components.values())
    return "/".join(names), ratio


def point_rows(scan_dir, label):
    rows = []
    for pdir in sorted(d for d in glob.glob(os.path.join(scan_dir, "p[0-9][0-9][0-9]_*")) if os.path.isdir(d)):
        try:
            with open(os.path.join(pdir, "config.json")) as f:
                cfg = json.load(f)["resolved"]
            with open(os.path.join(pdir, "gain_summary.json")) as f:
                summ = json.load(f)
        except FileNotFoundError:
            continue                          # point not run (yet) or failed
        if not cfg.get("hv"):
            sys.exit(f"{pdir}: config has no 'hv' block, so V_gem is not defined")
        fit = {}
        fp = os.path.join(pdir, "gain_fit.json")
        if os.path.exists(fp):
            with open(fp) as f:
                fit = json.load(f)
        use_fit = fit.get("ok") and "gain" in fit
        gas, ratio = gas_strings(cfg["gas"]["components"])
        r = cfg["gas"]["penning"]["r"]
        rows.append({
            "label": label.format(r=r if r is not None else "off"),
            "gas": gas, "gas_ratio": ratio,
            "pressure_psi": round(cfg["gas"]["pressure_bar"] / PSI, 1),
            "V_gem": cfg["hv"]["V_gem"],
            "gain": fit["gain"] if use_fit else summ["effective_gain"],
            "gain_err": fit["gain_err"] if use_fit else summ["effective_gain_err"],
            "gain_method": "zero-inflated Polya fit" if use_fit else "mean over all events",
            "mean_all": summ["effective_gain"], "mean_all_err": summ["effective_gain_err"],
            "penning_r": r, "events": summ["events"],
            "truncated_events": summ.get("truncated_events", 0),
            "point": os.path.relpath(pdir),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan_dirs", nargs="+", help="scan output folders (runs/scans/...)")
    ap.add_argument("--label", default="SimGEM (r = {r})",
                    help="curve label; {r} is replaced by the Penning r (default: %(default)s)")
    ap.add_argument("--out", help="CSV to write (default: print to stdout)")
    args = ap.parse_args()

    rows = []
    for d in args.scan_dirs:
        if not os.path.isdir(d):
            sys.exit(f"no such scan folder: {d}")
        rows += point_rows(d, args.label)
    if not rows:
        sys.exit("no finished points found")
    rows.sort(key=lambda r: (r["label"], r["pressure_psi"], r["V_gem"]))
    bad = [r for r in rows if r["truncated_events"]]
    if bad:
        print(f"WARNING: {len(bad)} point(s) have events that hit the avalanche size limit "
              "(gain underestimated); raise simulation.avalanche_size_limit", file=sys.stderr)

    out = open(args.out, "w", newline="") if args.out else sys.stdout
    out.write(f"# SimGEM gain curves exported {datetime.datetime.now().isoformat(timespec='seconds')} "
              f"from: {' '.join(args.scan_dirs)}\n")
    w = csv.DictWriter(out, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
    if args.out:
        out.close()
        print(f"{len(rows)} points -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
