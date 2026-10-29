"""Parallel gain runs: launch workers, show progress, merge, summarise."""
import csv
import glob
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time

from . import field
from .util import binary, ncpu

TORR_PER_BAR = 750.0616827


class WorkerError(RuntimeError):
    pass


def gas_args(cfg, penning=True):
    g = cfg["gas"]
    comps = ",".join(f"{k}:{v}" for k, v in g["components"].items())
    args = ["--gas", comps,
            "--pressure-torr", f"{g['pressure_bar'] * TORR_PER_BAR:.6f}",
            "--temperature-k", f"{g['temperature_K']}"]
    p = g["penning"]
    r = p["r"] if penning else None
    if r is None:
        args += ["--penning-r", "off"]
    else:
        args += ["--penning-r", "auto" if r == "auto" else repr(float(r))]
        args += ["--penning-lambda", repr(float(p["lambda"] or 0.0))]
        if p["gas"] and r != "auto":
            args += ["--penning-gas", p["gas"]]
    return args


def derive_seed(*parts):
    """Deterministic, well-separated 31-bit seed from arbitrary parts."""
    h = hashlib.sha1(":".join(str(p) for p in parts).encode()).hexdigest()
    return int(h[:8], 16) % 2147483646 + 1


def run_key(cfg, h, seed):
    """Identity of a gain run apart from its length: field, gas, settings, seed."""
    # settings that do not change the physics (length, seed handled separately,
    # watchdog) must not invalidate finished results
    sim = {k: v for k, v in cfg["simulation"].items()
           if k not in ("events", "seed", "stall_timeout_s")}
    blob = json.dumps({"field": h, "gas": cfg["gas"], "sim": sim, "seed": seed},
                      sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]


# ----------------------------------------------------------------- progress
def _render(console, statuses, first):
    if not first:
        console.write(f"\x1b[{len(statuses)}A")
    for label, cur, total, state in statuses:
        frac = cur / total if total else 0.0
        fill = int(round(22 * frac))
        bar = "#" * fill + "-" * (22 - fill)
        tag = state if state else f"{int(round(100 * frac)):3d}%"
        console.write("\r\x1b[K" + f" {label:>4} [{bar}] {cur:>6}/{total:<6} {tag}\n")
    console.flush()


STALL_TIMEOUT_S = 120     # no collision for this long -> the process is hung
MAX_STALLS_PER_JOB = 10   # give up after jobs x this many hangs in one run


def _read_int(path):
    try:
        with open(path) as f:
            return int(f.read().split()[0])
    except Exception:
        return None


def _count_rows(path):
    try:
        with open(path) as f:
            return max(sum(1 for _ in f) - 1, 0)
    except OSError:
        return 0


def _trim_partial_line(path):
    """A process killed mid-write can leave half a line: drop it."""
    with open(path, "rb") as f:
        data = f.read()
    if data and not data.endswith(b"\n"):
        with open(path, "wb") as f:
            f.write(data[:data.rfind(b"\n") + 1])


class _Chunk:
    """One worker process simulating n events into its own part file."""

    def __init__(self, idx, n, outdir, label, cmd, seed):
        self.idx, self.n = idx, n
        stem = os.path.join(outdir, f".{label}_c{idx:04d}")
        self.part, self.prog, self.heart = stem + ".csv", stem + ".progress", stem + ".heart"
        self.log = os.path.join(outdir, "logs", f"{label}_c{idx:04d}.log")
        with open(self.log, "w") as lf:
            self.p = subprocess.Popen(cmd + ["--events", str(n), "--seed", str(seed), "--out", self.part,
                                             "--progress", self.prog, "--heartbeat", self.heart],
                                      stdout=lf, stderr=subprocess.STDOUT)
        self.last_sig, self.last_change = None, time.time()

    def done(self):
        return _read_int(self.prog) or 0

    def stalled(self, timeout):
        sig = (_read_int(self.heart), _read_int(self.prog))
        if sig != self.last_sig:
            self.last_sig, self.last_change = sig, time.time()
            return False
        return time.time() - self.last_change > timeout

    def finish(self):
        """Rows actually written (complete lines); removes the helper files."""
        _trim_partial_line(self.part)
        for f in (self.prog, self.heart):
            if os.path.exists(f):
                os.remove(f)
        return _count_rows(self.part)


def _chunk_size(pending, jobs):
    """Guided scheduling: large chunks first, small ones at the end, so that no
    slot is left idle while another works through a long tail."""
    return max(1, min(pending, math.ceil(pending / (2 * jobs))))


def _run_pool(pending, jobs, start_idx, outdir, label, cmd, seeds, stall_timeout, console):
    """Simulate `pending` events on `jobs` slots. Returns (part files written,
    number of hung processes replaced)."""
    tty = getattr(console, "isatty", lambda: False)()
    total = pending
    done_total = 0
    active, parts = [], []
    idx, stalls = start_idx, 0
    first, last_print = True, 0.0
    try:
        while pending > 0 or active:
            while len(active) < jobs and pending > 0:
                n = _chunk_size(pending, jobs)
                active.append(_Chunk(idx, n, outdir, label, cmd, seeds(idx)))
                parts.append(active[-1].part)
                idx += 1
                pending -= n
            for c in list(active):
                rc = c.p.poll()
                if rc is None:
                    if c.stalled(stall_timeout):
                        c.p.kill()
                        c.p.wait()
                        rows = c.finish()
                        stalls += 1
                        if stalls > MAX_STALLS_PER_JOB * jobs:
                            raise WorkerError(f"{stalls} processes hung inside Garfield; giving up "
                                              f"(last log {c.log})")
                        msg = (f"  process c{c.idx}: no collision for {stall_timeout:.0f} s (hung inside "
                               f"Garfield); its {c.n - rows} remaining event(s) go back to the queue")
                        if tty:
                            console.write("\r\x1b[K" + msg + "\n")
                            first = True
                        else:
                            print(msg)
                        done_total += rows
                        pending += c.n - rows
                        active.remove(c)
                    continue
                if rc != 0:
                    with open(c.log) as f:
                        tail = "".join(f.readlines()[-8:])
                    raise WorkerError(f"process c{c.idx} exited with {rc} (log {c.log}):\n{tail}")
                rows = c.finish()
                done_total += rows
                pending += c.n - rows            # should be 0
                active.remove(c)
            current = done_total + sum(c.done() for c in active)
            if tty:
                lines = [(f"c{c.idx}", c.done(), c.n, "") for c in active]
                lines += [("", 0, 0, "idle")] * (jobs - len(active))
                lines.append(("all", current, total, ""))
                _render(console, lines, first)
                first = False
            elif time.time() - last_print > 60 or not (pending or active):
                print(f"  progress: {current}/{total} events")
                last_print = time.time()
            if pending or active:
                time.sleep(0.5)
    except BaseException:
        for c in active:
            if c.p.poll() is None:
                c.p.kill()
        raise
    return parts, stalls


# ----------------------------------------------------------------- stats
def summarise(rows):
    eff = [int(r["effective_gain"]) for r in rows]
    ab = [int(r["absolute_gain"]) for r in rows]
    n = len(eff)
    if n == 0:
        return {"events": 0}

    def ms(v):
        m = sum(v) / n
        sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
        return m, sd

    me, se = ms(eff)
    ma, sa = ms(ab)
    nz = [v for v in eff if v > 0]
    return {
        "events": n,
        "effective_gain": me, "effective_gain_err": se / math.sqrt(n), "effective_gain_rms": se,
        "absolute_gain": ma, "absolute_gain_err": sa / math.sqrt(n), "absolute_gain_rms": sa,
        "zero_fraction": 1 - len(nz) / n,
        "effective_gain_nonzero": (sum(nz) / len(nz)) if nz else 0.0,
        "truncated_events": sum(int(r.get("truncated", 0)) for r in rows),
    }


def format_summary(s, title="Stack gain"):
    if not s.get("events"):
        return f"{title}: no events"
    lines = [
        "=" * 64,
        f"{title}",
        f"  effective gain : {s['effective_gain']:.4g} +/- {s['effective_gain_err']:.2g}"
        f"   (RMS {s['effective_gain_rms']:.3g})",
        f"  absolute gain  : {s['absolute_gain']:.4g} +/- {s['absolute_gain_err']:.2g}"
        f"   (RMS {s['absolute_gain_rms']:.3g})",
        f"  events         : {s['events']}   zero-signal fraction {100 * s['zero_fraction']:.1f}%",
    ]
    if s.get("truncated_events"):
        lines.append(f"  WARNING: {s['truncated_events']} event(s) hit the avalanche size limit "
                     "(gain underestimated; raise simulation.avalanche_size_limit)")
    if s.get("stalled_events"):
        lines.append(f"  NOTE: {s['stalled_events']} event(s) hung inside Garfield and were replaced "
                     "by new ones (worker restarted with a new seed)")
    lines.append("=" * 64)
    return "\n".join(lines)


def read_events(path):
    with open(path) as f:
        return list(csv.DictReader(f))


# ----------------------------------------------------------------- the run
def run_gain(cfg, stack, h, outdir, events, jobs=None, seed=None, label="gain"):
    """Run `events` avalanches with `jobs` workers on field map `h`; write
    <outdir>/<label>_events.csv and <label>_summary.json; return the summary."""
    if not field.is_ready(h):
        raise RuntimeError(f"field map {h} is not built/valid")
    fdir = field.field_dir(h)
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    from . import penning           # imports this module
    run_cfg = penning.effective(cfg)
    sim = cfg["simulation"]
    seed = sim["seed"] if seed is None else seed
    jobs = max(1, min(jobs or ncpu(), events))
    exe = binary("simgem_gain")
    common = ["--field", fdir] + gas_args(run_cfg) + [
        "--z-start", f"{stack.z_start / 10:.8f}",
        "--z-readout", f"{stack.z_readout / 10:.8f}",
        "--size-limit", str(sim["avalanche_size_limit"]),
        "--primary-energy", str(sim["primary_energy_eV"])]

    # Every event a process finishes is in a part file .<label>_*.csv. If the
    # run is interrupted, a rerun of the SAME job (run key, events) keeps them
    # and simulates only the rest, whatever the number of workers. Anything else
    # in the folder is stale and removed.
    manifest = os.path.join(outdir, f".{label}_manifest.json")
    key = run_key(cfg, h, seed)
    resume = False
    if os.path.exists(manifest):
        try:
            with open(manifest) as f:
                m = json.load(f)
            resume = m.get("run_key") == key and m.get("events") == events
        except (OSError, ValueError):
            resume = False
    old_parts = sorted(glob.glob(os.path.join(outdir, f".{label}_*.csv")))
    if not resume:
        for f in glob.glob(os.path.join(outdir, f".{label}_*")):
            if not f.endswith("_manifest.json"):
                os.remove(f)
        old_parts = []
        with open(manifest, "w") as f:
            json.dump({"run_key": key, "events": events}, f)
    for f in glob.glob(os.path.join(outdir, f".{label}_*.progress")) + \
            glob.glob(os.path.join(outdir, f".{label}_*.heart")):
        os.remove(f)
    for f in old_parts:
        _trim_partial_line(f)
    kept = min(sum(_count_rows(f) for f in old_parts), events)
    start_idx = 1 + max([int(m.group(1)) for m in (re.search(r"_c(\d+)\.csv$", f) for f in old_parts) if m],
                        default=-1)
    stall = sim.get("stall_timeout_s") or STALL_TIMEOUT_S
    print(f">>> {events} events on {jobs} workers, field {h}, -> {outdir}/"
          + (f"  (resuming: {kept} events kept from an interrupted run)" if kept else ""))
    t0 = time.time()
    parts, stalls = _run_pool(events - kept, jobs, start_idx, outdir, label, [exe] + common,
                              lambda i: derive_seed(seed, outdir, label, "chunk", i), stall, sys.__stdout__)

    rows = []
    merged = os.path.join(outdir, f"{label}_events.csv")
    with open(merged, "w", newline="") as out:
        wr = None
        for k, part in enumerate(old_parts + parts):
            for r in read_events(part):
                if len(rows) == events:          # a resumed run may hold one extra
                    break
                r = {"worker": k, **r, "event": len(rows)}
                if wr is None:
                    wr = csv.DictWriter(out, fieldnames=list(r))
                    wr.writeheader()
                wr.writerow(r)
                rows.append(r)
    if len(rows) != events:
        raise WorkerError(f"expected {events} events, workers returned {len(rows)}")
    for f in old_parts + parts:
        os.remove(f)
    os.remove(manifest)
    s = summarise(rows)
    s.update({"field_hash": h, "seed": seed, "run_key": run_key(cfg, h, seed),
              "jobs": jobs, "processes": len(parts),
              "stalled_events": stalls,
              "resumed_events": kept,
              "wall_time_s": round(time.time() - t0, 1),
              "penning_r": cfg["gas"]["penning"]["r"],
              "penning_r_micro": run_cfg["gas"]["penning"]["r"] if run_cfg is not cfg else None})
    with open(os.path.join(outdir, f"{label}_summary.json"), "w") as f:
        json.dump(s, f, indent=2)
    print(f"  merged -> {merged}")
    return s
