"""
Penning transfer matched to the table formula (gas.penning.match = "table").

Published transfer probabilities r (Sahin et al., JINST 5 (2010) P05002 and
NIM A 768 (2014) 104; Garfield's r = "auto") are defined through the Townsend
coefficient of the table formula alpha (1 + r sum nu_exc / sum nu_ion), fitted
to measured gain curves. The microscopic avalanche, with the same r, gives a
lower alpha (1-6 % in Ar:CO2 / Ar:CH4), worth tens of percent in a GEM stack.
With matching on, `simgem_penning` finds the r_micro for which microscopic
tracking reproduces the table alpha at r over a set of fields, and the gain
runs use r_micro. Results are cached in gas/penning/ (SIMGEM_PENNING) per gas,
Penning settings and fields, so a voltage scan calibrates once. The
calibrations committed in the package's own gas/penning/ are found from any
working directory (read only, e.g. inside the container).
"""
import copy
import hashlib
import json
import os

from . import config as C
from .runner import gas_args
from .util import ROOT, binary, run_streamed

CACHE_DIR = os.path.abspath(os.environ.get("SIMGEM_PENNING", os.path.join("gas", "penning")))
SHIPPED_DIR = os.path.join(ROOT, "gas", "penning")
MATCH_VERSION = 1
DEFAULT_FIELDS_KV_CM = (30.0, 45.0, 60.0, 80.0)   # at 1 atm; scaled with pressure


def match_fields(cfg):
    """Fields (kV/cm) at which r_micro is calibrated."""
    p = cfg["gas"]["penning"]
    if p["match_fields_kV_cm"]:
        return [float(v) for v in p["match_fields_kV_cm"]]
    scale = cfg["gas"]["pressure_bar"] / 1.01325
    return [round(v * scale, 1) for v in DEFAULT_FIELDS_KV_CM]


def calibration_key(cfg):
    p = cfg["gas"]["penning"]
    blob = json.dumps({"v": MATCH_VERSION, "gas": C.gas_key(cfg),
                       "r": p["r"], "lambda": p["lambda"], "penning_gas": p["gas"],
                       "fields": match_fields(cfg), "events": p["match_events"]},
                      sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]


def calibrate(cfg):
    """The cached or freshly computed calibration (dict, see src/penning.cpp)."""
    key = calibration_key(cfg)
    path = os.path.join(CACHE_DIR, f"{key}.json")
    for cand in (path, os.path.join(SHIPPED_DIR, f"{key}.json")):
        if os.path.exists(cand):
            with open(cand) as f:
                return json.load(f)
    os.makedirs(CACHE_DIR, exist_ok=True)
    tmp = os.path.join(CACHE_DIR, f".{key}-{os.getpid()}.json")
    log = os.path.join(CACHE_DIR, f"{key}.log")
    p = cfg["gas"]["penning"]
    print(f">>> Penning matching for {C.gas_label(cfg)}, r = {p['r']} "
          f"(fields {match_fields(cfg)} kV/cm; a few minutes, once per gas)")
    run_streamed([binary("simgem_penning")] + gas_args(cfg) + [
        "--fields", ",".join(f"{1e3 * v:g}" for v in match_fields(cfg)),
        "--events", str(p["match_events"]), "--out", tmp],
        log_path=log, keep=lambda line: line.startswith(("Target", "  r =", "RESULT", "ERROR")))
    with open(tmp) as f:
        res = json.load(f)
    res["key"] = key
    res["gas"] = C.gas_label(cfg)
    with open(tmp, "w") as f:
        json.dump(res, f, indent=2)
    os.replace(tmp, path)
    if res.get("clamped"):
        print("  WARNING: the matching r_micro lies outside [0, 1] and was clamped")
    return res


def effective(cfg):
    """The config the avalanche runs with: unchanged without matching, else a
    copy whose Penning r is r_micro (on the component the table r applies to)."""
    p = cfg["gas"]["penning"]
    if not p.get("match"):
        return cfg
    res = calibrate(cfg)
    out = copy.deepcopy(cfg)
    q = out["gas"]["penning"]
    q["r_table"] = res["r_table"]
    q["r"] = res["r_micro"]
    q["lambda"] = res["lambda"]
    if res["penning_gas"]:
        q["gas"] = res["penning_gas"]
    return out
