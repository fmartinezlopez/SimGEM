# SimGEM

Microscopic simulation of a full GEM stack — any number of foils, any hole pattern and profile, any gas mixture — in one geometry and one field solve, with the gain computed by parallel microscopic avalanches (Garfield++) across multiple CPU cores.

```
config.json ─► gem.py ─► fields/<hash>/   field map (gmsh mesh + Elmer solve), cached & validated
                     └─► runs/<config>/   results: gain_events.csv, gain_summary.json, logs
```

Requirements: Garfield++ (with ROOT), Elmer (`ElmerGrid`, `ElmerSolver`) and the gmsh Python module. See the [installation guide](docs/install.md) (Linux, macOS, container). The [validation notes](docs/validation.md) document how far the outputs can be trusted.

```bash
cmake -S . -B build && cmake --build build -j      # builds the simgem_* programs
python gem.py check configs/baseline.json          # validate the config, print the stack
python gem.py run   configs/baseline.json          # field map (if needed) + gain
python -m unittest discover tests                  # fast checks of the Python side
```

## Commands

| command | what it does |
|---|---|
| `check` | validate a config; print electrode heights, potentials, cell, field-map hash |
| `field` | mesh + solve + validate the field map; cached by geometry hash (`--force`, `--vtu`, `--gui`, `--mesh-only`) |
| `gain`  | parallel microscopic stack gain on an existing field map |
| `run`   | `field` if needed, then `gain` |
| `scan`  | parameter scan from a scan file (`scans/*.json`), or a quick one-axis scan of a config |
| `analyze` | Polya fit of a run's gain distribution (`gain_fit.json`, `gain_fit.png`) |
| `mc`    | foil-misalignment Monte Carlo (needs a `misalignment_mc` block) |
| `transport` | drift velocity and diffusion in every gap, plus how much hole-pattern memory survives each transfer gap |
| `vis`   | electron trajectories of a few avalanches → CSV for plotting |
| `gastable` | *optional* Magboltz gas table (not used by the gain simulation) |
| `fields` | list (or `--remove`) cached field maps |

Every command takes `--set PATH=VALUE` (repeatable) to override any config value without editing the file, and `--out DIR`:

```bash
python gem.py run configs/baseline.json --set 'foils[*].dV_V=570' --set gas.penning.r=0.3
python gem.py run configs/baseline.json --set 'gas.components={"ar":95,"co2":5}'
python gem.py scan configs/baseline.json --param foils.1.dV_V --values 530 550 570
```

List elements are addressed as `gaps.1.thickness_mm` or `'gaps[1].thickness_mm'` (quote brackets: zsh treats them as a glob); `[*]` means every element. Results with overrides go to `runs/<config>/set_<overrides>/` unless `--out` is given.

## Configuration

Foils and gaps are listed top (drift side) → bottom (readout), there is one gap more than foils (drift, transfer 1, ..., induction). Units are in the key names. Defaults (in `simgem/config.py`) fill anything left out.

```jsonc
{
  "gas": {
    "components": {"ar": 98, "ch4": 2},      // up to 6, any Magboltz gas; fractions are normalised
    "pressure_bar": 8.0, "temperature_K": 293.15,
    "penning": {"r": 0.18, "lambda": 0.0, "gas": "ar",   // r: number, "auto" (binary mixtures), or null = off
                "match": null}                           // "table": r_micro matched to the table formula (see Gas)
  },
  "corrections": null,                       // null, "literature" or {"charging_up_per_foil": .., "gain_scale_per_foil": ..}
  "foil": {                                  // defaults for every foil
    "pattern": "hex",                        // "hex" (standard GEM) or "square"
    "pitch_um": 140, "kapton_thickness_um": 50, "copper_thickness_um": 5,
    "kapton_permittivity": 3.5,
    "offset_um": [0, 0],                     // lateral shift of the hole pattern
    "hole": {"profile": "double_conical",    // or "conical" (+ "wide_side": "top"/"bottom"), "cylindrical"
             "outer_diameter_um": 70, "inner_diameter_um": 50,
             "copper_diameter_um": null}     // null = same as the kapton hole (set larger for a rim)
  },
  "foils": [                                 // per-foil values override "foil"
    {"name": "gem1", "dV_V": 550},
    {"name": "gem2", "dV_V": 550, "offset_um": [70, 40.4]},
    {"name": "gem3", "dV_V": 550, "hole": {"inner_diameter_um": 55}}
  ],
  "gaps": [
    {"name": "drift",     "thickness_mm": 1.0, "field_kV_cm": 0.278},  // modeled slice of the drift volume
    {"name": "transfer1", "thickness_mm": 2.1, "field_kV_cm": 2.857},
    {"name": "transfer2", "thickness_mm": 2.1, "field_kV_cm": 2.381},
    {"name": "induction", "thickness_mm": 2.1, "field_kV_cm": 3.619}
  ],
  "mesh": {"fine_um": 6, "coarse_um": 60},   // optional: fine_distance_um, coarse_distance_um
  "simulation": {
    "events": 1000, "avalanche_size_limit": 500000, "seed": 1,
    "start_depth_um": null,                  // primaries start this far below the drift plane (default: mid drift gap)
    "readout_threshold_um": 10, "primary_energy_eV": 0.1
  },
  "misalignment_mc": { ... }                 // optional, see below
}
```

Gap thicknesses are measured between the outer copper surfaces. `gem.py` turns gap fields and foil voltages into electrode potentials (readout at 0 V, potentials decreasing upward). Foils may have different pitches/patterns as long as they are commensurate (e.g. 140 and 280 µm): the cell becomes a common supercell.

**Powering through a divider (`hv`).** A real detector is set by supply voltages, not per-foil values. With an `hv` block the foil voltages and gap fields are derived (and must not be given by hand):

```jsonc
"hv": {
  "V_gem": 2600,               // across the divider: GEM1 top sits at -V_gem, readout at 0
  "V_drift_gap": 500,          // V_drift - V_gem (or "V_drift": the absolute drift voltage)
  "drift_length_mm": 18.0,     // the TRUE drift gap: drift field = V_drift_gap / length
  "divider": {"gem1": 0.109, "transfer1": 0.2383, "gem2": 0.094,     // share of V_gem per
              "transfer2": 0.2383, "gem3": 0.082, "induction": 0.2384} // foil/gap; sums to 1
  // or "divider_resistance_kohm": {"gem1": 550, "transfer1": 1000, ...}  (shares = R / sum R)
}
```

The divider is taken as unloaded (the GEM electrodes draw no current), which is what the detector sees; a voltmeter on the chain reads lower.

A gain curve is then a scan of the single path `hv.V_gem`. `notes` (a list of strings) can record assumptions and sources in any config.

Example configs:
- `gorg`: the GORG setup, see below.
- `baseline`: hex foils, 8 bar Ar:CH₄.
- `baseline_square`: the original SimGEM geometry, for comparison.
- `misalign_gem2_worst`: a fixed offset.
- `misalign_mc`: the misalignment Monte Carlo.
- `ar_co2_thin_gaps`: Ar:CO₂ 50:50 with 0.5 mm transfer gaps, a case where
  alignment can matter.
- `validation/`: the benchmark set-ups of the [validation notes](docs/validation.md).

## Comparing with GORG

`configs/gorg.json` is the GORG setup: a TECHTRA (TTA) triple GEM in Ar:CH₄ 92:8 at 1 atm, with an 18 mm drift (1 mm simulated, true field applied), 2.1 mm transfer and induction gaps, and V_drift = V_GEM + 500 V. Foil names follow the TECHTRA board: GEM3 is the top foil (at −V_GEM, next to the drift), GEM1 the bottom one.

V_GEM is divided by the TECHTRA "GEM HV Divider Separate Drift cover V1.0", given as its resistances. The config reproduces that model's electrode potentials to 0.01 V:

| element | resistors | share of V_GEM |
|---|---|---|
| GEM3 (top) | 5 × 110 kΩ | 12.25 % |
| transfer32 | 5 × 200 kΩ | 22.27 % |
| GEM2 | 5 × 100 kΩ | 11.14 % |
| transfer21 | 5 × 200 kΩ | 22.27 % |
| GEM1 (bottom) | 4 × 84.5 kΩ + 102 kΩ | 9.80 % |
| induction | 5 × 200 kΩ | 22.27 % |

The foils are standard CERN GEM foils (140 µm hex pitch, 50 µm kapton, 5 µm Cu, 70/50 µm double-conical holes), and GORG pressures are absolute psi (14.5 psi is ambient). Penning transfer uses Garfield's Ar:CH₄ parameterisation (`"r": "auto"`, Şahin et al. 2010); the gas temperature is taken as 20 °C.

For the comparison with data, add Penning matching and the literature corrections ([validation notes](docs/validation.md), section 5):

```bash
python gem.py scan scans/gorg/gain_curve_14p5psi.json --set gas.penning.match=table \
       --set corrections=literature --out runs/scans/gorg/matched_14p5psi
python scripts/compare_gorg.py runs/scans/gorg/matched_14p5psi --corrected --adc-offset 18.5
```

```bash
export SIMGEM_JOBS=6 SIMGEM_PYTHON=~/miniconda3/envs/garfield/bin/python
scripts/run_scans.sh -n scans/gorg/*.json                  # what would run
scripts/run_scans.sh -d scans/gorg/gain_curve_14p5psi.json scans/gorg/gain_curve_15psi.json
python scripts/export_gorg_curves.py runs/scans/gorg/gain_curve_*psi \
       --out ../GORG/metadata/reference_curves/simgem.csv  # overlay on GORG's gain-curve figures
```

The exported gain is the zero-inflated-Polya gain, which is what GORG's photopeak-based gain measures. Cost grows with the gain (every avalanche electron is tracked): start with the 1 atm curves and use `--events 50` or `--only` for a first look at the high-gain points (up to ~10⁵ at 29.5 psi).

## Geometry and field

* The unit cell is a rectangle that is a translational period of every foil's hole pattern: p × p√3 for hex (one hole + four quarter holes), p × p for square. Side walls are translational-periodic (gmsh periodic mesh + Elmer periodic BC), so any lateral offset of any foil is represented exactly.
* Copper layers are real volumes removed from the domain; their surfaces are the electrodes, and electrons hitting them leave the field map (absorbed). `copper_thickness_um: 0` gives zero-thickness electrodes on the kapton.
* Field maps live in `fields/<hash>/`, where the hash covers exactly what determines the field (geometry, offsets, potentials, mesh). Configs, scan points and MC geometries that need the same field share it; a Penning or gas scan solves the field once.
* Every solve is validated before it can be used: Elmer's periodic projector must match all wall nodes, the potential must agree across opposite walls (to 10 mV; typically < 1 µV), and stay within the electrode potentials. A map that fails is kept in `fields/.failed-<hash>/` and nothing runs on it.

Look before you leap: `python gem.py field configs/x.json --gui` meshes and opens the gmsh GUI; `--vtu` also writes `case.vtu` for ParaView.

## Gain

Single primary electrons start at `start_depth_um` below the drift plane, uniformly over the cell, and are tracked microscopically through the whole stack. Per event (`gain_events.csv`):

* `absolute_gain`: electrons that were not attached, primary included (Garfield's avalanche size). It can be 0: in Ar:CH₄ at 8 bar about 2% of primaries are lost to CH₄ dissociative attachment inside the first hole before ionising.
* `effective_gain`: electrons that reach the readout (`< readout_threshold_um`); this folds in collection, multiplication, extraction and transfer losses of all foils. It is the quantity to compare with a measured gain.
* `n_ions`, `truncated` (the event hit `avalanche_size_limit`; reported as a warning since the gain is then underestimated)

`gain_summary.json` holds mean ± error, RMS, zero-signal fraction, the run key, the seed and the field hash. When comparing with a Fe-55 peak, use the mean over all events, not a Polya fit to non-zero events.

### Gain extraction (`gem.py analyze`)

```bash
python gem.py analyze runs/baseline                 # or a gain_events.csv, or a config
python gem.py analyze runs/baseline --threshold 3   # plain Polya above 3 electrons
```

The gain to compare with a measurement (readout current, or the position of a Fe-55 peak) is the mean effective gain over all events, zeros included: every primary electron counts, also those that give nothing at the readout.

By default the distribution is fitted with a **zero-inflated Polya**, using all events and no threshold: with probability 1 − ε a primary gives nothing (attached, not collected, all electrons lost), otherwise its gain follows a Polya(G0, θ). The likelihood factorises into the zero fraction and a Polya fit of the
events ≥ 1, done as an exact binned maximum-likelihood fit. The reported gain is the model's mean, ≈ ε·G0, the fitted counterpart of the mean over all events. G0 alone is the avalanche mean of the primaries that survive, not the gain. A single Polya fitted to all events including the zeros is biased (10–25% high in tests) and fits badly, because the zeros are a separate population.

`--threshold n` fits a plain Polya to the events ≥ n only, for comparing with single-electron spectra measured above an electronics threshold (convert the threshold to electrons). Output: G0 ± error, θ ± error, their correlation, χ²/ndof with p-value, ε (or the fraction above threshold), and the mean over all events. `gain_fit.json` and `gain_fit.png` are written next to the events.

## Parameter scans

A scan file (in `scans/`) names a base config and what to vary; each point is the base config plus overrides, with its own folder and full results.

```jsonc
{
  "description": "CH4 fraction x pressure, gap fields scaled with pressure",
  "base": "../configs/baseline.json",       // relative to the scan file
  "mode": "grid",                           // "grid": all combinations; "zip": paired values
  "axes": [
    {"name": "ch4", "path": "gas.components",
     "values": [{"ar": 98, "ch4": 2}, {"ar": 95, "ch4": 5}], "labels": ["2", "5"]},
    {"name": "p", "path": "gas.pressure_bar", "range": {"start": 4, "stop": 8, "step": 2}},
    {"name": "split", "set": {"foils[0].dV_V": [560, 540], "foils[2].dV_V": [520, 540]}}
  ],
  "derived": {"gaps[1].field_kV_cm": "base('gaps[1].field_kV_cm') * p / 8"},
  "set": {"simulation.avalanche_size_limit": 2000000},    // fixed, all points
  "points": [{"label": "nominal", "set": {}}],            // extra points, taken literally
  "events": 300,
  "analysis": {"polya": true, "column": "effective_gain"}   // "threshold": n -> plain Polya above n
}
```

* An axis is a `path` with `values` or a `range` (`start/stop/step`, or `start/stop/num`, `"log": true` for geometric spacing), or a `set` of several paths that move together (one axis). `labels` name the values in tables and folder names.
* `derived` values are computed per point from axis values (their names), the base config (`base('path')`) and math functions: e.g. keep reduced fields constant while scanning pressure, or keep the total stack voltage fixed.
* Any config path works: voltages, gap fields/thicknesses, gas composition, pressure, temperature, Penning r, pattern, pitch, hole diameters/profile, copper/kapton thickness, per-foil offsets, mesh, simulation settings.

```bash
python gem.py scan scans/gas_pressure.json --dry-run   # points, field maps needed, work to do
python gem.py scan scans/gas_pressure.json             # run (resumable)
python gem.py scan scans/gas_pressure.json --only 0,3-5 --events 100
python gem.py scan scans/gas_pressure.json --analyze-only   # refit + retabulate existing results
python gem.py scan configs/baseline.json --param gas.penning.r --values 0.1 0.3   # quick 1-axis scan
```

Output in `runs/scans/<scan name>/`: one folder per point (`gain_events.csv`, `gain_summary.json`, `gain_fit.json/.png`, `config.json`), `scan.csv` (axis values, derived values, mean gain ± error, fitted gain ± error, survival ε, Polya G0/θ, status, field hash), `scan.json` (the resolved plan), `scan.png` (gain vs the first axis, one line per value of the second) and `scan.log`.

How it runs: every point is validated before anything starts (all invalid points are reported at once); the field maps that are missing are built first, each only once (points that differ only in gas, Penning or simulation settings share a map); then the points run one after another, each on `--jobs` workers. Points with results for the same run key are skipped, so an interrupted scan just continues. A failing point is recorded (`status`, `<point>.error.txt`) and the scan goes on (`--stop-on-error` to abort); the exit code is 3 if any point failed.

**Queues.** `scripts/run_scans.sh` runs several scan files in sequence and prints a summary:

```bash
export SIMGEM_PYTHON=~/miniconda3/envs/garfield/bin/python   # a python with gmsh
scripts/run_scans.sh -j 6 scans/voltage.json scans/gas_pressure.json
scripts/run_scans.sh -d -j 6 scans/*.json      # detached (nohup); log in runs/scans/queue-*.log
scripts/run_scans.sh -n scans/*.json           # dry run of the whole queue
```

Example scans: `voltage` (same ΔV on all foils), `voltage_split` (GEM1↔GEM3 at constant total, `zip`), `gas_pressure` (CH₄ fraction × pressure, reduced fields fixed), `geometry` (hole diameter × transfer gap), `penning_r`.

## Gas

The microscopic avalanche needs no gas table: `MediumMagboltz` computes its collision rates from the Magboltz cross sections when tracking starts. A new mixture or pressure is just a config change. `gem.py gastable` still produces a Magboltz table (`gas`, settings in `gas.table`), but only for drift-line methods or cross-checks.

**Penning matching** (`"match": "table"`). Published r values (Sahin et al., Garfield's `"auto"`) are defined through the Townsend coefficient of the table formula alpha (1 + r nu_exc/nu_ion), fitted to measured gain curves. The microscopic avalanche gives 1-6 % less alpha with the same r, which compounds to 10-40 % in a GEM stack. With matching on, `simgem_penning` finds the r_micro that reproduces the table alpha (mean relative mismatch zero over `match_fields_kV_cm`) and the gain runs use it. It runs once per gas, Penning setting and field set, cached in `gas/penning/` (`SIMGEM_PENNING`), and recorded as `penning_r_micro` in `gain_summary.json`. The calibrations committed in the package's own `gas/penning/` are also found when running from another directory (e.g. in the container).

## Gain corrections

Garfield++ underestimates measured GEM gains by about x2 per foil, even with Penning matching. This is an open issue (RD51 WG4; CERN GDD) and SimGEM reproduces it against current-mode single-GEM data ([validation notes](docs/validation.md)). The optional `corrections` block multiplies the reported gains by `(charging_up_per_foil x gain_scale_per_foil) ^ n_foils`. The raw gains are unchanged, corrections never re-run a simulation, and scans add `gain_correction`, `gain_corrected` and `gain_corrected_err` columns. Use `scripts/compare_gorg.py --corrected` to compare with GORG.

`"corrections": "literature"` (needs `"match": "table"`) uses:
- `charging_up_per_foil` = 1.05 (range 1.02-1.10): measured surface charging of pre-polarised double-mask foils. Kapton polarisation is excluded, since the field map uses the static permittivity.
- `gain_scale_per_foil` = 1.85 x/÷ 1.3: current-mode single-GEM data in Ar:CO2 70:30 (GDD 2.09; Bachmann Fig. 24 ~2.2, Fig. 13 1.37) vs SimGEM with matching.

Both are stopgaps until an in-house calibration replaces them: single-GEM gain vs time over hours, in the mixtures, pressures and voltages of interest.

## Misalignment

Fixed offsets go in the foils (`offset_um`). For the statistical study:

```jsonc
"misalignment_mc": {
  "distribution": "gaussian",      // or "uniform" (random relative phase: foils with no alignment at all)
  "sigma_um": 20,
  "foils": null,                   // names to perturb; null = all
  "geometries": 20, "events_per_geometry": 500,
  "reference_events": null,        // nominal geometry, default geometries x events_per_geometry
  "seed": 1
}
```

`gem.py mc` first runs the transport check, then draws the geometries (each a validated, cached field map), runs the gain on each and on the nominal geometry, and reports:

* the mean over geometries vs the nominal geometry (shift ± error),
* the geometry-induced spread: the observed spread of the per-geometry means minus what the avalanche statistics alone produce, with a 95% upper limit, and a χ² test of "no geometry effect".

Interrupted MC runs resume: geometries with results for the same run key are not rerun (use `--force`).

**Does misalignment matter for my detector?** Run `gem.py transport`. Electrons leave a foil at its holes and diffuse by σ_T = D_T·√L across a transfer gap; the hole pattern's imprint on the arrival positions at the next foil falls as exp(−G²σ_T²/2) (G = 4π/(√3 p) for hex). If that "memory" is ≪ 1, arrivals are uniform over the next foil whatever its offset, and alignment cannot change the gain. Electron memory per transfer gap, from `gem.py transport` (140 µm hex foils):

| gas, transfer gap | σ_T | electron memory |
|---|---|---|
| Ar:CH₄ 98:2, 8 bar, 2.1 mm | 156 µm | 10⁻¹⁴ (irrelevant) |
| Ar:CO₂ 70:30, 1 bar, 1.0 mm | 80 µm | 2×10⁻⁴ |
| Ar:CO₂ 70:30, 1 bar, 0.5 mm | 56 µm | 0.01 |
| Ar:CO₂ 50:50, 1 bar, 0.5 mm | 42 µm | 0.10 |
| CO₂, 1 bar, 0.5 mm | 19 µm | 0.63 |

For ions (thermal diffusion) the memory is always large, so alignment matters for ion back-flow even where it cannot change the gain.

## Trajectories (`vis`)

`vis` runs the same avalanche single-threaded and writes where every electron went:

```
runs/<cfg>/vis_tracks.csv    event,track,point,x_cm,y_cm,z_cm,t_ns,e_eV
runs/<cfg>/vis_summary.csv   one row per electron: endpoints, status, weight
runs/<cfg>/vis_meta.json     cell, start/readout planes, electrode heights, hole positions
```

Join on `(event, track)`. A `track` in `vis_tracks.csv` with no summary row is a fragment of the preceding track; `track = -1` in the summary means that electron was drawn as part of the previous one (see `src/vis.cpp`). Points are thinned to `--min-step` (default 5 µm), only the first `--max-tracks` electrons per event get full paths, and `--size-limit` caps the avalanche.

## Speed

* Gain runs on all available cores by default (on an HTCondor grid slot, the CPUs requested for the job); limit with `--jobs N` or `export SIMGEM_JOBS=N` (each worker loads the field map).
* Build Garfield itself in release mode: the avalanche runs inside `libGarfield`.
* Field cost: ~10 s mesh + 20–40 s solve per geometry; field maps are ~40 MB each (`gem.py fields` lists them, `--remove` deletes).
* Microscopic tracking is collision-count bound (∝ pressure × path length). A thinner modeled drift slice and coarser `mesh.coarse_um` cut steps where nothing multiplies.
