# SimGEM

Full triple-GEM stack in one geometry and one field solve, with the gain computed
by **parallel microscopic avalanche** across all CPU cores. Everything is driven
by a JSON config, and **each config works in its own run folder** so you can keep
many configurations side by side.

```
configs/baseline.json  ─►  gem.py  ─►  runs/baseline/
                                         ├── stack.msh, stack/ (mesh + field)
                                         ├── dielectrics.dat, stack.sif
                                         ├── gain_events.csv       (merged, all workers)
                                         ├── gas.log field.log gain.log ...  (per stage)
                                         └── gain_worker_w*.log     (per worker)
```

Each `gem.py` stage writes its own log in the run folder (`gas.log`, `field.log`,
`gain.log`, `scan.log`). Gain results are merged into a single `gain_events.csv`
with columns `worker,event,absolute_gain,effective_gain` — where **absolute_gain**
is the full avalanche size from the one primary electron and **effective_gain** is
how many of those reach the readout (the number averaged for "gain"). `scan`
writes one `gain_events_r<r>.csv` per point plus `penning_scan.csv`.

## One CLI: `gem.py`

Each subcommand takes a config; the run folder is `runs/<config stem>/`.

```bash
python gem.py gas    configs/baseline.json                 # Magboltz gas table (once)
python gem.py field  configs/baseline.json                 # geometry + field map
python gem.py gain   configs/baseline.json                 # PARALLEL stack gain
python gem.py all    configs/baseline.json                 # gas + field + gain
python gem.py scan   configs/baseline.json --r 0.2 0.3 0.4 0.5   # Penning sweep
python gem.py vis    configs/baseline.json --events 3       # trajectories -> CSV
```

Gain always runs in parallel (independent microscopic avalanches, one worker per
core, aggregated). Override with `--events N`, `--jobs J`, `--r R`.

### Trajectories for plotting (`vis`)

`vis` runs the same avalanche as `gain` but single-threaded, and writes where
every electron went instead of just how many arrived:

```
runs/<cfg>/vis_tracks.csv    event,track,point,x_cm,y_cm,z_cm,t_ns,e_eV
runs/<cfg>/vis_summary.csv   one row per electron: endpoints, status, weight
runs/<cfg>/vis_meta.json     cell bounding box + start/readout planes
```

Join the two on `(event, track)`. A `track` in `vis_tracks.csv` with no row in
the summary is a fragment of the preceding track (see `src/vis_stack.cpp` for
why); `track = -1` in the summary means that electron was drawn as part of the
previous one.

Trajectories are enormous at 8 bar — one drift line is ~10^6 collisions — so
points are thinned to `--min-step` (default 5 µm) as they are written, and only
the first `--max-tracks` (default 200) electrons of an event get a full path.
Use `--size-limit` to cap the avalanche and `--dt-max` to tune how drift lines
are split apart.

### Typical first run

```bash
# build once
cmake -S . -B build \
  -DCMAKE_CXX_FLAGS="-D_LIBCPP_DISABLE_AVAILABILITY -D_LIBCPP_HAS_NO_VENDOR_AVAILABILITY_ANNOTATIONS"
cmake --build build -j"$(sysctl -n hw.ncpu)"

python gem.py all configs/baseline.json --events 1000
```

### Multiple configurations

Drop another JSON in `configs/` (e.g. `highfield.json` with different `dV_gem_*`).
It builds and runs in its **own** folder `runs/highfield/`, independent of others.
Configs that share the same `gas` block share one gas file (in `gas/`), so you
only generate it once.

```bash
python gem.py all configs/highfield.json
```

## What you edit in a config

- **gas**: composition, `pressure_bar`, `temperature_K`, `penning_r`, `gas_file`.
- **geometry**: `pitch_um`, `kapton_thickness_um`, hole diameters, `kapton_permittivity`.
- **stack**: the three real gaps + `drift_gap_modeled_mm` (the short slice of the
  18 mm drift actually meshed — the true drift field is applied across it).
- **fields**: gap fields (kV/cm) and ΔV across each GEM (V). GEM1 = top, GEM3 =
  bottom. `gem.py` turns these into absolute electrode potentials automatically.
- **mesh**: `lc_fine_um` near holes, `lc_coarse_um` in gaps (raise to go faster).
- **avalanche**: default `num_events`, `avalanche_size_limit`.

## Look before you leap

```bash
python gem.py field configs/baseline.json --mesh-only   # build mesh, stop
python gem.py field configs/baseline.json --gui         # open it in the gmsh GUI
```

Clip through the stack and confirm three foils, biconical holes, and a graded
mesh (fine at foils, coarse in gaps). After the solve, open
`runs/<name>/stack/case.vtu` in ParaView to check the potential steps across each
GEM and the field peaking in the holes.

## Calibrating Penning `r`

The default `r = 0.18` is the 1-bar / 5%-CH₄ literature value; at 8 bar the
effective transfer is higher. Find your value against the Fe-55 gain curve:

```bash
python gem.py scan configs/baseline.json --r 0.2 0.3 0.4 0.5 --events 1000
```

For a real calibration, repeat at a few `dV_gem_*` settings (separate configs) and
find the single `r` that matches the *slope* of measured gain-vs-voltage, not one
point.

## Speed

- **Parallelism is built in** — gain uses all cores. Set `--jobs` lower if memory
  is tight (each worker loads the field map).
- **Build Garfield optimized** — the avalanche runs inside `libGarfield`, so a
  Release build of Garfield itself is the biggest single win. See INSTALL_macos.md
  / INSTALL.md ("Optimized Garfield build").
- **Trim cost that carries no physics**: smaller `drift_gap_modeled_mm`, coarser
  `lc_coarse_um`. Neither approximates the avalanche; they cut steps in uniform-
  field regions where nothing multiplies. (Microscopic tracking at 8 bar is
  step-count-bound: ~8× the collisions per mm of a 1-bar gas.)

## Modeling choices (change if needed)

Holes aligned across the three foils, square cell, mirror-periodic walls; copper
as a zero-thickness electrode on the kapton faces; drift truncated to the modeled
slice; gain reported is **effective** (primaries below the drift plane, electrons
counted at the readout — folds in collection, all three avalanches, both
transfers, and extraction). Misaligned holes would need offset holes + true
periodic BCs, which this scaffold doesn't set up.

## Files

```
configs/*.json     one file per configuration
gem.py             CLI: gas / field / gain / all / scan / vis
stacklib.py        geometry + field-map construction (imported by gem.py)
src/generate_gas.cpp   gas table (args: ar% ch4% pressure_bar T outfile)
src/gain_stack.cpp     one worker: microscopic stack gain (seed + tag args)
src/vis_stack.cpp      same avalanche, but dumps trajectories to CSV
```

## Misalignment: rigid foil translations

Adding a `misalignment` block to a config switches the run from the legacy
aligned **mirror** cell to a translational-**periodic** cell whose three foils
can be shifted laterally. A rigid whole-foil shift preserves periodicity at the
original pitch, so the infinite misaligned detector is still represented
**exactly** on one pitch-sized cell — no super-cell needed. **Without the block,
nothing changes**: same mirror geometry, no periodic BCs, identical results.

Offsets are per-foil `(dx, dy)`; only the *relative* offsets affect the gain.
"Fully misaligned" = a diagonal half-pitch `(0.5, 0.5)` (a hole sitting over the
solid-kapton point equidistant from four reference holes).

**Fixed scenario** — e.g. GEM2 fully misaligned, GEM1 & GEM3 aligned. Runs like
any normal config (`field` → `gain`), one geometry, one field solve:

```json
"misalignment": {
  "mode": "fixed", "unit": "pitch_fraction",
  "offsets": { "gem1": [0,0], "gem2": [0.5,0.5], "gem3": [0,0] }
}
```
```bash
python gem.py all configs/mis_gem2_full.json --events 2000
```

**Random Monte Carlo** — each foil independently drawn from `N(0, σ)` (σ in
`sigma_um` *or* `sigma_pitch_fraction`). `gem.py mc` builds K geometries, throws
M events at each, pools them, runs a matched aligned baseline, and compares:

```json
"misalignment": {
  "mode": "random", "sigma_um": 20.0,
  "n_configs": 30, "events_per_config": 400, "seed": 12345
}
```
```bash
python gem.py gas configs/mis_random20um.json          # once
python gem.py mc  configs/mis_random20um.json          # K geoms + baseline, pooled
```
`mc` writes `runs/<cfg>/mc_manifest.csv` (per-config offsets + mean, for
gain-vs-offset plots), `mc_pooled.csv` (every event tagged by config+offsets),
and `mc_summary.txt` (pooled-misaligned vs aligned: mean shift, RMS change, and
the config-to-config spread the misalignment adds).
