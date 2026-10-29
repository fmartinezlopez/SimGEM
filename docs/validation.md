# Validation

How far SimGEM's gains can be trusted, checked against an independent field map, gas data and measured GEM gains (October 2026, Garfield++ 58d22389).

## Summary

- **Field and avalanche chain.** gmsh + Elmer reproduce the ANSYS field map of Garfield's `Examples/Gem` to a few percent, and identical avalanches on the two maps give the same gain. The mesh is converged.
- **Gas model.** Magboltz plus Garfield's Penning parameterisation reproduces measured proportional-counter gains (Ar:CO2) to a few percent.
- **Penning definition.** Microscopic tracking with a published r gives 1-6 % less Townsend coefficient than the table formula that r was fitted with. `"match": "table"` corrects this (`simgem_penning`).
- **Single-GEM deficit.** Against current-mode single-GEM data, SimGEM is low by about x2.2 per foil even with matching, as are Garfield++ simulations elsewhere: an open issue in the field (RD51 WG4, CERN GDD).
- **Corrections.** `"corrections": "literature"` (charging-up 1.05 and scale 1.85 per foil) is derived from those data alone.
- **GORG.** The triple-GEM shortfall (x7.5-9) is that per-foil deficit, compounded. With matching and the literature corrections, SimGEM is 1.3-1.7 above the GORG data, inside the corrections' uncertainty, with matching slopes.

## 1. Field and avalanche chain vs Garfield's `Examples/Gem`

`configs/validation/garfield_gem_example.json` reproduces the example:
- one standard GEM: 140 um hex pitch, 50 um kapton with eps = 4.0, 5 um Cu, 70/50 um holes, 80 um rim;
- drift and induction planes 1 mm from the kapton centre;
- potentials -200 / -150 / +150 / +350 V (dV = 300 V);
- Ar:CO2 80:20 at 760 Torr, Penning r = 0.51 on Ar.

Results:
- **Field.** Along the hole axis and 15/20 um off axis, SimGEM agrees with ANSYS point by point to a few percent (peak |E| on axis 37.1 vs 37.0 kV/cm). Halving the mesh size changes it by < 1 %.
- **Avalanches,** 300 events on each map, same gas and start positions:

| field map | absolute gain | effective gain | zero fraction |
|---|---|---|---|
| ANSYS (Garfield example) | 8.15 +- 0.41 | 3.42 +- 0.21 | 27 % |
| SimGEM (gmsh + Elmer) | 8.13 +- 0.39 | 3.26 +- 0.20 | 25 % |

Rerun this after any change to the geometry, mesh or field code. Run the example's `gem.C` avalanche on the ANSYS files of `Examples/Gem` and on a SimGEM map (`ComponentElmer`, as in `src/common.hh`) with the same seed.

## 2. Gas model

**Proportional counters (Sahin et al., NIM A 768 (2014) 104).** The method is the paper's: Magboltz table + Garfield's Penning parameterisation, integrated along the field of the 24 um wire (r_c = 1.25 cm), with `simgem_counter`. It reproduces the measured current-based gains of Ar:CO2 50:50 to a few percent:

| p | V | SimGEM | measured |
|---|---|---|---|
| 0.79 atm | 1.2-2.4 kV | 14.8 / 74 / 467 / 3545 / 31300 | ~14 / 72 / 460 / 3700 / 33000 |
| 1.19 atm | 1.6-2.8 kV | 22 / 160 / 1560 / 19700 | ~22 / 165 / 1800 / 23000 |

At 1.19 atm the last two measured points are where Sahin reports photon feedback.

**Ar:CH4 at high pressure (Z. Ye et al., as fitted in Sahin et al. 2010, Fig. 11).** Ar:CH4 90:10, r_tube = 1.5 cm, r_anode = 12.5 um, Garfield's r (0.248 at 3.5 atm, 0.275 at 7 atm); data read off the figure.

| p | V | G, no Penning | G, Penning | data | sim/data |
|---|---|---|---|---|---|
| 3.5 atm | 1690 / 1883 / 2059 / 2340 | 13 / 27 / 55 / 196 | 44 / 117 / 306 / 1612 | 82 / 261 / 628 / 3880 | 0.54 / 0.45 / 0.49 / 0.42 |
| 7 atm | 2493 / 2790 / 3239 | 10.8 / 24 / 101 | 48 / 159 / 1221 | 70 / 386 / 3760 | 0.69 / 0.41 / 0.32 |

Sahin's own no-Penning curves are x1.33-1.6 above ours. That is about the normalisation g = 1.25 he applied to this data set, plus a few percent lower alpha in today's Magboltz at the highest fields. His per-pressure fit also needs more Penning transfer than the global r(p, c) gives: r ~ 0.30-0.33 instead of 0.248 at 3.5 atm. Near 1 atm the Ar:CH4 model is good to that level; at 7-10 bar it is the largest gas-model uncertainty.

**Microscopic vs table Townsend coefficient (uniform field).**

| gas, field | r | table alpha (/cm) | microscopic alpha (/cm) |
|---|---|---|---|
| Ar:CO2 70:30, 30 / 45 / 60 / 80 kV/cm | 0 | 132 / 382 / 690 / 1146 | 128 / 374 / 696 / 1117 (+-2 %) |
| same | 0.57 | 210 / 550 / 931 / 1470 | 199 / 511 / 872 / 1382 |
| Ar:CO2 90:10, 30 / 45 / 60 kV/cm | 0.427 | 353 / 752 / 1170 | 336 / 705 / 1120 |
| Ar:CH4 92:8, 40 / 60 kV/cm | 0 | 463 / 948 | 468 / 942 |
| same | 0.205 | 609 / 1160 | 578 / 1143 |
| same | 0.4 | 749 / 1361 | 740 / 1307 |

Without Penning transfer the two agree. With it, microscopic tracking gives less extra ionisation than the table formula alpha (1 + r nu_exc/nu_ion), which is how published r values are defined (Sahin's eq. 3.1). The deficit is about 6 % in Ar:CO2 70:30 and 1-5 % in Ar:CH4. This is a difference of definition, not a wrong r. `simgem_penning` finds the microscopic r that removes it:
- Ar:CH4 92:8, 1 atm: 0.205 -> 0.235-0.238;
- Ar:CH4 92:8, 22 psi: 0.210 -> 0.245;
- Ar:CO2 70:30, 1 atm: 0.574 -> 0.660.

Its ionisation-slope estimate of alpha agrees with the table to 0.3-1.4 % without Penning. The gain-ratio estimate in the table above, alpha = ln(G_2d / G_d) / d with 400 avalanches per distance, was used before `simgem_penning` existed.

## 3. Single GEM vs current-mode data

These measurements take the gain from currents (G = I / (e n R)), so no pulser calibration is involved. All are standard CERN GEMs in Ar:CO2 70:30 at 1 atm.

**Bachmann et al., NIM A 438 (1999) 376** (`scans/validation/bachmann1999.json`, 70/55 um holes; values read off the figures):

| point | SimGEM effective | measured | sim/data | SimGEM absolute | measured real |
|---|---|---|---|---|---|
| 500 V, E_D 1, E_I 3 kV/cm (Fig. 13) | 220 +- 12 | ~350 | ~0.63 | 659 +- 35 | ~900 |
| 400 V, E_D 3, E_I 4.4 kV/cm (Fig. 24) | 26.8 +- 1.9 | ~75 | ~0.36 | 66 +- 5 | - |
| 450 V, E_D 3, E_I 4.4 kV/cm (Fig. 24) | 88 +- 5 | ~210 | ~0.42 | 228 +- 12 | - |

- **Electron losses are modelled correctly.** The simulated effective/absolute ratio (0.33 at 500 V) matches Bachmann's "about 1/3" and his electrode currents (I_S/I_TOT ~ 0.31).
- **Slopes agree:** 23.8 +- 1.8 vs 22.4-23.7 per kV.
- **The two data sets disagree.** Fig. 24 comes from another set-up and sits 1.9-2.4x above Fig. 13, more than the field change explains (~1.5x in his currents, 1.3x in SimGEM). The GDD data side with Fig. 24.

**CERN GDD (D. Janssens, CERN summer student report 2019, repository.cern/records/a5zd2-n0x93)** (`scans/validation/gdd_janssens2019*.json`):
- double-mask foil, 70/50 um holes, E_D = 2, E_I = 3 kV/cm;
- Fe-55, Keithley anode current;
- fit G_eff = exp(0.02162 V - 4.503), measured after charging-up.

| dV | SimGEM | SimGEM, Penning matched | GDD | sim/data (matched) |
|---|---|---|---|---|
| 350 V | 9.6 +- 0.6 | 11.5 +- 0.9 | 21.4 | 0.54 |
| 400 V | 28.8 +- 1.8 | 27.5 +- 1.7 | 63.1 | 0.44 |
| 450 V | 70 +- 4 | 85 +- 5 | 186 | 0.46 |
| 500 V | 212 +- 11 | 233 +- 14 | 548 | 0.43 |
| mean | 0.41 (x2.44) | | | 0.456 +- 0.015 (x2.19) |

Slopes: 20.4 +- 0.5 vs 21.6 +- 0.3 /kV. Matching raises the single-GEM gain by x1.11.

The deficit is not specific to SimGEM:
- the GDD's own Garfield++ (ANSYS map, r = 0.56) is low by x2.25 +- 0.26;
- Croci's thesis data show the same;
- RD51 WG4 lists "simulation underestimates the gain in a single GEM foil by a factor ~2; this propagates also to multi-GEM structures" as open;
- BESIII scales Garfield++ triple-GEM gains by x4.5-6.8 (arXiv:1904.06142).

The GDD study ruled out surface currents, secondary emission and the straight-line stepping (3%), and found no single cause.

## 4. Gain corrections (`"corrections": "literature"`)

**Charging-up, per foil:** gain after saturation over gain at switch-on.
- Hauer et al., NIM A 976 (2020) 164205: double-mask, Ar:CO2 90:10, pre-polarised, 1.04-1.09 (Fe-55: 1.01-1.03).
- Correia et al., JINST 9 (2014) P07025: data 1.07; their simulation 1.10-1.15.
- Chatterjee et al. (arXiv:2007.11444, 2107.00890, 2206.10876): pre-polarised double-mask triple 1.12 in total (~1.04 per foil); single-mask triples 1.03-1.16.
- Kapton polarisation after switch-on is larger and inconsistent between groups: Kumar et al. see x1.6-1.85 per foil, Chatterjee et al. a dip of -5 to -10 %, and Janssens' x1.43 includes it. The field map uses the static permittivity, i.e. the steady polarised state, so polarisation is not corrected for.

Choice: **1.05 (1.02-1.10)**.

**Per-foil scale.** Measured / simulated effective gain per foil, after matching and charging-up:
- GDD: 2.19 / 1.05 = 2.09;
- Bachmann Fig. 24: ~2.2;
- Bachmann Fig. 13: 1.37.

Choice: **log-mean 1.85, scatter x/÷1.3**. The BESIII tuning (1.65-1.9 per foil) agrees. Together the two factors give 1.94 per foil and x7.3 on a triple GEM, uncertain to x/÷ ~2.2. GORG was deliberately not used.

Both are placeholders for an in-house calibration: single-GEM gain against time over hours, in the mixtures, pressures and voltages of interest.

## 5. GORG (TECHTRA triple GEM)

**Uncorrected:** SimGEM is low by a voltage-independent factor (data ADC-offset corrected):

| gas | sim/data | slopes, SimGEM vs GORG (/kV) |
|---|---|---|
| Ar:CH4 92:8, 14.5 psi | 0.134 +- 0.006 | 8.32 +- 0.52 vs 8.13 +- 0.15 |
| Ar:CH4 92:8, 22 psi | 0.157 +- 0.006 | 7.14 +- 0.16 vs 7.33 +- 0.13 |
| Ar:CO2 90:10, 15 psi | 0.111 +- 0.006 | 7.8 +- 0.7 vs 6.3 +- 0.9 (3 points) |

The gas model is not the cause:
- Scaling the Townsend coefficient by k multiplies ln G (and the slope) by k. Closing x7.5 needs k = 1.2-1.35, which would raise the 22 psi slope by 20-35 %, against a measured sim/data of 0.97 +- 0.03.
- The same factor appears in two gases.

Per foil the factor is 0.48-0.51, the single-GEM deficit of section 3.

**With matching and literature corrections** (x7.33; data ADC-offset corrected, one run per voltage, `scripts/compare_gorg.py --corrected --data-at-sim`):

| V_GEM (14.5 psi) | 2400 | 2450 | 2500 | 2550 | 2600 | 2650 |
|---|---|---|---|---|---|---|
| sim/data | 1.10 | 1.26 | 1.44 | 1.51 | 1.43 | 1.39 |

- 14.5 psi: mean 1.34 +- 0.06; slopes 8.88 +- 0.50 vs 8.13 +- 0.15 /kV.
- 22 psi (3000, 3100 V): 1.75, 1.69.
- Matching alone raises the raw triple-GEM gain by x1.36 (14.5 psi) to x1.5 (22 psi).

As a calibration, GORG would want 1.62-1.77 per foil instead of 1.94. That is within the uncertainty. Possible reasons:
- gas and field dependence of the per-foil scale;
- the charging-up state;
- GORG's pulser charge.

## Reproducing

```bash
python gem.py scan scans/validation/gdd_janssens2019_match.json       # section 3
python gem.py scan scans/validation/bachmann1999.json
build/simgem_counter --gas ar:90,ch4:10 --pressure-torr 2660 --temperature-k 293.15 \
    --penning-r auto --anode-radius-cm 0.00125 --cathode-radius-cm 1.5 \
    --voltages 1690,1883,2059,2340                                  # section 2 (Ye, 3.5 atm)
python gem.py scan scans/gorg/gain_curve_14p5psi.json --set gas.penning.match=table \
    --set corrections=literature --out runs/scans/gorg/matched_14p5psi
python scripts/compare_gorg.py runs/scans/gorg/matched_14p5psi --corrected --data-at-sim --adc-offset 18.5
```
