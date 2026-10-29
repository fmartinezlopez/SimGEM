#!/usr/bin/env python3
"""
SimGEM -- microscopic simulation of a full GEM stack (gmsh + Elmer + Garfield++).

  python gem.py check     configs/baseline.json      # validate, show the stack
  python gem.py field     configs/baseline.json      # mesh + solve + validate (cached)
  python gem.py gain      configs/baseline.json      # parallel stack gain
  python gem.py run       configs/baseline.json      # field (if needed) + gain
  python gem.py scan      scans/voltage.json         # parameter scan from a scan file
  python gem.py scan      configs/baseline.json --param gas.penning.r --values 0.1 0.2 0.3
  python gem.py analyze   runs/baseline              # Polya fit of the gain distribution
  python gem.py mc        configs/misalign_mc.json   # foil-misalignment Monte Carlo
  python gem.py transport configs/baseline.json      # diffusion + misalignment sensitivity
  python gem.py vis       configs/baseline.json      # trajectories -> CSV
  python gem.py gastable  configs/baseline.json      # optional Magboltz table
  python gem.py fields                               # list cached field maps

Run `python gem.py <command> -h` for options. See README.md.
"""
from simgem.cli import main

if __name__ == "__main__":
    main()
