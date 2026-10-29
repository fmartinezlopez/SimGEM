# Installation

SimGEM needs four tools:

- **ROOT**: Garfield's framework dependency (C++17, GSL, gfortran).
- **Garfield++**: the avalanche engine, which bundles Magboltz. Build it from source, in Release mode. SimGEM is validated with commit `58d22389`.
- **gmsh**, as a Python module: meshing.
- **Elmer** (`ElmerGrid`, `ElmerSolver`): the electrostatics solver.

Beyond the standard library, the Python side needs only gmsh. The plots, the GORG scripts and `studies/` also use numpy, pandas, matplotlib and openpyxl.

Platforms:
- **Ubuntu 22.04 / 24.04** (natively or under WSL2): the simplest route.
- **macOS:** see the macOS section.
- **Container:** see the Container section, which builds everything.

## Linux (Ubuntu / WSL2)

### 0. System build tools

```bash
sudo apt update
sudo apt install -y git cmake build-essential gfortran \
                    libgsl-dev libboost-dev python3-pip
```

### 1. ROOT (via conda — least painful, version-agnostic)

Install Miniforge if you don't have conda, then:

```bash
conda create -n garfield -c conda-forge root gmsh
conda activate garfield
```

That single env gives you both **ROOT** and **Gmsh**. Verify:

```bash
root-config --version     # e.g. 6.34/xx
gmsh --version            # e.g. 4.13.x
```

(Alternative to conda: download a prebuilt ROOT tarball for your exact Ubuntu version from root.cern and `source /path/to/root/bin/thisroot.sh`. Gmsh can then come from `pip install gmsh` or `sudo apt install gmsh`.)

### 2. Elmer (Ubuntu PPA)

```bash
sudo add-apt-repository ppa:elmer-csc-ubuntu/elmer-csc-ppa
sudo apt update
sudo apt install -y elmerfem-csc
```

Verify:

```bash
ElmerGrid            # prints usage
ElmerSolver --help   # prints usage
```

If the PPA has no build for your Ubuntu release, use the Docker route instead, or the official Elmer Docker image.

### 3. Garfield++ (from source)

With the conda env **active** (so ROOT is on PATH):

```bash
cd ~
git clone https://gitlab.cern.ch/garfield/garfieldpp.git
cd garfieldpp
git checkout 58d22389        # optional: the version SimGEM is validated with
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/garfieldpp/install ..
make -j"$(nproc)"
make install
```

This compiles Magboltz (Fortran) and the Heed database into the install tree. After it finishes, the setup script lives at `$HOME/garfieldpp/install/share/Garfield/setupGarfield.sh`.

> **`-DCMAKE_BUILD_TYPE=Release` matters for speed.** The avalanche runs inside
> `libGarfield`, so an unoptimized library is slow regardless of how your own code
> is built. Release (`-O3 -DNDEBUG`) is the single biggest speedup for microscopic
> tracking. If you already built without it, reconfigure the build dir with that
> flag (`rm -rf *` first) and `make install` again.

> Common build errors and fixes:
> - `RootNewMacros.cmake not found` → ROOT env not active. `conda activate garfield`
>   (or source `thisroot.sh`) **before** running cmake.
> - `heed++/... file not found` at runtime → you didn't source `setupGarfield.sh`
>   (it sets `HEED_DATABASE`).
> - Fortran/OpenMP errors → `sudo apt install gfortran libgomp1`.

### 4. Point the project at your install

Edit `env.sh` in this project so `GARFIELD_INSTALL` matches, and uncomment `conda activate garfield` (or your `thisroot.sh` line). Then, in every new shell:

```bash
conda activate garfield      # ROOT + Gmsh
source env.sh                # Garfield + sanity check
```

You should see a line for `root-config`, `gmsh`, `ElmerGrid`, `ElmerSolver`, and a Garfield path — none saying `MISSING`.

---

## macOS (Apple Silicon and Intel)

Starting point: you already have a conda env with **ROOT + Gmsh**. This covers the two remaining pieces — **Elmer** and **Garfield++** — plus the macOS-specific gotchas (Fortran, OpenMP, Apple Silicon, dylib paths). Works on both Apple Silicon and Intel. Building and running SimGEM is the same on both platforms (see Build and test SimGEM).

Confirm the env first:

```bash
conda activate garfield          # whatever you named it
root-config --version
gmsh --version
```

Install Xcode command-line tools if you haven't:

```bash
xcode-select --install
```

---

### 1. Elmer

No apt PPA on macOS. Pick one of these — **MacPorts is the least fussy** if you have it; the Homebrew tap compiles from source and can need an Apple-Silicon workaround.

**Option A — MacPorts (simplest):**

```bash
sudo port install elmerfem
```

**Option B — Homebrew tap (compiles from source, takes a while):**

```bash
brew install gcc cmake open-mpi          # gfortran + build deps
brew tap elmercsc/homebrew-elmerfem
brew install elmer --HEAD
```

On Apple Silicon the tap has historically needed a small patch (tracked in elmerfem issue #313). If `brew install elmer --HEAD` fails on an M-series Mac, fall back to Option C or the Docker image.

**Option C — Docker (escape hatch, or if the above fight you):**
use the official Elmer image, or the project `Dockerfile` (Linux container; runs fine on Mac under Docker Desktop). You'd then run the whole pipeline inside the container.

Verify whichever you chose:

```bash
ElmerGrid            # prints usage
ElmerSolver --help
```

---

### 2. Garfield++ (from source)

The macOS wrinkle is compilers: Magboltz is Fortran (Apple clang can't compile it) and Garfield uses OpenMP (Apple clang doesn't ship it). The clean fix is to build with conda-forge's compilers inside your env — they're consistent with the conda ROOT build and bring OpenMP (`llvm-openmp`) along, so you avoid the Apple-clang/libomp dance entirely.

```bash
conda activate garfield
conda install -c conda-forge compilers cmake make gsl   # + GSL (see note below)
```

Activating the env now sets `CC`, `CXX`, and `FC` to the conda toolchain. Build:

```bash
cd ~
git clone https://gitlab.cern.ch/garfield/garfieldpp.git
cd garfieldpp
mkdir build && cd build
cmake -DGARFIELD_WITH_GSL=ON \
      -DCMAKE_CXX_FLAGS="-D_LIBCPP_DISABLE_AVAILABILITY -D_LIBCPP_HAS_NO_VENDOR_AVAILABILITY_ANNOTATIONS" \
      -DCMAKE_INSTALL_PREFIX="$HOME/garfieldpp/install" ..
make -j"$(sysctl -n hw.ncpu)"
make install
```

> ✓ **Confirmed working** on Apple Silicon with the conda-forge LLVM-20 toolchain
> (clang 20 / libc++ 20, gfortran 14, ROOT + Gmsh from conda-forge). The three
> flags above — `GARFIELD_WITH_GSL`, and the two `_LIBCPP_..._AVAILABILITY`
> defines — are all required together on this setup; drop any one and the build
> fails at neBEM or at the first `<charconv>` include.

**Full copy-paste recipe** (from an activated conda env with root + gmsh):

```bash
conda activate garfield
conda install -c conda-forge compilers cmake make gsl
cd ~
git clone https://gitlab.cern.ch/garfield/garfieldpp.git
cd garfieldpp && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DGARFIELD_WITH_GSL=ON \
      -DCMAKE_CXX_FLAGS="-D_LIBCPP_DISABLE_AVAILABILITY -D_LIBCPP_HAS_NO_VENDOR_AVAILABILITY_ANNOTATIONS" \
      -DCMAKE_INSTALL_PREFIX="$HOME/garfieldpp/install" ..
make -j"$(sysctl -n hw.ncpu)" && make install
```

> **`-DCMAKE_BUILD_TYPE=Release` matters for speed.** The avalanche runs *inside*
> `libGarfield`, so an unoptimized Garfield library is slow no matter how your own
> code is compiled. Release gives `-O3 -DNDEBUG` and is the single biggest speedup
> for microscopic tracking — well worth the compile time. If you already built
> Garfield without it, just reconfigure the existing build dir with the line above
> (add `rm -rf *` first) and `make install` again.

> **Why the `_LIBCPP_..._AVAILABILITY` flags?** conda-forge's newest (LLVM-20)
> toolchain ships its own libc++ that *has* floating-point `std::from_chars`, but
> the header wrongly applies Apple **system**-libc++ availability markup and
> expands to a malformed attribute (`__attribute__((macos, strict, introduced =
> 99.0, …))`). Any file that includes `<chrono>`/`<charconv>`
> (e.g. `AvalancheMicroscopic.cc`, `ComponentGrid.cc`) then fails with
> `expected ')'`. These two defines switch the spurious markup off — safe, because
> you link conda's own libc++, which has the symbol. Note that
> `-DCMAKE_OSX_SYSROOT=…` does **not** help: conda's clang wrapper injects its own
> `-isysroot` from the `CONDA_BUILD_SYSROOT` env var, overriding the CMake flag. If
> you'd rather fix the SDK than disable availability, set that env var instead:
> `export CONDA_BUILD_SYSROOT="$(xcrun --show-sdk-path)"` (aim for SDK 15.x). As a
> heavier fallback, drop the toolchain below LLVM 20:
> `conda install -c conda-forge "clang=19" "clangxx=19" "libcxx=19"`.

> **`-DGARFIELD_WITH_GSL=ON` is effectively mandatory on macOS.** Garfield's neBEM
> component (`Isles.c`) calls the C++17 special-math function `std::comp_ellint_1`,
> which macOS/conda **clang's libc++ does not implement** (only GNU libstdc++ does).
> Without the flag the build dies at `Isles.c` with
> `no member named 'comp_ellint_1' in namespace 'std'`. The flag routes those
> elliptic integrals through GSL instead, so it compiles. (Fix confirmed by the
> Garfield++ maintainer on the ROOT forum.) If you re-run cmake after a failed
> attempt, `rm -rf *` in the build dir first — the cache remembers the old config.

The setup script ends up at `$HOME/garfieldpp/install/share/Garfield/setupGarfield.sh`.

> **Alternative (system clang + Homebrew), if you'd rather not add conda compilers:**
> ```bash
> brew install gcc libomp cmake
> cmake -DCMAKE_INSTALL_PREFIX=$HOME/garfieldpp/install \
>       -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
>       -DOpenMP_CXX_LIB_NAMES="omp" \
>       -DOpenMP_omp_LIBRARY="$(brew --prefix libomp)/lib/libomp.dylib" ..
> ```
> This works but is the thing that most often needs fiddling; the conda-compilers
> route is more reliable, especially on Apple Silicon.

---

### 3. Point the project at your install

Edit `env.sh`:
- uncomment `conda activate garfield` (on macOS with conda ROOT you do **not** need `thisroot.sh` — conda activation sets ROOT up);
- set `GARFIELD_INSTALL="$HOME/garfieldpp/install"`.

Then, in every new shell:

```bash
source env.sh
```

You should see non-MISSING lines for `root-config`, `gmsh`, `ElmerGrid`, `ElmerSolver`, and the Garfield path. On macOS `setupGarfield.sh` sets `DYLD_LIBRARY_PATH` so the Garfield dylib is found at runtime — that's why sourcing it every session matters.

---

## Build and test SimGEM

From the project root, with the environment set up (`source env.sh`):

```bash
cmake -S . -B build && cmake --build build -j
# builds simgem_gain, simgem_vis, simgem_transport, simgem_gastable, simgem_penning, simgem_counter
# macOS + conda LLVM 20: add
#   -DCMAKE_CXX_FLAGS="-D_LIBCPP_DISABLE_AVAILABILITY -D_LIBCPP_HAS_NO_VENDOR_AVAILABILITY_ANNOTATIONS"
python -m unittest discover tests            # fast checks of the Python side
python gem.py check configs/baseline.json    # validate a config, print the stack
python gem.py run   configs/baseline.json --events 200   # field map + gain
```

Use the Python that has `gmsh` (the conda env), not the system one. Outputs go to `runs/` and `fields/` under the working directory (`SIMGEM_RUNS`, `SIMGEM_FIELDS`). The programs are found in `build/`, or in `SIMGEM_BIN`, or on `PATH`.

---

## First-run smoke test (before your real geometry)

To confirm the whole chain is wired up, run the **official** GEM example that ships with Garfield (its field map is included, so no Gmsh/Elmer needed):

```bash
cp -r $GARFIELD_INSTALL/share/Garfield/Examples/Gem ~/gem_test
cd ~/gem_test && mkdir build && cd build
cmake .. && make -j && ./gem
```

If that runs and reports an avalanche, ROOT + Garfield are healthy and any remaining problems are in the Gmsh/Elmer geometry, not the toolchain.

---

## Container

`Dockerfile` builds a runtime image on AlmaLinux 9 (x86_64). Everything is compiled in a first stage; the final image has only what running SimGEM needs, with no compilers, headers, sources or docs:
- ROOT from EPEL, only the libraries Garfield++ links;
- Garfield++ at the validated commit 58d22389;
- Elmer (`ElmerGrid`, `ElmerSolver`), release 26.2.1, serial, without GUI or MPI;
- the headless gmsh shared library and Python module (no GUI, so `gem.py field --gui` does not work in the container);
- Python 3.12 with the standard library only. `gem.py` needs nothing else; plots are skipped without matplotlib, so make them (and run `scripts/` and `studies/`) outside the container;
- SimGEM, cloned from GitHub at `SIMGEM_REF` (a release tag; default `main`): the package, configs, scans, docs and tests, the programs in `/opt/SimGEM/build`, and a `simgem` command that runs `gem.py`. `/opt/SimGEM/VERSION` records the SimGEM tag and commit and the Garfield++, Elmer, gmsh and ROOT versions.

```bash
docker build -t simgem:v1.0.0 --build-arg SIMGEM_REF=v1.0.0 .   # Apple Silicon: add --platform linux/amd64 (slow)
docker run --rm -it -v "$PWD":/work -w /work simgem:v1.0.0 bash
simgem run /opt/SimGEM/configs/baseline.json                      # inside the container
```

Only pushed tags can be built: commit, tag and push first. The image is read only, so run from a writable directory; `runs/`, `fields/` and `gas/` are created there. The Penning calibrations shipped in `/opt/SimGEM/gas/penning/` are found from any directory.

`.github/workflows/build-sif.yml` does this on GitHub. Pushing a tag `v*` (or running the workflow by hand from the Actions tab, for any tag or branch):
1. builds the image for that tag;
2. smoke-tests it: `simgem check`, then a short `simgem run` of `configs/validation/garfield_gem_example.json` (gmsh mesh, Elmer solve, field-map check, 4 avalanches);
3. pushes it to `ghcr.io/fmartinezlopez/simgem:<tag>`, plus `latest` for tags;
4. converts it to `simgem-<tag>.sif`, kept as a workflow artifact for 30 days.

---

## Troubleshooting

### Linux

| Symptom | Cause / fix |
|---|---|
| `cmake` can't find Garfield | `source env.sh`; check `GARFIELD_INSTALL`. |
| `RootNewMacros.cmake` missing | ROOT not active when building Garfield. Activate conda / source `thisroot.sh`, rebuild. |
| `HeedParticle.h` / heed db errors | `source setupGarfield.sh` (env.sh does this). |
| `ElmerGrid: command not found` | PPA not installed, or new shell without PATH. Re-run apt install. |
| `could not load the field map` | The field build failed or was removed; rerun `gem.py field <config> --force` and read `fields/.failed-<hash>/elmersolver.log`. |
| `field map FAILED validation` | Periodic walls or potentials wrong; the failing map is kept in `fields/.failed-<hash>/` for inspection. |
| `Cannot find 'simgem_gain'` | Build with cmake, or set `SIMGEM_BIN` to the folder holding the programs. |
| Avalanche never finishes | Gain too high / runaway. Lower `simulation.avalanche_size_limit` (truncated events are flagged). |

### macOS

| Symptom | Fix |
|---|---|
| `expected ')'` in `from_chars_floating_point.h` / `introduced = 99.0` | conda LLVM-20 libc++ emits bad availability markup. Add `-DCMAKE_CXX_FLAGS="-D_LIBCPP_DISABLE_AVAILABILITY -D_LIBCPP_HAS_NO_VENDOR_AVAILABILITY_ANNOTATIONS"`. (`-DCMAKE_OSX_SYSROOT` won't help — set `CONDA_BUILD_SYSROOT` env var instead, or downgrade to `clang/clangxx/libcxx=19`.) |
| `no member named 'comp_ellint_1' in namespace 'std'` (neBEM/Isles.c) | libc++ lacks C++17 special-math functions. Reconfigure with `-DGARFIELD_WITH_GSL=ON` after `conda install -c conda-forge gsl` and `rm -rf *` in the build dir. |
| `No CMAKE_Fortran_COMPILER` / Magboltz won't compile | No Fortran compiler. `conda install -c conda-forge compilers` (or `brew install gcc`), reconfigure. |
| `Could NOT find OpenMP` | Use the conda-compilers route, or pass the explicit `libomp` flags shown above. |
| `ModuleNotFoundError: No module named 'gmsh'` | You are on the system Python; `conda activate garfield` (or call its `python`). |
| `library not loaded: libGarfield...dylib` at runtime | `source env.sh` (sets `DYLD_LIBRARY_PATH` via setupGarfield.sh) in the shell you run from. |
| `RootNewMacros.cmake` missing when building Garfield | conda env not active during cmake. `conda activate garfield`, rebuild. |
| Elmer `brew install` fails on M1/M2/M3 | Use MacPorts (`sudo port install elmerfem`) or the Docker image. |
| Downloaded Elmer binary "cannot be opened, unidentified developer" | `xattr -dr com.apple.quarantine <path>`, or right-click → Open once. |
| `ranlib: file ... has no symbols` during Elmer build | Harmless warning; ignore. |
| ABI/link errors mixing brew-gcc Garfield with conda ROOT | Don't mix. Build Garfield's C++ with clang (conda or Apple), use gfortran only for the Fortran; the conda-compilers route keeps everything consistent. |
