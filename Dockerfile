# syntax=docker/dockerfile:1
# ---------------------------------------------------------------------------
# SimGEM runtime image: AlmaLinux 9 (x86_64), gmsh + Elmer + Garfield++ + SimGEM.
#
# Two stages: everything is compiled in `build`; the final image gets only what
# running SimGEM needs (no compilers, headers, sources or docs):
#   ROOT     EPEL RPMs, only the libraries Garfield++ links (+ their deps)
#   Garfield GARFIELD_REF, the commit SimGEM is validated with (docs/validation.md)
#   Elmer    ELMER_REF, serial, without GUI or MPI
#   gmsh     the headless (nox) SDK: shared library and Python module only
#   SimGEM   SIMGEM_REF (a release tag) from GitHub: package, configs, scans,
#            docs, tests and the compiled programs
# Python is the system python3.12 with the standard library only: gem.py
# needs nothing else, and plots are skipped when matplotlib is missing.
#
# Build with:
#   docker build -t simgem:v1.0.0 --build-arg SIMGEM_REF=v1.0.0 .
# On Apple Silicon add --platform linux/amd64 (emulated, slow).
# ---------------------------------------------------------------------------
ARG ALMA=9

# ===========================================================================
FROM almalinux:${ALMA} AS build
SHELL ["/bin/bash", "-c"]

ARG GARFIELD_REF=58d223891c06337b5014970ccc2595d7f19d5617
ARG ELMER_REF=release-26.2.1
ARG GMSH_TGZ=gmsh-4.15.2-nox-Linux64-sdk.tgz
ARG SIMGEM_REPO=https://github.com/fmartinezlopez/SimGEM.git
ARG SIMGEM_REF=main

RUN dnf install -y dnf-plugins-core epel-release \
    && dnf config-manager --set-enabled crb \
    && dnf install -y --setopt=install_weak_deps=False \
        gcc-c++ gcc-gfortran make cmake git tar gzip \
        openblas-devel eigen3-devel \
        root-core root-io root-tree root-mathcore root-matrix root-hist \
        root-graf root-graf-gpad root-graf3d root-geom root-gdml \
    && dnf clean all

# ---- Garfield++ ----
RUN git clone --filter=blob:none https://gitlab.cern.ch/garfield/garfieldpp.git /src/garfield \
    && cd /src/garfield && git checkout ${GARFIELD_REF} \
    && cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/opt/garfield -DCMAKE_BUILD_TYPE=Release \
             -DGARFIELD_WITH_EXAMPLES=OFF -DGARFIELD_WITH_TESTS=OFF -DGARFIELD_WITH_CUDA=OFF \
    && cmake --build build -j"$(nproc)" && cmake --install build

# ---- Elmer (ElmerGrid + ElmerSolver) ----
RUN git clone --depth 1 --branch ${ELMER_REF} https://github.com/ElmerCSC/elmerfem.git /src/elmer \
    && cmake -S /src/elmer -B /src/elmer/build -DCMAKE_INSTALL_PREFIX=/opt/elmer \
             -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=OpenBLAS \
             -DWITH_MPI=OFF -DWITH_OpenMP=OFF -DWITH_ELMERGUI=OFF -DWITH_ElmerIce=OFF \
             -DWITH_LUA=OFF -DWITH_Mumps=OFF -DWITH_Hypre=OFF -DWITH_Trilinos=OFF \
             -DWITH_Zoltan=OFF -DWITH_MMG=OFF -DWITH_NETCDF=OFF \
    && cmake --build /src/elmer/build -j"$(nproc)" && cmake --install /src/elmer/build \
    && rm -rf /opt/elmer/share/doc /opt/elmer/include

# ---- gmsh (headless): shared library + Python module ----
RUN curl -fsSL https://gmsh.info/bin/Linux/${GMSH_TGZ} | tar -xz -C /src \
    && mkdir -p /opt/gmsh/lib && cp -a /src/gmsh-*-sdk/lib/libgmsh.so* /src/gmsh-*-sdk/lib/gmsh.py /opt/gmsh/lib/

# ---- SimGEM: tracked files at SIMGEM_REF + the compiled programs ----
RUN source /opt/garfield/share/Garfield/setupGarfield.sh \
    && git clone --depth 1 --branch ${SIMGEM_REF} ${SIMGEM_REPO} /src/simgem \
    && cmake -S /src/simgem -B /src/simgem/build && cmake --build /src/simgem/build -j"$(nproc)" \
    && mkdir -p /opt/SimGEM/build \
    && git -C /src/simgem archive HEAD | tar -x -C /opt/SimGEM \
    && cp /src/simgem/build/simgem_* /opt/SimGEM/build/ \
    && { echo "SimGEM ${SIMGEM_REF} $(git -C /src/simgem rev-parse HEAD)"; \
         echo "Garfield++ ${GARFIELD_REF}"; echo "Elmer ${ELMER_REF}"; echo "gmsh ${GMSH_TGZ}"; \
         rpm -q root-core; } > /opt/SimGEM/VERSION \
    && rm -rf /opt/garfield/include /opt/garfield/lib*/cmake /opt/garfield/share/Garfield/Examples

# ===========================================================================
FROM almalinux:${ALMA}-minimal
LABEL maintainer="frmart@iu.edu" \
      description="SimGEM runtime: gmsh + Elmer + Garfield++ + SimGEM"

# runtime libraries only: ROOT (pulls its own deps), Fortran/OpenMP/BLAS for
# Garfield and Elmer; the headless libgmsh needs only zlib and libgomp
RUN microdnf install -y epel-release \
    && microdnf install -y --enablerepo=crb --nodocs --setopt=install_weak_deps=0 \
        python3.12 libgfortran libgomp openblas-serial \
        root-core root-io root-tree root-mathcore root-matrix root-hist \
        root-graf root-graf-gpad root-graf3d root-geom root-gdml zlib \
    && microdnf clean all \
    && ln -s /usr/bin/python3.12 /usr/local/bin/python3

COPY --from=build /opt /opt

ENV GARFIELD_INSTALL=/opt/garfield \
    HEED_DATABASE=/opt/garfield/share/Heed/database \
    ELMER_HOME=/opt/elmer \
    SIMGEM_BIN=/opt/SimGEM/build \
    PYTHONPATH=/opt/gmsh/lib \
    PATH=/opt/SimGEM/build:/opt/elmer/bin:/opt/garfield/bin:/usr/local/bin:/usr/bin:/bin \
    LD_LIBRARY_PATH=/opt/garfield/lib64:/opt/garfield/lib:/opt/elmer/lib/elmersolver:/opt/gmsh/lib:/usr/lib64/root

RUN printf '#!/bin/sh\nexec python3 /opt/SimGEM/gem.py "$@"\n' > /usr/local/bin/simgem \
    && chmod +x /usr/local/bin/simgem \
    && cd /opt/SimGEM && python3 -m unittest discover -s tests

# Outputs (runs/, fields/, gas/) go to the working directory: run from a
# writable directory, e.g. the grid job's scratch area.
WORKDIR /work
CMD ["/bin/bash"]
