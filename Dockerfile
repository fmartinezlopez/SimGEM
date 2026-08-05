# syntax=docker/dockerfile:1
# ---------------------------------------------------------------------------
# GEM simulation stack: Gmsh + Elmer + Garfield++ + SimGEM
# Base = Ubuntu 24.04
# Build with:
#   docker build -t simgem:latest --build-arg SIMGEM_REPO=<url> .
# ---------------------------------------------------------------------------
FROM ubuntu:24.04

LABEL maintainer="frmart@iu.edu" \
      description="GEM simulation: Gmsh + Elmer + Garfield++ + SimGEM"

# use bash so we can `source` ROOT/Garfield setup scripts inside RUN steps
SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ARG ROOT_TGZ=root_v6.40.02.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz
ARG GMSH_TGZ=gmsh-4.15.2-Linux64-sdk.tgz
ARG SIMGEM_REPO=https://github.com/fmartinezlopez/SimGEM.git

# ---- system + Garfield/ROOT runtime deps ----
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git wget ca-certificates gnupg curl \
        gfortran python3 python3-pip \
        libgsl-dev libx11-dev libxpm-dev libxft-dev libxext-dev \
        libssl-dev libtbb-dev libomp-dev \
        libvdt-dev libgif-dev libtiff-dev libgl1-mesa-dev \
        libglu1-mesa \
        software-properties-common \
    && rm -rf /var/lib/apt/lists/* \
    && git config --global http.sslverify false

# ---- Elmer (CSC PPA) ----
RUN add-apt-repository -y ppa:elmer-csc-ubuntu/elmer-csc-ppa \
    && apt-get update \
    && apt-get install -y elmerfem-csc \
    && rm -rf /var/lib/apt/lists/*

# ---- ROOT (precompiled binary) ----
WORKDIR /opt
RUN wget -q https://root.cern/download/${ROOT_TGZ} \
    && tar -xzf ${ROOT_TGZ} && rm ${ROOT_TGZ}
ENV ROOTSYS=/opt/root

# ---- Gmsh SDK ----
RUN wget -q https://gmsh.info/bin/Linux/${GMSH_TGZ} \
    && tar -xzf ${GMSH_TGZ} && mv gmsh-*-sdk gmsh && rm ${GMSH_TGZ}
ENV GMSH_ROOT=/opt/gmsh

# ---- Garfield++ (build against ROOT) ----
RUN source ${ROOTSYS}/bin/thisroot.sh \
    && git clone --depth 1 https://gitlab.cern.ch/garfield/garfieldpp.git \
    && cd garfieldpp && sed -i 's/add_subdirectory(Examples)/# add_subdirectory(Examples)/g' CMakeLists.txt \
    && rm -rf build && mkdir build && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/garfieldpp/install -DCMAKE_BUILD_TYPE=Release .. \
    && make -j"$(nproc)" && make install
ENV GARFIELD_INSTALL=/opt/garfieldpp/install
ENV HEED_DATABASE=/opt/garfieldpp/install/share/Heed/database

# ---- SimGEM ----
RUN source ${ROOTSYS}/bin/thisroot.sh \
    && source ${GARFIELD_INSTALL}/share/Garfield/setupGarfield.sh \
    && git clone ${SIMGEM_REPO} /opt/SimGEM \
    && cd /opt/SimGEM && mkdir build && cd build \
    && cmake .. && make -j"$(nproc)"

# ---- runtime environment ----
ENV PATH=${ROOTSYS}/bin:${GMSH_ROOT}/bin:${GARFIELD_INSTALL}/bin:/opt/SimGEM/build:${PATH}
ENV LD_LIBRARY_PATH=${ROOTSYS}/lib:${GMSH_ROOT}/lib:${GARFIELD_INSTALL}/lib:${GARFIELD_INSTALL}/lib64:${LD_LIBRARY_PATH:-/lib:/usr/lib}

WORKDIR /work
CMD ["/bin/bash"]