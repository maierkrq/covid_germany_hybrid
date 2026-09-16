# syntax=docker/dockerfile:1
#
# Two-stage build:
#   build   - full compiler toolchain + trimmed Kaskade dependencies, compiles
#             the covid_germany_abm_pde_ode tutorial binary.
#   runtime - slim image containing only the compiled binary and the shared
#             libraries it needs at runtime.
#
# The Kaskade dependency tree (./deps/) must be downloaded manually beforehand
# from lakeFS onto this machine - see README/plan for the lakectl commands.
# No lakeFS credentials or network access are required to build or run either
# stage. Confidential model input data (work/input_data) and run output
# (work/output) are never part of the build context; they are supplied to the
# running container via bind mounts.

ARG UBUNTU_VERSION=20.04

# Alternative to PROGRAM=covid which is the default model binary built from the tutorial source code. Other options are: no_zero_covid, create_trajectories.
ARG PROGRAM=covid

########################################################################
# build stage
########################################################################
FROM ubuntu:${UBUNTU_VERSION} AS build
ARG PROGRAM

ENV DEBIAN_FRONTEND=noninteractive

ENV PROJECT_ROOT=/root/covid_germany_hybrid
ENV KASKADE_ROOT=/root/covid_germany_hybrid/kaskade7_test

RUN apt-get update && apt-get install -y \
    libyaml-cpp-dev \
    gcc-10 \
    g++-10 \
    gfortran-10 \
    libgfortran5 \
    build-essential \
    make \
    cmake \
    pkg-config \
    git \
    libopenblas-dev \
    liblapack-dev \
    libopenmpi-dev \
    openmpi-bin \
    libboost-program-options-dev \
    libboost-timer-dev \
    libboost-thread-dev \
    libboost-chrono-dev \
    libboost-system-dev \
    libtirpc-dev \
    libnuma-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR ${PROJECT_ROOT}

# Trimmed Kaskade dependency tree (downloaded manually from lakeFS beforehand,
# see plan/README - only the paths actually referenced by
# installed/Makefile.Local are kept, ~2.5GB instead of the full ~25GB tree).
COPY deps/KaskadeDependencies/ ${KASKADE_ROOT}/work/input/KaskadeDependencies/
COPY deps/MKL/ ${KASKADE_ROOT}/work/input/KaskadeDependencies/MKL/

# Kaskade7 source
COPY kaskade7_test/ ${KASKADE_ROOT}/

RUN find kaskade7_test \
    -name "Makefile*" -type f \
    -exec sed -i \
    "s|/data/numerik/software/KaskadeDependencies|${KASKADE_ROOT}/work/input/KaskadeDependencies|g" \
    {} +

RUN find kaskade7_test \
    -name "Makefile*" -type f \
    -exec sed -i \
    's|CXX = .*/installed/bin/g++|CXX = /usr/bin/g++-10|g' \
    {} +

RUN find kaskade7_test \
    -name "Makefile*" -type f \
    -exec sed -i \
    's|BOOSTLIB = -L.*/installed/lib|BOOSTLIB =|g' \
    {} +

RUN mkdir -p /opt/yaml-cpp-kaskade/lib /opt/yaml-cpp-kaskade/include && \
    ln -s /usr/lib/x86_64-linux-gnu/libyaml-cpp.a /opt/yaml-cpp-kaskade/lib/libyaml-cpp.a && \
    ln -s /usr/include/yaml-cpp /opt/yaml-cpp-kaskade/include/yaml-cpp

WORKDIR ${KASKADE_ROOT}

ENV LD_LIBRARY_PATH=${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib:${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib64:${KASKADE_ROOT}/work/input/KaskadeDependencies/MKL/mkl/lib/intel64

RUN sed -i \
    's|^FLAGS = |FLAGS = -no-pie |' \
    Makefile.Local

RUN grep -RIlZ \
    -E '/root/kaskade7_master/kaskade7_test' \
    "${KASKADE_ROOT}" \
    | xargs -0 -r sed -i \
        -e "s|/root/kaskade7_master/kaskade7_test|${KASKADE_ROOT}|g" \
        -e "s|/root/covid_germany_hybrid/kaskade7_test|${KASKADE_ROOT}|g"

RUN if grep -Rns \
    -E '/root/kaskade7_master/kaskade7_test' \
    "${KASKADE_ROOT}"; then \
        echo "Ein alter Kaskade-Pfad ist noch vorhanden."; \
        exit 1; \
    fi

RUN chmod +x "${KASKADE_ROOT}/tools/gccmakedep"

RUN make kasklib || \
    { \
        [ ! -f gccerr.txt ] || cat gccerr.txt; \
        exit 1; \
    }

RUN make FILE=${PROGRAM} build-covid_germany_abm_pde_ode-tutorial || \
    { \
        [ ! -f gccerr.txt ] || cat gccerr.txt; \
        exit 1; \
    }

########################################################################
# runtime stage
########################################################################
FROM ubuntu:${UBUNTU_VERSION} AS runtime
ARG PROGRAM

ENV DEBIAN_FRONTEND=noninteractive

ENV PROJECT_ROOT=/root/covid_germany_hybrid
ENV KASKADE_ROOT=/root/covid_germany_hybrid/kaskade7_test
ENV TZ=Europe/Berlin

RUN apt-get update && apt-get install -y \
    libnuma1 \
    gnuplot-nox \
    tzdata \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone \
    && rm -rf /var/lib/apt/lists/*

# Runtime shared libraries only (.so*) - headers, static archives (.a) and the
# build toolchain from the build stage are not needed to run the compiled
# binary. Kept at the same absolute paths the binary's baked-in rpath expects.
COPY --from=build /root/covid_germany_hybrid/kaskade7_test/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib/ ${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib/
COPY --from=build /root/covid_germany_hybrid/kaskade7_test/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib64/ ${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib64/
COPY --from=build /root/covid_germany_hybrid/kaskade7_test/work/input/KaskadeDependencies/MKL/mkl/lib/intel64/ ${KASKADE_ROOT}/work/input/KaskadeDependencies/MKL/mkl/lib/intel64/

RUN find "${KASKADE_ROOT}/work/input" -type f \( -name "*.a" -o -name "*.la" \) -delete

ENV LD_LIBRARY_PATH=${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib:${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib64:${KASKADE_ROOT}/work/input/KaskadeDependencies/MKL/mkl/lib/intel64

# Compiled model binary
COPY --from=build ${KASKADE_ROOT}/tutorial/covid_germany_abm_pde_ode/${PROGRAM} ${KASKADE_ROOT}/tutorial/covid_germany_abm_pde_ode/model

# Mesh/domain data (triangulated state boundaries etc.) - non-confidential,
# part of the source tree, read by the binary relative to its working
# directory at runtime (e.g. domain/domain_Berlin.1.node).
COPY kaskade7_test/tutorial/covid_germany_abm_pde_ode/domain/ ${KASKADE_ROOT}/tutorial/covid_germany_abm_pde_ode/domain/

# Parameter file for the model binary.
COPY run-config.yaml ${PROJECT_ROOT}/run-config.yaml

# Confidential input data (work/input_data) and run output (work/output) are
# supplied at `docker run` time via bind mounts, e.g.:
#   docker run \
#     -v /local/path/to/input_data:${KASKADE_ROOT}/work/input_data \
#     -v /local/path/to/output:${KASKADE_ROOT}/work/output \
#     <image>

WORKDIR ${KASKADE_ROOT}/tutorial/covid_germany_abm_pde_ode

# work/output is a bind mount supplied at `docker run` time; the model
# expects this subdirectory to already exist.
ENV CONFIG=${PROJECT_ROOT}/run-config.yaml
ENV LOG_DIR=${KASKADE_ROOT}/work/output
CMD ["sh", "-c", "LOG_FILE=\"$LOG_DIR/run_$(date +'%Y-%m-%d_%H-%M-%S').log\"; exec ./model --config \"$CONFIG\" > \"$LOG_FILE\" 2>&1"]