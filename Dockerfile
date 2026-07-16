FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG LAKECTL_VERSION=1.83.0

# ENV PROJECT_ROOT=/root/kaskade7_master
# ENV KASKADE_ROOT=/root/kaskade7_master/kaskade7_test
ENV PROJECT_ROOT=/root/covid_germany_abm_pde_ode
ENV KASKADE_ROOT=/root/covid_germany_abm_pde_ode/kaskade7_test
ENV LAKEFS_ENDPOINT=https://lake-episerve.zib.de/api/v1

RUN apt-get update && apt-get install -y \
    curl \
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
    ca-certificates \
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

RUN curl -fsSL \
    "https://github.com/treeverse/lakeFS/releases/download/v${LAKECTL_VERSION}/lakeFS_${LAKECTL_VERSION}_Linux_x86_64.tar.gz" \
    -o /tmp/lakefs.tar.gz \
    && tar -xzf /tmp/lakefs.tar.gz -C /tmp \
    && install -m 0755 /tmp/lakectl /usr/local/bin/lakectl \
    && rm -f /tmp/lakefs.tar.gz /tmp/lakectl \
    && lakectl --version

WORKDIR ${PROJECT_ROOT}

# Projekt kopieren
COPY kaskade7_test ./kaskade7_test
# COPY run.sh /run.sh

# # ohne diese doppelung funktioniert es nicht
# COPY "${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib/libboost_program_options.so.1.73.0" \
# "${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib/libboost_program_options.so.1.73.0"

RUN --mount=type=secret,id=lakefs_access_key,required=true \
    --mount=type=secret,id=lakefs_secret_key,required=true \
    if [ ! -d "${KASKADE_ROOT}/work/input" ] || \
       [ ! -d "${KASKADE_ROOT}/work/output" ]; then \
        echo "work-Ordner fehlt oder ist unvollständig; lade ihn aus lakeFS herunter."; \
        rm -rf "${KASKADE_ROOT}/work"; \
        mkdir -p "${KASKADE_ROOT}/work"; \
        export LAKECTL_CREDENTIALS_ACCESS_KEY_ID="$(cat /run/secrets/lakefs_access_key)"; \
        export LAKECTL_CREDENTIALS_SECRET_ACCESS_KEY="$(cat /run/secrets/lakefs_secret_key)"; \
        export LAKECTL_SERVER_ENDPOINT_URL="${LAKEFS_ENDPOINT}"; \
        lakectl \
            fs download \
            lakefs://sandbox/main/RAW/work/ \
            "${KASKADE_ROOT}/work" \
            --recursive \
            --no-progress \
            || { echo "lakectl download FAILED"; exit 1; }; \
    else \
        echo "work/input und work/output sind bereits vorhanden."; \
    fi

RUN test -d "${KASKADE_ROOT}/work/input" \
    && test -d "${KASKADE_ROOT}/work/output" \
    || { echo "Der work-Ordner wurde nicht korrekt heruntergeladen."; exit 1; }

# RUN mkdir -p "${KASKADE_ROOT}/work/output"
RUN mkdir -p "${KASKADE_ROOT}/work/output/graph_abm_pde_ode"

# Dependencies kopieren
# RUN mkdir -p "${KASKADE_ROOT}/work/input/KaskadeDependencies"

# COPY "${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2" \
#      "${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2"

# # COPY KaskadeDependencies/Kaskade7.5Dependencies-10.2 \
# #      "${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2"

# COPY "${KASKADE_ROOT}/work/input/KaskadeDependencies/MKL" \
#      "${KASKADE_ROOT}/work/input/KaskadeDependencies/MKL"

RUN grep -R "/data/numerik/software/KaskadeDependencies" -n kaskade7_test || true

# RUN find kaskade7_test -type f \
#     -exec sed -i \
#     's|/data/numerik/software/KaskadeDependencies|${KASKADE_ROOT}/work/input/KaskadeDependencies|g' \
#     {} +
RUN find kaskade7_test \
    -name "Makefile*" -type f \
    -exec sed -i \
    "s|/data/numerik/software/KaskadeDependencies|${KASKADE_ROOT}/work/input/KaskadeDependencies|g" \
    {} +

# # GCC 10 aus Ubuntu verwenden
# RUN find kaskade7_test -type f \
#     -exec sed -i \
#     's|CXX = .*/installed/bin/g++|CXX = /usr/bin/g++-10|g' \
#     {} +
RUN find kaskade7_test \
    -name "Makefile*" -type f \
    -exec sed -i \
    's|CXX = .*/installed/bin/g++|CXX = /usr/bin/g++-10|g' \
    {} +

# RUN find kaskade7_test -type f \
#     -exec sed -i \
#     's|BOOSTLIB = -L.*/installed/lib|BOOSTLIB = |g' \
#     {} +
RUN find kaskade7_test \
    -name "Makefile*" -type f \
    -exec sed -i \
    's|BOOSTLIB = -L.*/installed/lib|BOOSTLIB =|g' \
    {} +

RUN grep -R "CXX =" -n "${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/Makefile.Local"

# RUN find kaskade7_test -type f \
#     -exec sed -i \
#     's|CXX = .*/installed/bin/g++|CXX = /usr/bin/g++|g' \
#     {} +
    
RUN grep -R "/data/numerik/software/KaskadeDependencies" -n kaskade7_test || true

# # Pfade anpassen
# RUN find kaskade7_test \
#     -type f \
#     -exec sed -i \
#     's|/data/numerik/software/KaskadeDependencies|${KASKADE_ROOT}/work/input/KaskadeDependencies|g' \
#     {} +

RUN mkdir -p /opt/yaml-cpp-kaskade/lib /opt/yaml-cpp-kaskade/include && \
    ln -s /usr/lib/x86_64-linux-gnu/libyaml-cpp.a /opt/yaml-cpp-kaskade/lib/libyaml-cpp.a && \
    ln -s /usr/include/yaml-cpp /opt/yaml-cpp-kaskade/include/yaml-cpp

WORKDIR ${KASKADE_ROOT}

ENV LD_LIBRARY_PATH=${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib:${KASKADE_ROOT}/work/input/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib64:${KASKADE_ROOT}/work/input/KaskadeDependencies/MKL/mkl/lib/intel64

RUN /usr/bin/g++-10 --version

RUN sed -i \
    's|^FLAGS = |FLAGS = -no-pie |' \
    Makefile.Local

# RUN chmod +x /run.sh 

# ENTRYPOINT ["/run.sh"]

# CMD ["make", "clean"]
RUN find . -name "*.o" -delete && \
    find . -name "*.d" -delete && \
    rm -f libs/libkaskade.a

RUN make kasklib || (cat gccerr.txt && false)

CMD ["make", "tutorial"]