FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG LAKECTL_VERSION=1.83.0

ENV PROJECT_ROOT=/root/covid_germany_hybrid
ENV KASKADE_ROOT=/root/covid_germany_hybrid/kaskade7_test

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
COPY kaskade7_test/ ${KASKADE_ROOT}/

RUN --mount=type=secret,id=lakefs_access_key,required=true \
    --mount=type=secret,id=lakefs_secret_key,required=true \
    --mount=type=secret,id=lakefs_endpoint,required=true \
    if [ ! -d "${KASKADE_ROOT}/work/input" ] || \
       [ ! -d "${KASKADE_ROOT}/work/output" ]; then \
        echo "work-Ordner fehlt oder ist unvollständig; lade ihn aus lakeFS herunter."; \
        rm -rf "${KASKADE_ROOT}/work"; \
        mkdir -p "${KASKADE_ROOT}/work"; \
        export LAKECTL_CREDENTIALS_ACCESS_KEY_ID="$(cat /run/secrets/lakefs_access_key)"; \
        export LAKECTL_CREDENTIALS_SECRET_ACCESS_KEY="$(cat /run/secrets/lakefs_secret_key)"; \
        export LAKECTL_SERVER_ENDPOINT_URL="$(tr -d '\r\n' < /run/secrets/lakefs_endpoint)"; \
        mkdir -p "${KASKADE_ROOT}/work/input" "${KASKADE_ROOT}/work/output"; \
        lakectl \
            fs download \
            lakefs://sandbox/main/RAW/work/input/ \
            "${KASKADE_ROOT}/work/input" \
            --recursive \
            --no-progress \
            || { echo "lakectl download FAILED (input)"; exit 1; }; \
        lakectl \
            fs download \
            lakefs://sandbox/main/RAW/work/output/ \
            "${KASKADE_ROOT}/work/output" \
            --recursive \
            --no-progress \
            || { echo "lakectl download FAILED (output)"; exit 1; }; \
    else \
        echo "work/input und work/output sind bereits vorhanden."; \
    fi

RUN test -d "${KASKADE_ROOT}/work/input" \
    && test -d "${KASKADE_ROOT}/work/output" \
    || { echo "Der work-Ordner wurde nicht korrekt heruntergeladen."; exit 1; }

RUN mkdir -p "${KASKADE_ROOT}/work/output/graph_abm_pde_ode"

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

RUN find . -name "*.o" -delete && \
    find . -name "*.d" -delete && \
    rm -f libs/libkaskade.a

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

CMD ["make", "tutorial"]