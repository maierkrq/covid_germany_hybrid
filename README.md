# covid_germany_abm_pde_ode Docker image

This repo builds a Docker image for the Kaskade7-based `covid_germany_abm_pde_ode`
(Berlin) model. The image is built once by a maintainer and pushed to a
registry; anyone who wants to run the model pulls it and supplies the
(confidential) model input data via a local mount at run time.

The `Dockerfile` is a two-stage build:

- **`build`** — full compiler toolchain (gcc-10, cmake, Boost, etc.) plus the
  Kaskade dependency tree, compiles the `covid` model binary.
- **`runtime`** — only the compiled binary and the shared libraries it
  actually needs (~3GB). No compiler, no headers, no source, no confidential
  data.

The dependency download is a manual step a maintainer runs locally, and
confidential model input/output are never baked into the image at all.

## For maintainers: building and publishing the image

### 1. Download the Kaskade dependencies (one-time, or when they change)

The build only needs a small subset of `lakefs://sandbox/main/RAW/work/input/`
— the parts actually referenced by
`kaskade7_test/Makefile.Local` → `installed/Makefile.Local` (~2.5GB out of the
full ~25GB tree). Requires `lakectl` configured with lakeFS credentials
(`~/.lakectl.yaml` or equivalent).

```bash
cd model_kristina

BASE="lakefs://sandbox/main/RAW/work/input/KaskadeDependencies"
DEST="deps"

lakectl fs download "$BASE/Kaskade7.5Dependencies-10.2/installed/Makefile.Local" \
  "$DEST/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/"
lakectl fs download "$BASE/Kaskade7.5Dependencies-10.2/installed/lib/" \
  "$DEST/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib/" --recursive
lakectl fs download "$BASE/Kaskade7.5Dependencies-10.2/installed/lib64/" \
  "$DEST/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/lib64/" --recursive
lakectl fs download "$BASE/Kaskade7.5Dependencies-10.2/installed/include/" \
  "$DEST/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/include/" --recursive
lakectl fs download "$BASE/MKL/mkl/lib/intel64/" \
  "$DEST/MKL/mkl/lib/intel64/" --recursive
lakectl fs download "$BASE/MKL/mkl/include/" \
  "$DEST/MKL/mkl/include/" --recursive
```

This populates `./deps/` (gitignored — never committed, ~2.5GB). You only
need to redo this when the Kaskade dependency versions actually change, not
on every build.

> **Large-file downloads may fail with `stream error: ... CANCEL`.** Some of
> the larger objects (e.g. `libmkl_core.so`, ~664MB) can hit HTTP/2
> stream-reset errors at `lakectl`'s default parallelism. If a download
> fails this way, retry with `-p 1` (single-threaded):
> `lakectl fs download ... --recursive -p 1`.

### 2. Build the image

```bash
docker build -t ghcr.io/the-episerve-consortium/kaskade-covid-berlin:<tag> .
```

`deps/` and `kaskade7_test/` are the only build inputs; no secrets, no
network access to lakeFS needed at this step.

### 3. Push to the registry

```bash
docker push ghcr.io/the-episerve-consortium/kaskade-covid-berlin:<tag>
```

Requires `docker login ghcr.io` with an account/token that has write access
to the `the-episerve-consortium` org.

## For end users: running the model

### 1. Pull the image

```bash
docker pull ghcr.io/the-episerve-consortium/kaskade-covid-berlin:<tag>
```

Available tags/versions:
https://github.com/orgs/the-episerve-consortium/packages/container/package/kaskade-covid-berlin

### 2. Get the model input data

The model needs `work/input_data` (confidential — e.g. `Berlin_V_*.bin`,
case data by state) available locally. This is **not** included in the
image; obtain it separately (e.g. from lakeFS at
`lakefs://sandbox/main/RAW/work/input_data/`) and note its local path.

```bash
lakectl fs download "lakefs://sandbox/main/RAW/work/input_data/" \
  /local/path/to/input_data --recursive
```

If this fails with a `stream error: ... CANCEL` on a large file, retry with
`-p 1` — see the note in the maintainer section above.

### 3. Run it

```bash
docker run --rm \
  -v /local/path/to/input_data:/root/covid_germany_hybrid/kaskade7_test/work/input_data \
  -v /local/path/to/output:/root/covid_germany_hybrid/kaskade7_test/work/output \
  ghcr.io/the-episerve-consortium/kaskade-covid-berlin:<tag>
```

- `input_data` mount: read-only is fine, e.g. add `:ro` to that `-v` flag.
- `output` mount: needs to be writable; results land under
  `output/output_data_optim_<N>/...` on the host.
