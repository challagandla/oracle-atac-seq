#!/usr/bin/env bash
# Activation-free ORACLE pipeline runner.
set -Eeuo pipefail

PIPELINE_NAME="oracle-atac-seq"
ENV_NAME="oracle-atac-runner"
SNAKEFILE="workflow/Snakefile"
CORES=4
DRY_RUN=0
EXTRA_ARGS=()

usage() {
    cat <<EOF
Run ${PIPELINE_NAME} with its tested Snakemake environments.

Usage:
  cp config/samples.tsv config/samples.project.tsv
  cp config/project.example.yaml config/project_overrides.yaml
  # Edit config/project_overrides.yaml and config/samples.project.tsv, then:
  bash run.sh --dry-run --configfile config/config.yaml config/project_overrides.yaml
  bash run.sh --cores 8 qc_review \
    --configfile config/config.yaml config/project_overrides.yaml
  bash run.sh --cores 8 --configfile config/config.yaml config/project_overrides.yaml

Any unrecognised options are forwarded to Snakemake.
Set MINIFORGE_HOME=/custom/path if setup.sh installed Miniforge there.
EOF
}

die() { printf '[ERROR] %s\n' "$1" >&2; exit 1; }

CONDA_BIN=""

locate_conda() {
    local candidate
    local -a candidates=("${CONDA_EXE:-}")

    if [[ -n "${MINIFORGE_HOME:-}" ]]; then
        candidates+=("${MINIFORGE_HOME%/}/bin/conda")
    fi
    candidates+=(
        "$(type -P conda 2>/dev/null || true)"
        "$HOME/miniforge3/bin/conda"
        "$HOME/miniconda3/bin/conda"
    )

    for candidate in "${candidates[@]}"; do
        if [[ -n "$candidate" && -x "$candidate" ]]; then
            CONDA_BIN="$candidate"
            return 0
        fi
    done
    return 1
}

# Expose discovery to focused tests without starting a workflow when sourced.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    return 0
fi

while (($#)); do
    case "$1" in
        --cores|-c)
            (($# >= 2)) || die "$1 requires a value"
            CORES="$2"
            shift 2
            ;;
        --dry-run|-n)
            DRY_RUN=1
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        --)
            shift
            EXTRA_ARGS+=("$@")
            break
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$ROOT_DIR/.snakemake/cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$XDG_CACHE_HOME/matplotlib}"
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
[[ -f "$SNAKEFILE" ]] || die "Missing $SNAKEFILE"

locate_conda || die "Conda is missing; run: bash setup.sh"
"$CONDA_BIN" run --name "$ENV_NAME" python -c 'import sys' >/dev/null 2>&1 || \
    die "Runner $ENV_NAME is missing; run: bash setup.sh"

SNAKEMAKE_VERSION="$("$CONDA_BIN" run --name "$ENV_NAME" snakemake --version | tail -n 1)"
SNAKEMAKE_MAJOR="${SNAKEMAKE_VERSION%%.*}"
[[ "$SNAKEMAKE_MAJOR" =~ ^[0-9]+$ ]] || die "Could not parse Snakemake version: $SNAKEMAKE_VERSION"
if ((SNAKEMAKE_MAJOR >= 8)); then
    DEPLOYMENT_ARGS=(--software-deployment-method conda)
else
    DEPLOYMENT_ARGS=(--use-conda)
fi
if ((DRY_RUN)); then
    EXTRA_ARGS=(--dry-run "${EXTRA_ARGS[@]}")
fi

printf '[INFO] %s | Snakemake %s | cores=%s\n' "$PIPELINE_NAME" "$SNAKEMAKE_VERSION" "$CORES" >&2
exec env -u R_LIBS -u R_LIBS_USER -u PYTHONPATH \
    R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null \
    "$CONDA_BIN" run --no-capture-output --name "$ENV_NAME" \
    snakemake --snakefile "$SNAKEFILE" --cores "$CORES" \
    "${DEPLOYMENT_ARGS[@]}" --rerun-incomplete --printshellcmds \
    "${EXTRA_ARGS[@]}"
