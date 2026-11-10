#!/bin/bash
set -euo pipefail

usage() {
    printf 'Usage: %s OUTPUT_DIRECTORY [--chia-manifest FILE] [--chia-genes FILE] [--only-chia | --final-only]\n' "$0" >&2
}

if [[ $# -eq 0 || -z "$1" ]]; then
    usage
    exit 1
fi
TARGET="$1"
shift
CHIA_MANIFEST=""
CHIA_GENES=""
MODE="full"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --chia-manifest|--chia-genes)
            if [[ $# -lt 2 || -z "$2" ]]; then usage; exit 1; fi
            if [[ "$1" == --chia-manifest ]]; then CHIA_MANIFEST="$2"; else CHIA_GENES="$2"; fi
            shift 2
            ;;
        --only-chia|--final-only)
            if [[ "$MODE" != full ]]; then usage; exit 1; fi
            MODE="${1#--}"; shift ;;
        *) usage; exit 1 ;;
    esac
done
if [[ "$MODE" == only-chia && -z "$CHIA_MANIFEST" ]]; then
    printf "ERROR: --only-chia requires --chia-manifest\n" >&2; exit 1
fi
if [[ -n "$CHIA_GENES" && -z "$CHIA_MANIFEST" ]]; then
    printf 'ERROR: --chia-genes requires --chia-manifest\n' >&2
    exit 1
fi

for program in bedtools; do
    if ! command -v "$program" >/dev/null 2>&1; then
        printf 'ERROR: Required command not found: %s\n' "$program" >&2
        exit 1
    fi
done

REPO_DIR="$(CDPATH= cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Automatically include the verified ChIA-PET manifest in full runs.
# Preserve explicit --chia-manifest overrides and other execution modes.
if [[ "$MODE" == full && -z "$CHIA_MANIFEST" ]]; then
    CHIA_MANIFEST="$REPO_DIR/chia_pet/samples.tsv"
fi

if [[ "$TARGET" != /* ]]; then TARGET="$PWD/$TARGET"; fi
if [[ -n "$CHIA_MANIFEST" && "$CHIA_MANIFEST" != /* ]]; then CHIA_MANIFEST="$PWD/$CHIA_MANIFEST"; fi
if [[ -n "$CHIA_GENES" && "$CHIA_GENES" != /* ]]; then CHIA_GENES="$PWD/$CHIA_GENES"; fi
if [[ -z "$CHIA_GENES" ]]; then CHIA_GENES="$REPO_DIR/chia_pet/hg19genes.txt"; fi
PYTHON="$REPO_DIR/venv/bin/python"
if [[ ! -x "$PYTHON" ]]; then
    printf 'ERROR: Existing virtual environment missing: %s\n' "$PYTHON" >&2
    printf 'Set up the verified environment separately before running this workflow.\n' >&2
    exit 1
fi
mkdir -p "$TARGET"
cd "$REPO_DIR"

if [[ "$MODE" == final-only ]]; then
    if [[ -n "$CHIA_MANIFEST" ]]; then
        printf 'ERROR: --final-only accepts existing ChIA output, not --chia-manifest; use --only-chia first.\n' >&2
        exit 1
    fi
    CHIA_ARGS=()
    if [[ -f "$TARGET/chia_pet/linksDBchia" ]]; then
        CHIA_ARGS=(--chia-links "$TARGET/chia_pet/linksDBchia")
    fi
    "$PYTHON" "$REPO_DIR/tad/cutdowntad.py" "$TARGET/tad" "${CHIA_ARGS[@]}"
    "$PYTHON" "$REPO_DIR/tad/orderlinks.py" "$TARGET/tad"
    exit 0
fi

if [[ -n "$CHIA_MANIFEST" ]]; then
    "$PYTHON" "$REPO_DIR/chia_pet/ChIA_PET.py" "$TARGET/chia_pet" \
        --manifest "$CHIA_MANIFEST" --genes "$CHIA_GENES" --check-inputs
fi

if [[ "$MODE" == only-chia ]]; then
    "$PYTHON" "$REPO_DIR/chia_pet/ChIA_PET.py" "$TARGET/chia_pet" \
        --manifest "$CHIA_MANIFEST" --genes "$CHIA_GENES"
    exit 0
fi

"$PYTHON" "$REPO_DIR/eQTL/eqtl.py" "$TARGET/eQTL"
echo 'eQTL done'
"$PYTHON" "$REPO_DIR/heirarchicalTAD/heirarchicalTAD.py" "$TARGET/heirarchicalTAD"
echo 'heirarchical TAD done'
TAD_ARGS=("$TARGET/tad")
if [[ -n "$CHIA_MANIFEST" ]]; then
    "$PYTHON" "$REPO_DIR/chia_pet/ChIA_PET.py" "$TARGET/chia_pet" \
        --manifest "$CHIA_MANIFEST" --genes "$CHIA_GENES"
    echo 'ChIA-PET done'
    TAD_ARGS+=(--chia-links "$TARGET/chia_pet/linksDBchia")
fi
bash -e -o pipefail "$REPO_DIR/tad/tad.sh" "${TAD_ARGS[@]}"
echo 'TAD done'
