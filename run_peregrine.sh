#!/usr/bin/env bash
# Run the existing PEREGRINE stages with explicit ChIA-PET integration.
set -euo pipefail

usage() {
    printf 'Usage: %s OUTPUT_DIRECTORY [--only-chia | --final-only] [--check-inputs] [--chia-manifest FILE] [--chia-genes FILE]\n' "$0" >&2
}

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

if [[ $# -eq 0 || -z "$1" || "$1" == --* ]]; then
    usage
    exit 1
fi

TARGET="$1"
shift
MODE="full"
MODE_SEEN=0
CHECK_INPUTS=0
MANIFEST_SEEN=0
GENES_SEEN=0
CHIA_MANIFEST=""
CHIA_GENES=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --only-chia|--final-only)
            if [[ "$MODE_SEEN" -ne 0 ]]; then
                usage
                fail 'Supply only one mode option, once.'
            fi
            MODE="${1#--}"
            MODE_SEEN=1
            shift
            ;;
        --check-inputs)
            if [[ "$CHECK_INPUTS" -ne 0 ]]; then
                usage
                fail 'Duplicate --check-inputs option.'
            fi
            CHECK_INPUTS=1
            shift
            ;;
        --chia-manifest|--chia-genes)
            if [[ $# -lt 2 || -z "$2" || "$2" == --* ]]; then
                usage
                fail "Missing file path for $1."
            fi
            if [[ "$1" == --chia-manifest ]]; then
                [[ "$MANIFEST_SEEN" -eq 0 ]] || fail 'Duplicate --chia-manifest option.'
                CHIA_MANIFEST="$2"
                MANIFEST_SEEN=1
            else
                [[ "$GENES_SEEN" -eq 0 ]] || fail 'Duplicate --chia-genes option.'
                CHIA_GENES="$2"
                GENES_SEEN=1
            fi
            shift 2
            ;;
        *)
            usage
            fail "Unknown argument: $1"
            ;;
    esac
done

if [[ "$MODE" == final-only && ( "$MANIFEST_SEEN" -ne 0 || "$GENES_SEEN" -ne 0 ) ]]; then
    fail '--final-only consumes existing output/chia_pet/linksDBchia; omit --chia-manifest and --chia-genes.'
fi

# Caller-supplied relative paths must be resolved before changing directories.
if [[ "$TARGET" != /* ]]; then TARGET="$PWD/$TARGET"; fi
if [[ -n "$CHIA_MANIFEST" && "$CHIA_MANIFEST" != /* ]]; then CHIA_MANIFEST="$PWD/$CHIA_MANIFEST"; fi
if [[ -n "$CHIA_GENES" && "$CHIA_GENES" != /* ]]; then CHIA_GENES="$PWD/$CHIA_GENES"; fi

REPO_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -z "$CHIA_MANIFEST" ]]; then CHIA_MANIFEST="$REPO_DIR/chia_pet/samples.tsv"; fi
if [[ -z "$CHIA_GENES" ]]; then CHIA_GENES="$REPO_DIR/chia_pet/hg19genes.txt"; fi

PYTHON="$REPO_DIR/venv/bin/python"
if [[ ! -x "$PYTHON" ]]; then
    if [[ "$MODE" == full ]]; then
        fail "Full runs require the existing environment: $PYTHON"
    fi
    PYTHON="$(command -v python3)" || fail 'Python 3 is required.'
fi

case "$MODE" in
    full) SCOPE="full" ;;
    only-chia) SCOPE="chia" ;;
    final-only) SCOPE="final" ;;
esac

VERIFY_ARGS=(--repo "$REPO_DIR" --scope "$SCOPE" --output "$TARGET")
if [[ "$MODE" != final-only ]]; then
    VERIFY_ARGS+=(--manifest "$CHIA_MANIFEST" --genes "$CHIA_GENES")
fi
# Every execution validates its real inputs; the wrapper's flag means check and exit.
VERIFY_ARGS+=(--check-inputs)

cd -- "$REPO_DIR"
"$PYTHON" "$REPO_DIR/verify_peregrine.py" "${VERIFY_ARGS[@]}"
if [[ "$CHECK_INPUTS" -eq 1 ]]; then exit 0; fi

case "$MODE" in
    full)
        bash "$REPO_DIR/steps.sh" "$TARGET" \
            --chia-manifest "$CHIA_MANIFEST" --chia-genes "$CHIA_GENES"
        ;;
    only-chia)
        "$PYTHON" "$REPO_DIR/chia_pet/ChIA_PET.py" "$TARGET/chia_pet" \
            --manifest "$CHIA_MANIFEST" --genes "$CHIA_GENES"
        printf 'ChIA-PET done\n'
        ;;
    final-only)
        "$PYTHON" "$REPO_DIR/tad/cutdowntad.py" "$TARGET/tad" \
            --chia-links "$TARGET/chia_pet/linksDBchia"
        printf 'cutdowntad done\n'
        "$PYTHON" "$REPO_DIR/tad/orderlinks.py" "$TARGET/tad"
        printf 'orderlinks done\n'
        printf 'TAD done\n'
        ;;
esac
