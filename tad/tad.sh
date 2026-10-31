#!/bin/bash
set -euo pipefail

if [[ $# -eq 0 || -z "$1" ]]; then
    printf 'Usage: %s OUTPUT_DIRECTORY [--chia-links FILE]\n' "$0" >&2
    exit 1
fi
TARGET="$1"
shift
CHIA_LINKS=""
if [[ $# -gt 0 ]]; then
    if [[ $# -ne 2 || "$1" != --chia-links || -z "$2" ]]; then
        printf 'Usage: %s OUTPUT_DIRECTORY [--chia-links FILE]\n' "$0" >&2
        exit 1
    fi
    CHIA_LINKS="$2"
fi

REPO_DIR="$(CDPATH= cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ "$TARGET" != /* ]]; then TARGET="$PWD/$TARGET"; fi
if [[ -n "$CHIA_LINKS" && "$CHIA_LINKS" != /* ]]; then CHIA_LINKS="$PWD/$CHIA_LINKS"; fi
if [[ -n "$CHIA_LINKS" && ! -f "$CHIA_LINKS" ]]; then
    printf 'ERROR: ChIA links file not found: %s\n' "$CHIA_LINKS" >&2
    exit 1
fi

PYTHON="$REPO_DIR/venv/bin/python"
if [[ ! -x "$PYTHON" ]]; then
    printf 'ERROR: Virtual environment Python not found: %s\n' "$PYTHON" >&2
    exit 1
fi
cd "$REPO_DIR"

for stage in \
    TADorder \
    bedtoolIntersect \
    split \
    bTADlinks \
    concatenate_linksbTAD \
    tTADlinks \
    linksbTAD_tissueReplace \
    final_concatenate \
    cutdowntad \
    orderlinks
do
    STAGE_ARGS=("$TARGET")
    if [[ "$stage" == cutdowntad && -n "$CHIA_LINKS" ]]; then
        STAGE_ARGS+=(--chia-links "$CHIA_LINKS")
    fi
    "$PYTHON" "$REPO_DIR/tad/$stage.py" "${STAGE_ARGS[@]}"
    echo "$stage done"
done
