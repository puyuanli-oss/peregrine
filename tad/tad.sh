#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 || -z "$1" ]]; then
    printf 'Usage: %s OUTPUT_DIRECTORY\n' "$0" >&2
    exit 1
fi

REPO_DIR="$(CDPATH= cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="$1"
if [[ "$TARGET" != /* ]]; then
    TARGET="$PWD/$TARGET"
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
    "$PYTHON" "$REPO_DIR/tad/$stage.py" "$TARGET"
    echo "$stage done"
done