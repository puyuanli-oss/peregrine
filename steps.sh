#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 || -z "$1" ]]; then
    printf 'Usage: %s OUTPUT_DIRECTORY\n' "$0" >&2
    exit 1
fi

for program in python3 bedtools; do
    if ! command -v "$program" >/dev/null 2>&1; then
        printf 'ERROR: Required command not found: %s\n' "$program" >&2
        exit 1
    fi
done

REPO_DIR="$(CDPATH= cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET="$1"
if [[ "$TARGET" != /* ]]; then
    TARGET="$PWD/$TARGET"
fi
mkdir -p "$TARGET"
cd "$REPO_DIR"

python3 -m venv "$REPO_DIR/venv"
source "$REPO_DIR/venv/bin/activate"
PYTHON="$REPO_DIR/venv/bin/python"
"$PYTHON" -m pip install "pip==25.3"
"$PYTHON" -m pip install \
    --build-constraint "$REPO_DIR/build-constraints.txt" \
    -r "$REPO_DIR/requirements.txt"
echo 'venv done'

"$PYTHON" "$REPO_DIR/eQTL/eqtl.py" "$TARGET/eQTL"
echo 'eQTL done'

"$PYTHON" "$REPO_DIR/heirarchicalTAD/heirarchicalTAD.py" "$TARGET/heirarchicalTAD"
echo 'heirarchical TAD done'

bash -e -o pipefail "$REPO_DIR/tad/tad.sh" "$TARGET/tad"
echo 'TAD done'

# ChIA-PET remains disabled pending workflow review.
deactivate
