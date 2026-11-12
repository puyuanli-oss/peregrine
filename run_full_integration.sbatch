#!/bin/bash
#SBATCH --job-name=peregrine_full_chia
#SBATCH --partition=main
#SBATCH --output=logs/full_integration_%j.out
#SBATCH --error=logs/full_integration_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT="${PEREGRINE_OUTPUT:-/scratch1/${USER}/peregrine_full_${SLURM_JOB_ID}}"

cd "$ROOT"

module load gcc/13.3.0
module load bedtools2/2.31.1

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "=== START ==="
date
echo "Commit: $(git rev-parse HEAD)"
echo "Output: $OUTPUT"

if [ -e "$OUTPUT" ]; then
    echo "ERROR: Output directory already exists: $OUTPUT" >&2
    exit 1
fi

./run_peregrine.sh "$OUTPUT" --check-inputs
./run_peregrine.sh "$OUTPUT"

echo "=== FULL PEREGRINE RUN COMPLETED ==="
date
