#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# This script runs a pipeline using scaProcessMSA.py and scaCore.py for each
# *.fa file found in <input_dir>, storing outputs in <output_dir>, and
# removing the *processed.fasta afterwards.
#
# Usage:
#   ./run_sca_pipeline.sh <script_dir> <input_dir> <output_dir>
#
# Example:
#   ./run_sca_pipeline.sh ./scripts/pySCA ./input_fa ./output
#
# Requirements:
#   - scaProcessMSA.py and scaCore.py must be located in <script_dir>.
#   - Python environment with scaTools etc. installed.
#   - This script must have execution permission (chmod +x).
# ------------------------------------------------------------------------------

set -euo pipefail

# 1) Parse arguments
if [ $# -lt 3 ]; then
  echo "Usage: $0 <script_dir> <input_dir> <output_dir>"
  echo "Example: $0 ./scripts/pySCA ./input_fa ./output"
  exit 1
fi

script_dir="$1"
input_dir="$2"
output_dir="$3"

# 2) Ensure output directory exists
mkdir -p "$output_dir"

# 3) Locate all .fa files in input_dir
fa_files=("$input_dir"/*.fa)
if [ ${#fa_files[@]} -eq 0 ]; then
  echo "[ERROR] No .fa files found in $input_dir"
  exit 2
fi

echo "[INFO] Found ${#fa_files[@]} .fa files in $input_dir. Results go to $output_dir."

# 4) Loop over each .fa file
for fa_file in "${fa_files[@]}"; do
    fa_basename=$(basename "$fa_file")  # e.g. "SpCas9_af3.fa"
    prefix="${fa_basename%.fa}"        # e.g. "SpCas9_af3"

    echo "-------------------------------------------------"
    echo "[INFO] Processing: $fa_file"

    # Step A: Run scaProcessMSA.py -> generates prefix.db & prefixprocessed.fasta
    python "$script_dir/scaProcessMSA.py" \
        "$fa_file" \
        -d "$output_dir" \
        -i 0

    db_file="$output_dir/${prefix}.db"
    processed_fasta="$output_dir/${prefix}processed.fasta"

    if [ ! -f "$db_file" ]; then
        echo "[ERROR] Expected DB file not found: $db_file"
        exit 3
    fi

    # Step B: Run scaCore.py on the newly created .db
    echo "[INFO] Running scaCore.py on $db_file ..."
    python "$script_dir/scaCore.py" "$db_file"

    # Step C: Remove the *processed.fasta if it exists
    if [ -f "$processed_fasta" ]; then
        rm -f "$processed_fasta"
        echo "[INFO] Removed $processed_fasta"
    fi

    echo "[INFO] Done with $prefix."
done

echo "-------------------------------------------------"
echo "[INFO] All done! Final outputs are in: $output_dir"
