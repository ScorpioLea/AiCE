#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# A wrapper pipeline to run com_mut_prediction.py twice for each .fa file in
# <input_dir>. It checks if <number_or_list> is an integer or a string to
# decide whether to call "-n <number>" or "-l <listfile>".
#
# Usage:
#   ./run_com_mut.sh <script_dir> <input_dir> <number-or-list> <output_dir>
#
# Arguments:
#   script_dir        Directory containing com_mut_prediction.py
#   input_dir         Directory containing *.fa (and their matching .sca_matrix.tsv,
#                    .ld, .comb, .vcf files)
#   number-or-list    If a numeric value, use -n <number>. If not, use -l <list>.
#   output_dir        Output directory for .result files
#
# Example:
#   ./run_com_mut.sh scripts/pySCA ./output 2 ./results
#   ./run_com_mut.sh scripts/pySCA ./output aa.list ./results
#
# This script expects files like:
#   <prefix>.fa
#   <prefix>.sca_matrix.tsv
#   <prefix>.ld
#   <prefix>.comb
#   <prefix>.vcf
# in <input_dir>.
# ------------------------------------------------------------------------------

set -euo pipefail

if [ $# -ne 4 ]; then
  echo "Usage: $0 <script_dir> <input_dir> <number-or-list> <output_dir>"
  exit 1
fi

script_dir="$1"
input_dir="$2"
param="$3"           # either an integer or a list filename
output_dir="$4"

# Create output directory if needed
mkdir -p "$output_dir"

# Check if 'param' is an integer using a simple regex
# (You could also use: if [[ $param =~ ^[0-9]+$ ]]; then ...)
if [[ "$param" =~ ^[0-9]+$ ]]; then
  is_number=true
else
  is_number=false
fi

# Loop over all .fa files in input_dir
for fa_file in "$input_dir"/*.fa; do
  # Get filename prefix (e.g. "SpCas9_af3" from "SpCas9_af3.fa")
  prefix="$(basename "$fa_file" .fa)"

  # Construct related file paths
  sca_file="$input_dir/${prefix}.sca_matrix.tsv"
  comb_file="$input_dir/${prefix}.comb"
  ld_file="$input_dir/${prefix}.ld"
  vcf_file="$input_dir/${prefix}.vcf"

  # Build output names
  sca_result="$output_dir/${prefix}.sca.result"
  ld_result="$output_dir/${prefix}.ld.result"

  echo "--------------------------------------------------"
  echo "[INFO] Processing prefix: $prefix"
  echo "      -> SCA file: $sca_file"
  echo "      -> COMB file: $comb_file"
  echo "      -> LD file: $ld_file"
  echo "      -> VCF file: $vcf_file"

  if [ "$is_number" = true ]; then
    # param is an integer => use -n <param>
    echo "[INFO] Using -n $param"

    # 1) Run com_mut_prediction.py for sca_matrix
    cmd_sca=(
      python "$script_dir/com_mut_prediction.py"
      -i "$sca_file"
      --comb "$comb_file"
      -n "$param"
      -o "$sca_result"
    )
    echo "[CMD] ${cmd_sca[*]}"
    "${cmd_sca[@]}"

    # 2) Run com_mut_prediction.py for ld (with --vcf)
    cmd_ld=(
      python "$script_dir/com_mut_prediction.py"
      -i "$ld_file"
      --comb "$comb_file"
      -n "$param"
      -o "$ld_result"
      --vcf "$vcf_file"
    )
    echo "[CMD] ${cmd_ld[*]}"
    "${cmd_ld[@]}"

  else
    # param is a string => use -l <param>
    echo "[INFO] Using -l $param"

    # 1) Run com_mut_prediction.py for sca_matrix
    cmd_sca=(
      python "$script_dir/com_mut_prediction.py"
      -i "$sca_file"
      --comb "$comb_file"
      -l "$param"
      -o "$sca_result"
    )
    echo "[CMD] ${cmd_sca[*]}"
    "${cmd_sca[@]}"

    # 2) Run com_mut_prediction.py for ld (with --vcf)
    cmd_ld=(
      python "$script_dir/com_mut_prediction.py"
      -i "$ld_file"
      --comb "$comb_file"
      -l "$param"
      -o "$ld_result"
      --vcf "$vcf_file"
    )
    echo "[CMD] ${cmd_ld[*]}"
    "${cmd_ld[@]}"
  fi

done

echo "--------------------------------------------------"
echo "[INFO] All done! Results are in $output_dir"
