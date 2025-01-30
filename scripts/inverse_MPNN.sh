#!/bin/bash

# Parameter validation
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 input_folder [output_folder]"
    echo "Example: $0 ./pdbs ../custom_output"
    exit 1
fi

# Input parameters
input_folder="$1"
output_dir="${2:-../output}"  # Default output directory

# Validate input folder existence
if [ ! -d "$input_folder" ]; then
    echo "Error: Input folder does not exist: $input_folder"
    exit 2
fi

# Create output directory (including parent directories)
mkdir -p "$output_dir" || {
    echo "Error: Failed to create output directory: $output_dir"
    exit 3
}

# Build intermediate file path
path_for_parsed_chains="${output_dir}/parsed_pdbs.jsonl"

# Step 1: Parse PDB files
python ../scripts/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
    --input_path="$input_folder" \
    --output_path="$path_for_parsed_chains" || {
    echo "Error: PDB parsing failed"
    exit 4
}

# Step 2: Run ProteinMPNN
python ../scripts/ProteinMPNN/protein_mpnn_run.py \
    --jsonl_path "$path_for_parsed_chains" \
    --out_folder "$output_dir" \
    --num_seq_per_target 1000 \
    --sampling_temp "0.5" \
    --seed 37 \
    --batch_size 16 || {
    echo "Error: ProteinMPNN execution failed"
    exit 5
}

echo "Processing completed! Results saved to: $output_dir"
