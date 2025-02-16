#!/bin/bash

# ------------------------------------------------------------------------------
# This script automates a multi-step protein analysis pipeline.
# It first invokes 'inverse_MPNN.sh' to perform inverse folding on all PDB files
# in the input folder, generating .seq/.fa files. Then, for each .fa file,
# it runs an awk command to reformat FASTA headers,
# then does secondary structure prediction, residue frequency calculation,
# and mutation prediction (auto or manual beta/gamma) via predicted_single_HF_mutations.py.
#
# Usage:
#   ./run_pipeline.sh <scripts_dir> <input_folder> [output_folder]
#
# Example:
#   ./run_pipeline.sh ../scripts ./input_pdbs ../custom_output
# ------------------------------------------------------------------------------

# 1. Parameter check
if [ $# -lt 2 ]; then
    echo "Usage: $0 <scripts_dir> <input_folder> [output_folder]"
    echo "Example: $0 ../scripts ./input_pdbs ../custom_output"
    exit 1
fi

# 2. Parse input arguments
scripts_dir="$1"
input_dir="$2"
output_dir="${3:-../output}"

# 3. Validate scripts directory
if [ ! -d "$scripts_dir" ]; then
    echo "Error: Scripts directory not found: $scripts_dir"
    exit 2
fi

# 4. Create output directory if it does not exist
mkdir -p "$output_dir" || {
    echo "Error: Failed to create output directory: $output_dir"
    exit 3
}

# 5. Run inverse folding (MPNN) on all PDB files
echo "Running inverse folding (MPNN) on all PDB files in: $input_dir ..."
sh "$scripts_dir/inverse_MPNN.sh" "$input_dir" "$output_dir" || {
    echo "Error: Inverse folding failed"
    exit 4
}

# 6. Locate generated *.seq or *.fa files
seq_dir="$output_dir/seqs"
if [ ! -d "$seq_dir" ]; then
    echo "Error: No 'seq' directory found under $output_dir. Check MPNN output."
    exit 5
fi

shopt -s nullglob
seq_files=("$seq_dir"/*.fa)
shopt -u nullglob

if [ ${#seq_files[@]} -eq 0 ]; then
    echo "Error: No *.fa files found in $seq_dir. Check MPNN output."
    exit 6
fi

echo "Found ${#seq_files[@]} .fa files. Proceeding with reformat and subsequent predictions..."

# 7. Process each *.fa file
for seq_file in "${seq_files[@]}"; do
    file_prefix=$(basename "$seq_file" .fa)
    echo "-----------------------------------------"
    echo "[INFO] Processing sequence file: $file_prefix.fa"

    # Move the .fa file to the main output directory
    mv "$seq_file" "$output_dir/${file_prefix}.fa" 2>/dev/null

    # --------------------------------------------------------------------------
    # Reformat FASTA headers (awk) to make the first ID >ref, others >1, >2, ...
    # --------------------------------------------------------------------------
    tmp_file="$output_dir/${file_prefix}.tmp.fa"
    awk 'BEGIN{c=0}; /^>/{c++; if(c==1){print ">ref"} else {print ">" c-1}; next}; {print}' \
        "$output_dir/${file_prefix}.fa" > "$tmp_file"

    # Replace original .fa with reformatted version
    mv "$tmp_file" "$output_dir/${file_prefix}.fa"
    echo "  [Done] FASTA header reformatted -> $file_prefix.fa"

    # Identify corresponding PDB (assuming same prefix)
    pdb_file="$input_dir/${file_prefix}.pdb"
    if [ ! -f "$pdb_file" ]; then
        echo "Warning: No corresponding PDB file found for prefix '$file_prefix' in $input_dir."
        echo "Skipping secondary structure prediction for $file_prefix."
        continue
    fi

    # Step A: Predict secondary structure
    echo "[Step A] Predicting secondary structure for $file_prefix..."
    python "$scripts_dir/predict_dssp.py" \
        "$pdb_file" \
        "$output_dir/${file_prefix}.ss" || {
        echo "Error: DSSP prediction failed for $file_prefix"
        exit 7
    }

    # Step B: Calculate residue frequencies
    echo "[Step B] Calculating residue frequencies for $file_prefix..."
    python "$scripts_dir/count_residue_freq.py" \
        --msa "$output_dir/${file_prefix}.fa" \
        --count "$output_dir/${file_prefix}.freq" \
        --out "$output_dir/${file_prefix}.re.freq" || {
        echo "Error: Frequency calculation failed for $file_prefix"
        exit 8
    }

    # Step C: Generate mutation predictions
    echo "[Step C] Generating mutation predictions for $file_prefix..."

    # If you want purely automatic (model-based) thresholds for each file:
    python "$scripts_dir/predicted_single_HF_mutations.py" \
        -freq "$output_dir/${file_prefix}.re.freq" \
        -dssp "$output_dir/${file_prefix}.ss" \
        -comb "$output_dir/${file_prefix}.comb" \
        -mut "$output_dir/${file_prefix}.mut" \
        -model_path "$scripts_dir" || {
        echo "Error: Mutation prediction failed for $file_prefix"
        exit 9
    }

    # Alternatively, if you want to specify certain manual thresholds,
    # you could do:
    # python "$scripts_dir/predicted_single_HF_mutations.py" \
    #    -freq "$output_dir/${file_prefix}.re.freq" \
    #    -dssp "$output_dir/${file_prefix}.ss" \
    #    -beta 0.3 \
    #    -gama 0.2 \
    #    -comb "$output_dir/${file_prefix}.comb" \
    #    -mut "$output_dir/${file_prefix}.mut"

    echo "Done with $file_prefix!"
done

# 8. Optional: Cleanup any intermediate files if needed
rm -f "$output_dir/parsed_pdbs.jsonl"
rm -rf "$seq_dir"

echo "-----------------------------------------"
echo "Pipeline completed successfully!"
echo "All output files are located in: $output_dir"
ls -lh "$output_dir"
