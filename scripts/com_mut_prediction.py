#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage Example:
  1) Without --vcf (uses percentile threshold):
     python script.py -i my_sca_matrix.tsv -o out.txt -n 2
  2) With --vcf (uses numeric threshold):
     python script.py -i dummy.csv --vcf path/to/mydata.vcf -o out.txt -n 3

Key points:
  - If --vcf is used, we switch to a direct numeric threshold (default = 0.5 if -t is not given).
  - If --vcf and -n are both provided:
      * We parse the VCF's prefix (basename minus extension),
      * read prefix.sca_matrix.tsv from the same directory,
      * take that matrix's size => length => generate combinations of size -n.
  - This example includes map_positions(...) or merge_three_results(...) placeholders that you can adapt
    for your "three-file merging" or mapping logic.
"""

import sys
import os
import argparse
import numpy as np
import pandas as pd
from itertools import combinations

###############################################################################
# Placeholder for advanced logic: AAA->DNA->VCF mapping, merges, etc.
###############################################################################

def map_positions(aa_position_file, vcf_file, output_mapped_file):
    """
    Example placeholder for amino acid pos => DNA => VCF line mapping.
    Implementation is domain-specific. You could adapt your earlier code or
    combine logic from map_positions(...) that you posted before.
    """
    print(f"[Placeholder] map_positions => {aa_position_file} => {vcf_file} => {output_mapped_file}")
    # TODO: implement your domain logic here
    pass

def merge_three_results(file1, file2, file3, output_file):
    """
    Example placeholder for merging 3 results into a single integrated file.
    You can adapt from your earlier partial implementation. 
    """
    print(f"[Placeholder] merge_three_results => {file1}, {file2}, {file3} => {output_file}")
    # TODO: implement your domain logic
    pass

###############################################################################
# Utility functions for reading / threshold
###############################################################################

def flatten_and_sort_desc(matrix_2d):
    """
    Flatten a 2D matrix into 1D and sort in descending order.
    """
    return np.sort(matrix_2d.flatten())[::-1]

def pick_percentile_value(sorted_desc, percentile=0.9):
    """
    Return the value at the requested percentile from a descending-sorted array.
    percentile=0.9 => top 10%.
    """
    total = len(sorted_desc)
    if total < 1:
        return np.nan
    rank = int((1 - percentile) * total)
    rank = max(1, min(rank, total))
    return sorted_desc[rank - 1]

###############################################################################
# Optional parse_comb_file for "refPosAlt" formatting (if needed)
###############################################################################

def parse_comb_file(comb_file):
    """
    Parse a .comb file with columns like:
        pos   ref_aa   highest_freq_aa  frequency   SS
    Return a dictionary: comb_dict[pos] = (ref_aa, alt_aa)
    """
    df = pd.read_csv(comb_file, delim_whitespace=True)
    comb_dict = {}
    for _, row in df.iterrows():
        pos = int(row['pos'])
        ref_aa = str(row['ref_aa'])
        alt_aa = str(row['highest_freq_aa'])
        comb_dict[pos] = (ref_aa, alt_aa)
    return comb_dict

def format_positions(positions, comb_dict=None):
    """
    Convert positions [p1, p2, ...] to "refPalt" if comb_dict is available,
    otherwise just numeric strings.
    E.g., p=1 => M1A if comb_dict[1] = (M, A).
    """
    out_list = []
    for p in positions:
        if comb_dict and p in comb_dict:
            ref_aa, alt_aa = comb_dict[p]
            out_list.append(f"{ref_aa}{p}{alt_aa}")
        else:
            out_list.append(str(p))
    return ",".join(out_list)

###############################################################################
# Core calculations (LD, SCA) - placeholders
###############################################################################

def calculate_mean_ld(ld_matrix):
    n = ld_matrix.shape[0]
    if n < 2:
        return np.nan
    vals = ld_matrix[np.triu_indices(n,1)]
    return vals.mean() if len(vals) else np.nan

def calculate_multilocus_ld(ld_matrix):
    n = ld_matrix.shape[0]
    if n < 2:
        return np.nan
    vals = ld_matrix[np.triu_indices(n,1)]
    return np.prod(vals) if len(vals) else np.nan

def calculate_sca_values(ld_matrix):
    n = ld_matrix.shape[0]
    sca = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            cov = ld_matrix[i,j] - ld_matrix[i,i]*ld_matrix[j,j]
            if ld_matrix[i,i] == 0 or ld_matrix[j,j] == 0:
                w_i = 0
                w_j = 0
            else:
                w_i = (ld_matrix[i,i]*(1-ld_matrix[j,j])) / ((1-ld_matrix[i,i])*ld_matrix[j,j])
                w_j = (ld_matrix[j,j]*(1-ld_matrix[i,i])) / ((1-ld_matrix[j,j])*ld_matrix[i,i])
            sca_val = w_i * w_j * cov
            sca[i,j] = sca_val
            sca[j,i] = sca_val
    return sca

def calculate_specific_combinations(ld_matrix, zero_based_positions):
    sub_mat = ld_matrix[np.ix_(zero_based_positions, zero_based_positions)]
    if np.isnan(sub_mat).any():
        # Return all nan if submatrix is invalid
        return sub_mat, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    log_sub_mat = np.log(sub_mat + 1e-10)

    mean_ld = calculate_mean_ld(sub_mat)
    mean_ld_log = calculate_mean_ld(log_sub_mat)
    multi_ld = calculate_multilocus_ld(sub_mat)
    multi_ld_log = calculate_multilocus_ld(log_sub_mat)
    sca_mat = calculate_sca_values(sub_mat)
    sca_mat_log = calculate_sca_values(log_sub_mat)

    return (sub_mat, log_sub_mat, mean_ld, mean_ld_log,
            multi_ld, multi_ld_log, sca_mat, sca_mat_log)

###############################################################################
# process_combinations
###############################################################################

def process_combinations(ld_matrix_file,
                         output_file,
                         combination_file=None,
                         combos=None,
                         ld_threshold=None,
                         output_matrices=False,
                         comb_dict=None):
    """
    - combos: a list of integer-lists. If provided, we skip reading from combination_file.
    - Otherwise, read combos from combination_file. Each line: "1 2 5" (1-based).
    - For each combo, compute submatrix metrics. If meanLD >= threshold => label=1 else 0.
    - If comb_dict is provided, convert positions to 'refPosAlt' style instead of "pos".
    """
    ld_mat = pd.read_csv(ld_matrix_file, header=None).values
    n = ld_mat.shape[0]

    if combos is not None:
        combo_data = combos
    else:
        if combination_file is None:
            raise ValueError("[ERROR] No combos nor combination_file provided.")
        combo_data = []
        with open(combination_file, 'r') as f:
            for line in f:
                tokens = line.strip().split()
                if not tokens:
                    combo_data.append(None)
                    continue
                try:
                    arr = [int(x) for x in tokens]
                except ValueError:
                    combo_data.append(None)
                    continue
                # check range
                if any((p<1 or p>n) for p in arr):
                    combo_data.append(None)
                else:
                    combo_data.append(arr)

    with open(output_file, 'w') as out:
        for positions in combo_data:
            if positions is None:
                out.write("None\tNaN\n")
                continue

            zero_idx = [p-1 for p in positions]
            (sub_m, sub_m_log,
             mean_ld, mean_ld_log,
             multi_ld, multi_ld_log,
             sca_m, sca_m_log) = calculate_specific_combinations(ld_mat, zero_idx)

            pos_str = format_positions(positions, comb_dict=comb_dict)
            if np.isnan(mean_ld):
                out.write(f"{pos_str}\tNaN\n")
                continue

            label = 0
            if ld_threshold is not None and mean_ld >= ld_threshold:
                label = 1

            out.write(f"{pos_str}\tMean Pairwise score: {mean_ld:.4f}\tLog Mean Pairwise score: {mean_ld_log:.4f}\t{label}\n")

            if output_matrices:
                out.write("Sub Matrix:\n")
                np.savetxt(out, sub_m, fmt="%.6f")
                out.write("Log Sub Matrix:\n")
                np.savetxt(out, sub_m_log, fmt="%.6f")
                out.write(f"Multilocus LD: {multi_ld:.6f}\n")
                out.write(f"Log Multilocus LD: {multi_ld_log:.6f}\n")
                out.write("SCA Matrix:\n")
                np.savetxt(out, sca_m, fmt="%.6f")
                out.write("Log SCA Matrix:\n")
                np.savetxt(out, sca_m_log, fmt="%.6f")
                out.write("\n")

###############################################################################
# MAIN
###############################################################################

def main():
    parser = argparse.ArgumentParser(description="LD/SCA script with optional VCF and combos logic.")
    parser.add_argument("-i", "--input", required=True,
                        help="LD/SCA matrix file to parse (CSV/TSV, no header).")
    parser.add_argument("-o", "--output", required=True, help="Path to output file.")
    parser.add_argument("-l", "--list", help="File with combos of positions (1-based).")
    parser.add_argument("--vcf", help="If provided => numeric threshold mode. If -n => read prefix.sca_matrix.tsv for length.")
    parser.add_argument("--comb", help="Path to .comb file, for output in refPosAlt style.")
    parser.add_argument("-n", "--number", type=int,
                        help="Generate combos of size N from [1..protein_length]. Mutually exclusive with -l.")
    parser.add_argument("--output_matrices", action="store_true",
                        help="If set, output sub-matrices and SCA details.")
    parser.add_argument("-t", "--threshold", type=float,
                        help="If --vcf => direct numeric cutoff (default=0.5), else percentile (default=0.9).")

    args = parser.parse_args()

    # Check mutual exclusivity
    if args.number is not None and args.list is not None:
        sys.exit("[ERROR] -n/--number and -l/--list cannot both be specified.")

    # If user wants to do "three-file merges" or "pos->DNA->VCF mapping",
    # you may call 'map_positions' or 'merge_three_results' here or after combos are generated.

    # Load .comb file if provided
    comb_dict = None
    if args.comb:
        if not os.path.isfile(args.comb):
            sys.exit(f"[ERROR] .comb file not found: {args.comb}")
        comb_dict = parse_comb_file(args.comb)
        print(f"[INFO] Loaded comb_dict with {len(comb_dict)} entries.")

    use_vcf_mode = (args.vcf is not None)
    if use_vcf_mode:
        # numeric threshold
        if args.threshold is None:
            args.threshold = 0.5
        print(f"[INFO] Using direct numeric threshold => {args.threshold}")
        ld_threshold = args.threshold

        # If also using -n => parse prefix.sca_matrix.tsv to get length
        if args.number is not None:
            vcf_dir = os.path.dirname(args.vcf) or "."
            vcf_prefix = os.path.splitext(os.path.basename(args.vcf))[0]
            sca_file = os.path.join(vcf_dir, vcf_prefix + ".sca_matrix.tsv")
            if not os.path.isfile(sca_file):
                sys.exit(f"[ERROR] Could not find matching sca_matrix => {sca_file}")
            sca_data = pd.read_csv(sca_file, header=None).values
            length = sca_data.shape[0]
            print(f"[INFO] Found sca_matrix of size {length} => generating combos of size {args.number}")
            all_c = list(combinations(range(1, length+1), args.number))
            combos = [list(x) for x in all_c]
            print(f"[INFO] Generated {len(combos)} combos from sca_matrix.")
        else:
            combos = None

    else:
        # percentile mode
        if args.threshold is None:
            args.threshold = 0.9
        print(f"[INFO] Using percentile approach => {args.threshold:.2f} => top {(1-args.threshold)*100:.1f}%")

        # read main input for flatten
        main_mat = pd.read_csv(args.input, header=None).values
        sorted_desc = flatten_and_sort_desc(main_mat)
        numeric_cut = pick_percentile_value(sorted_desc, args.threshold)
        print(f"[INFO] Computed cutoff => {numeric_cut}")
        ld_threshold = numeric_cut

        if args.number is not None:
            # If user gave -n but no vcf => we use the -i file dimension
            length = main_mat.shape[0]
            print(f"[INFO] Generating combos from 1..{length}, size={args.number}")
            all_c = list(combinations(range(1, length+1), args.number))
            combos = [list(x) for x in all_c]
            print(f"[INFO] Generated {len(combos)} combos.")
        else:
            combos = None

    # Now run process_combinations
    process_combinations(
        ld_matrix_file=args.input,
        output_file=args.output,
        combination_file=args.list,
        combos=combos,
        ld_threshold=ld_threshold,
        output_matrices=args.output_matrices,
        comb_dict=comb_dict
    )

    # Optional: do final merges or mapping
    # map_positions(...)? merge_three_results(...)? etc. 
    # left as placeholders.

    print("[DONE] =>", args.output)


if __name__ == "__main__":
    main()
