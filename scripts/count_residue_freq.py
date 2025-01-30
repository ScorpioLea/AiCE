#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse

def count_and_frequency(input_file, output_file):
    """
    Read a given MSA file (input_file),
    count the amino acid frequencies at each position for all sequences
    except the first one (assumed to be the reference),
    and write the results to output_file.
    """
    try:
        with open(input_file, 'r') as file:
            lines = file.readlines()

        # Parse MSA
        sequences = []
        current_seq = ''
        for line in lines:
            if line.startswith('>'):
                # Encounter a new sequence header; store the previous sequence
                if current_seq:
                    sequences.append(current_seq.upper())
                current_seq = ''
            else:
                current_seq += line.strip()
        # At the end of file, don't forget to store the last sequence
        if current_seq:
            sequences.append(current_seq.upper())

        # Ensure at least 2 sequences (1 reference + 1 other)
        if len(sequences) < 2:
            raise ValueError("The MSA file must contain at least 2 sequences (including 1 reference).")

        # Common amino acids plus the gap character
        amino_acids = "ACDEFGHIKLMNPQRSTVWY-X"
        # Reference sequence length
        ref_len = len(sequences[0])

        # Initialize counts and frequencies for each amino acid at each position
        counts = {aa: [0] * ref_len for aa in amino_acids}
        frequencies = {aa: [0] * ref_len for aa in amino_acids}
        highest_freq_amino_acids = [''] * ref_len

        # Start counting from the 2nd sequence (skip the reference)
        for seq_index in range(1, len(sequences)):
            seq = sequences[seq_index]
            # If a sequence is shorter, only take up to ref_len
            for pos in range(ref_len):
                aa = seq[pos] if pos < len(seq) else '-'
                if aa in amino_acids:
                    counts[aa][pos] += 1

        # Compute frequencies and find the highest-frequency amino acid
        for pos in range(ref_len):
            total = sum(counts[aa][pos] for aa in amino_acids)
            highest_freq = 0
            highest_freq_aa = ''
            for aa in amino_acids:
                if total > 0:
                    freq = counts[aa][pos] / total
                    frequencies[aa][pos] = freq
                    if freq > highest_freq:
                        highest_freq = freq
                        highest_freq_aa = aa
            highest_freq_amino_acids[pos] = highest_freq_aa

        # Write statistics to output_file
        with open(output_file, 'w') as f:
            # 1) Write the first row: position labels (the reference sequence characters)
            f.write("Position: " + "\t".join(sequences[0]) + "\n")

            # 2) Write raw counts
            for aa in amino_acids:
                f.write(aa + ": " + "\t".join(str(count) for count in counts[aa]) + "\n")

            # 3) Write frequencies
            for aa in amino_acids:
                f.write(aa + "_freq: " + "\t".join(f"{freq:.2%}" for freq in frequencies[aa]) + "\n")

            # 4) Write the highest-frequency amino acids
            f.write("Highest_freq_AA: " + "\t".join(highest_freq_amino_acids) + "\n")

    except FileNotFoundError:
        raise FileNotFoundError(f"Cannot find file: {input_file}")
    except Exception as e:
        raise RuntimeError(f"Error processing MSA file: {e}")


def transpose_and_pickmax(count_file, out_file):
    """
    Simulate the original AWK logic:
      1) Transpose count_file
      2) From the transposed lines, look at columns 23..43 (after removing '%')
         and find the maximum. Then output column 1, column 44, and that maximum.
    """
    # Step 1: Read count_file into rows
    with open(count_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # Split each line into tokens
    rows = [line.split() for line in lines]

    # If the file is empty or has only blank lines, just create an empty out_file
    if not rows:
        with open(out_file, 'w') as fw:
            fw.write("")
        return

    # Transpose: rows => columns
    max_cols = max(len(r) for r in rows)
    columns = [[] for _ in range(max_cols)]
    for row in rows:
        for col_index, val in enumerate(row):
            columns[col_index].append(val)

    transposed_lines = []
    for col_vals in columns:
        transposed_lines.append("\t".join(col_vals))

    # Step 2: For each transposed line:
    #    * split into tokens
    #    * parse tokens[22..42] => find max
    #    * output tokens[0], tokens[45], max + "%"
    result_lines = []
    for line in transposed_lines:
        tokens = line.split("\t")

        # We need tokens[0], tokens[45], and the max from tokens[23..45]
        # Make sure the line has enough columns
        if len(tokens) < 46:
            # If not enough columns, skip or raise an error
            continue

        # Find max in columns [23..45]; remove trailing '%' if present
        max_val = 0
        for i in range(23, 46):  # 22..42 in AWK-based indexing
            val_str = tokens[i].rstrip('%')
            try:
                val_num = float(val_str)
            except:
                val_num = 0
            if val_num > max_val:
                max_val = val_num

        # Construct the result line: "$1 \t $46 \t max_val%"
        # Where $1 => tokens[0], $46 => tokens[45]
        result_line = f"{tokens[0]}\t{tokens[45]}\t{max_val}%"
        result_lines.append(result_line)

    with open(out_file, 'w') as fw:
        for line in result_lines:
            fw.write(line + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="A Python script to parse, count frequencies, transpose, and find max, replacing AWK usage."
    )
    parser.add_argument("--msa", required=True, help="Input MSA file")
    parser.add_argument("--count", required=True, help="Intermediate statistic file (formerly 'out')")
    parser.add_argument("--out", required=True, help="Final output file (formerly 'out.freq')")
    args = parser.parse_args()

    # Step 1: count and frequency from MSA, then write to args.count
    count_and_frequency(args.msa, args.count)

    # Step 2: transpose the intermediate file, pick max, and write to args.out
    transpose_and_pickmax(args.count, args.out)

    print(f"[INFO] Done! freq_file -> {args.count}, max_freq -> {args.out}")


if __name__ == "__main__":
    main()
