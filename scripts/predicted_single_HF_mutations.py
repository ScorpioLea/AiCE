#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def parse_out_freq(freq_file):
    """
    Parse out.freq file.
    Each (non-header) line has 3 columns:
        ref_aa, highest_freq_aa, freq%
    If highest_freq_aa == '-', it means a gap (kept as '-').
    Return list of (ref_aa, highest_freq_aa, freq_float).
    """
    freq_data = []
    with open(freq_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # Skip header if it contains "Position:"
    idx = 0
    if lines and "Position:" in lines[0]:
        idx = 1

    for line in lines[idx:]:
        parts = line.split()
        if len(parts) < 3:
            continue

        ref_aa = parts[0]
        highest_freq_aa = parts[1]  # could be '-', or actual AA
        freq_str = parts[2].rstrip('%')  # e.g. "67.34%"

        try:
            freq_val = float(freq_str) / 100.0  # Convert percentage -> 0..1
        except ValueError:
            freq_val = 0.0

        freq_data.append((ref_aa, highest_freq_aa, freq_val))

    return freq_data


def parse_out_txt(dssp_file):
    """
    Parse out.txt (DSSP results).
    Each line: pos   ref_aa   SS
    Return list of (pos_int, ref_aa, SS).
    """
    dssp_data = []
    with open(dssp_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                pos = int(parts[0])
            except ValueError:
                continue
            ref_aa = parts[1]
            ss = parts[2]
            dssp_data.append((pos, ref_aa, ss))
    return dssp_data


def merge_data(freq_data, dssp_data):
    """
    Merge freq_data & dssp_data, taking into account 'X' lines in freq_data.

    We'll keep two indices:
      i for freq_data
      j for dssp_data

    We'll also maintain an output pos_counter that increments for every line of freq_data.

    - If freq_data[i].ref_aa == 'X':
        -> This line does NOT consume dssp_data[j].
        -> We'll produce (pos_counter, 'X', highest_freq_aa, freq_val, 'NA') in merged list.
    - Else:
        -> We consume one line from dssp_data[j].
        -> Suppose dssp_data[j] = (pos_dssp, dssp_ref, ss).
        -> We produce (pos_counter, dssp_ref, highest_freq_aa, freq_val, ss).
        -> Then j += 1

    i always goes from 0..(len(freq_data)-1)
    j only advances when ref_aa != 'X'.

    Return merged list of (pos_counter, ref_aa_final, highest_freq_aa, freq_val, ss).
    """
    merged = []
    j = 0  # pointer for dssp_data
    pos_counter = 1

    for i in range(len(freq_data)):
        ref_aa_freq, highest_freq_aa, freq_val = freq_data[i]

        if ref_aa_freq == 'X':
            # Unresolved region: we don't consume dssp_data
            merged.append((pos_counter, 'X', highest_freq_aa, freq_val, 'NA'))
            pos_counter += 1
        else:
            # We should consume one line from dssp_data (if available)
            if j < len(dssp_data):
                _, dssp_ref, ss = dssp_data[j]
                # We decide to use dssp_ref as the final reference
                merged.append((pos_counter, dssp_ref, highest_freq_aa, freq_val, ss))
                j += 1
            else:
                # If there's no more dssp line, we can only store partial info
                merged.append((pos_counter, ref_aa_freq, highest_freq_aa, freq_val, 'NA'))

            pos_counter += 1

    return merged


def main():
    parser = argparse.ArgumentParser(description="Merge out.freq & out.txt (DSSP), handle 'X' lines, then filter by beta/gama.")
    parser.add_argument("-freq", required=True, help="Path to out.freq file (may contain 'X' lines)")
    parser.add_argument("-dssp", required=True, help="Path to out.txt file (DSSP result)")
    parser.add_argument("-beta", required=True, type=float, help="Frequency threshold if SS != 'C'")
    parser.add_argument("-gama", required=True, type=float, help="Frequency threshold if SS == 'C'")
    parser.add_argument("-comb", required=True, help="Output merged file")
    parser.add_argument("-mut", required=True, help="Output filtered file")
    args = parser.parse_args()

    # 1) Parse freq
    freq_data = parse_out_freq(args.freq)

    # 2) Parse dssp
    dssp_data = parse_out_txt(args.dssp)

    # 3) Merge with 'X' handling
    merged_list = merge_data(freq_data, dssp_data)

    # 4) Write combined file
    with open(args.comb, 'w') as fout:
        fout.write("pos\tref_aa\thighest_freq_aa\tfrequency\tSS\n")
        for (pos_counter, ref_aa, hfreq_aa, freq_val, ss) in merged_list:
            fout.write(f"{pos_counter}\t{ref_aa}\t{hfreq_aa}\t{freq_val:.4f}\t{ss}\n")

    # 5) Filter by beta/gama
    beta = args.beta
    gama = args.gama

    filtered_set = set()

    for (pos_counter, ref_aa, hfreq_aa, freq_val, ss) in merged_list:
        # Typically skip 'X' lines from filter. If you want to keep them, remove this check:
        if ref_aa == 'X':
            continue

        # Only keep if ref_aa != highest_freq_aa AND freq_val >= threshold
        if ref_aa != hfreq_aa:
            if ss == 'C':
                if freq_val >= gama:
                    filtered_set.add((pos_counter, ref_aa, hfreq_aa, freq_val, ss))
            else:
                if freq_val >= beta:
                    filtered_set.add((pos_counter, ref_aa, hfreq_aa, freq_val, ss))

    filtered_list = sorted(filtered_set, key=lambda x: x[0])

    # 6) Write filtered file
    with open(args.mut, 'w') as fout:
        fout.write("pos\tref_aa\thighest_freq_aa\tfrequency\tSS\n")
        for (pos_counter, ref_aa, hfreq_aa, freq_val, ss) in filtered_list:
            fout.write(f"{pos_counter}\t{ref_aa}\t{hfreq_aa}\t{freq_val:.4f}\t{ss}\n")

    print(f"[INFO] Combined file saved to {args.comb}")
    print(f"[INFO] Filtered mutations saved to {args.mut}")

if __name__ == "__main__":
    main()
