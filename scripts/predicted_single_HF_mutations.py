#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import joblib
import numpy as np

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


def round_to_0_1_increment(x):
    """
    Round a float x to the nearest 0.1 increment.
    Example:
      0.03 -> 0.0
      0.05 -> 0.1
      0.28 -> 0.3
      0.34 -> 0.3
      0.35 -> 0.4
    Mathematically: floor(x/0.1 + 0.5) * 0.1
    """
    return int(x / 0.1 + 0.5) * 0.1


def main():
    parser = argparse.ArgumentParser(
        description="Merge out.freq & out.txt (DSSP), then filter by beta/gama. "
                    "Optionally auto-predict beta/gama from pretrained models."
    )
    parser.add_argument("-freq", required=True, help="Path to out.freq file (may contain 'X' lines)")
    parser.add_argument("-dssp", required=True, help="Path to out.txt file (DSSP result)")
    parser.add_argument("-comb", required=True, help="Output merged file")
    parser.add_argument("-mut", required=True, help="Output filtered file")

    # Two ways: either user passes in -beta, -gama or uses -model_path to auto-predict.
    parser.add_argument("-beta", type=float, help="Frequency threshold if SS != 'C'")
    parser.add_argument("-gama", type=float, help="Frequency threshold if SS == 'C'")

    parser.add_argument("-model_path", help="Path with best_model_a.pkl/best_model_b.pkl/scaler.pkl for auto-predict")

    args = parser.parse_args()

    # 1) Parse freq
    freq_data = parse_out_freq(args.freq)

    # 2) Parse dssp
    dssp_data = parse_out_txt(args.dssp)

    # 3) Merge
    merged_list = merge_data(freq_data, dssp_data)

    # 4) Determine beta, gama either from user or auto-model
    beta = args.beta
    gama = args.gama

    if (beta is None or gama is None) and args.model_path:
        # If user didn't provide either beta or gama, but gave a model_path => auto-predict
        model_a_path = f"{args.model_path}/best_model_a.pkl"  # for predicting beta
        model_b_path = f"{args.model_path}/best_model_b.pkl"  # for predicting gama
        scaler_path  = f"{args.model_path}/scaler.pkl"

        # from the dssp file => protein_size, flex_ratio
        with open(args.dssp, 'r') as f:
            lines = [l.strip() for l in f if l.strip()]
        protein_size = len(lines)
        c_count = 0
        for line in lines:
            parts = line.split()
            if len(parts) >= 3 and parts[2] == 'C':
                c_count += 1
        flex_ratio = c_count / protein_size if protein_size > 0 else 0.0

        # load models
        model_a = joblib.load(model_a_path)  # best_model_a => beta
        model_b = joblib.load(model_b_path)  # best_model_b => gama
        scaler  = joblib.load(scaler_path)

        new_data = np.array([[protein_size, flex_ratio]], dtype=float)
        new_data_scaled = scaler.transform(new_data)

        pred_a = model_a.predict(new_data_scaled)
        pred_b = model_b.predict(new_data_scaled)

        # round to nearest 0.1
        auto_beta = round_to_0_1_increment(pred_a[0])
        auto_gama = round_to_0_1_increment(pred_b[0])

        # fill in missing
        if beta is None:
            beta = auto_beta
        if gama is None:
            gama = auto_gama

        print(f"[INFO] Raw predicted => beta={pred_a[0]:.4f}, gama={pred_b[0]:.4f}")
        print(f"[INFO] Rounded to 0.1 => beta={beta:.1f}, gama={gama:.1f}")

    if beta is None or gama is None:
        raise ValueError("Need either -beta/-gama or -model_path to predict them automatically.")

    print(f"[INFO] Final beta={beta}, gama={gama}")

    # 5) Write combined file
    with open(args.comb, 'w') as fout:
        fout.write("pos\tref_aa\thighest_freq_aa\tfrequency\tSS\n")
        for (pos_counter, ref_aa, hfreq_aa, freq_val, ss) in merged_list:
            fout.write(f"{pos_counter}\t{ref_aa}\t{hfreq_aa}\t{freq_val:.4f}\t{ss}\n")

    # 6) Filter by beta/gama
    filtered_set = set()

    for (pos_counter, ref_aa, hfreq_aa, freq_val, ss) in merged_list:
        if ref_aa == 'X':
            continue

        # Must differ from ref_aa, and freq_val >= threshold
        if ref_aa != hfreq_aa:
            if ss == 'C':
                # coil => use gama
                if freq_val >= gama:
                    filtered_set.add((pos_counter, ref_aa, hfreq_aa, freq_val, ss))
            else:
                # non-coil => use beta
                if freq_val >= beta:
                    filtered_set.add((pos_counter, ref_aa, hfreq_aa, freq_val, ss))

    filtered_list = sorted(filtered_set, key=lambda x: x[0])

    # 7) Write filtered file
    with open(args.mut, 'w') as fout:
        fout.write("pos\tref_aa\thighest_freq_aa\tfrequency\tSS\n")
        for (pos_counter, ref_aa, hfreq_aa, freq_val, ss) in filtered_list:
            fout.write(f"{pos_counter}\t{ref_aa}\t{hfreq_aa}\t{freq_val:.4f}\t{ss}\n")

    print(f"[INFO] Combined file saved to {args.comb}")
    print(f"[INFO] Filtered mutations saved to {args.mut}")


if __name__ == "__main__":
    main()
