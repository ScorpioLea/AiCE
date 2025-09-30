#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
msa_to_pssm.py
Generate PSSM (Position-Specific Scoring Matrix) from MSA in FASTA format.
- The first sequence is treated as wild-type (WT) and is not included in frequency counting.
- Only 20 standard amino acids (ACDEFGHIKLMNPQRSTVWY) are counted; gaps and non-standard symbols are ignored.
- Pseudocounts are used for smoothing; calculates log2(P(a|j) / q(a)).
- Background frequency options: uniform or robinson.
- Output CSV with columns: pos, wt, A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y

Usage:
    python msa_to_pssm.py --msa input.fasta --out pssm.csv
Optional:
    --pseudocount 1.0
    --background uniform|robinson
    --ignore-wt-gaps  (if specified, skip positions where WT is '-')
    --verbose        (display detailed statistics)
    --min-valid-seqs (minimum valid sequences per column, default 1)
"""

import argparse
import sys
import math
from collections import Counter
from pathlib import Path

AA20 = "ACDEFGHIKLMNPQRSTVWY"

# Robinson & Robinson (1991) background frequencies (in AA20 order)
ROBINSON_BG = {
    'A': 0.07805, 'C': 0.01925, 'D': 0.05364, 'E': 0.06321, 'F': 0.03856,
    'G': 0.07415, 'H': 0.02616, 'I': 0.05368, 'K': 0.05816, 'L': 0.09019,
    'M': 0.02243, 'N': 0.04487, 'P': 0.05203, 'Q': 0.04264, 'R': 0.05129,
    'S': 0.0712,  'T': 0.05841, 'V': 0.06441, 'W': 0.0133,  'Y': 0.03216
}


def read_fasta(path):
    """
    Read multiple sequence alignment in FASTA format.
    
    Args:
        path (str): Path to FASTA file
        
    Returns:
        tuple: (list of sequence names, list of sequences)
        
    Raises:
        FileNotFoundError: File does not exist
        ValueError: File is empty, format is incorrect, or sequences have inconsistent lengths
    """
    if not Path(path).exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    names, seqs = [], []
    name, buf = None, []
    
    with open(path, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seq = ''.join(buf)
                    if not seq:
                        raise ValueError(f"Sequence '{name}' is empty (line {line_num})")
                    names.append(name)
                    seqs.append(seq)
                name = line[1:].strip()
                if not name:
                    name = f"seq_{len(names) + 1}"
                buf = []
            else:
                buf.append(line.strip())
        
        # Process the last sequence
        if name is not None:
            seq = ''.join(buf)
            if not seq:
                raise ValueError(f"Sequence '{name}' is empty")
            names.append(name)
            seqs.append(seq)
    
    if not seqs:
        raise ValueError("FASTA is empty or format is incorrect.")
    
    # Check length consistency
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        length_info = '\n'.join([f"  {name}: {len(seq)}" for name, seq in zip(names, seqs)])
        raise ValueError(f"MSA sequences have inconsistent lengths:\n{length_info}")
    
    return names, seqs


def get_background(mode):
    """
    Get amino acid background frequencies.
    
    Args:
        mode (str): 'uniform' or 'robinson'
        
    Returns:
        dict: Mapping from amino acid to frequency
    """
    if mode == "uniform":
        return {aa: 1.0/20.0 for aa in AA20}
    elif mode == "robinson":
        return dict(ROBINSON_BG)
    else:
        raise ValueError(f"background only supports uniform or robinson, got: {mode}")


def safe_log2(x):
    """Safe log2 calculation, avoiding x <= 0."""
    if x <= 0:
        return float('-inf')
    return math.log2(x)


def analyze_msa(msa_seqs, verbose=False):
    """
    Analyze basic statistics of the MSA.
    
    Args:
        msa_seqs (list): List of sequences
        verbose (bool): Whether to output detailed information
        
    Returns:
        dict: Statistical information
    """
    stats = {
        'num_seqs': len(msa_seqs),
        'length': len(msa_seqs[0]) if msa_seqs else 0,
        'gap_positions': 0,
        'total_gaps': 0,
        'total_non_standard': 0
    }
    
    if len(msa_seqs) < 2:
        return stats
    
    wt = msa_seqs[0]
    others = msa_seqs[1:]
    
    for j in range(stats['length']):
        col_gaps = 0
        col_non_standard = 0
        
        for s in others:
            if s[j] == '-':
                col_gaps += 1
            elif s[j] not in AA20:
                col_non_standard += 1
        
        stats['total_gaps'] += col_gaps
        stats['total_non_standard'] += col_non_standard
        
        if wt[j] == '-':
            stats['gap_positions'] += 1
    
    if verbose:
        print(f"[Stats] Total sequences: {stats['num_seqs']} (WT + {stats['num_seqs']-1} others)", file=sys.stderr)
        print(f"[Stats] Sequence length: {stats['length']}", file=sys.stderr)
        print(f"[Stats] Gap positions in WT: {stats['gap_positions']}", file=sys.stderr)
        print(f"[Stats] Total gaps in non-WT sequences: {stats['total_gaps']}", file=sys.stderr)
        print(f"[Stats] Total non-standard residues: {stats['total_non_standard']}", file=sys.stderr)
    
    return stats


def msa_to_pssm(msa_seqs, pseudocount=1.0, background="uniform", 
                ignore_wt_gaps=False, min_valid_seqs=1, verbose=False):
    """
    Generate PSSM from MSA.
    
    Args:
        msa_seqs (list[str]): List of sequences, first one is WT
        pseudocount (float): Pseudocount value
        background (str): Background frequency mode
        ignore_wt_gaps (bool): Whether to ignore positions where WT is a gap
        min_valid_seqs (int): Minimum valid sequences per column
        verbose (bool): Whether to display detailed information
        
    Returns:
        list[list]: Each row is [pos(1-based), wt, 20 PSSM scores for each AA]
        
    Raises:
        ValueError: Invalid parameters or insufficient data
    """
    if pseudocount <= 0:
        raise ValueError(f"Pseudocount must be greater than 0, got: {pseudocount}")
    
    if min_valid_seqs < 1:
        raise ValueError(f"min_valid_seqs must be >= 1, got: {min_valid_seqs}")
    
    bg = get_background(background)
    L = len(msa_seqs[0])
    
    if len(msa_seqs) < 2:
        raise ValueError("MSA has only one WT sequence, cannot compute statistics.")

    wt = msa_seqs[0]
    others = msa_seqs[1:]
    
    rows = []
    skipped_positions = []
    low_coverage_positions = []
    
    for j in range(L):
        wt_res = wt[j]
        
        # Check if we should ignore WT gap positions
        if ignore_wt_gaps and wt_res == '-':
            skipped_positions.append(j + 1)
            continue

        # Collect valid amino acids in this column (excluding WT)
        col = []
        for s in others:
            aa = s[j]
            if aa in AA20:  # Only count 20 standard amino acids
                col.append(aa)
            # Ignore gaps '-' and non-standard residues (e.g., B, Z, X)

        # Check valid sequence count
        valid_count = len(col)
        if valid_count < min_valid_seqs:
            if verbose:
                low_coverage_positions.append((j + 1, valid_count))
            # Still calculate even with insufficient valid sequences (using pseudocounts)
            # If there are no valid sequences at all, PSSM will be entirely determined by pseudocounts

        counts = Counter(col)

        # Pseudocount + normalization to get P(a|j)
        probs = {}
        total = sum(counts.get(aa, 0) for aa in AA20) + pseudocount * 20.0
        for aa in AA20:
            probs[aa] = (counts.get(aa, 0) + pseudocount) / total

        # Calculate PSSM = log2( P(a|j) / q(a) )
        pssm_scores = []
        for aa in AA20:
            p = probs[aa]
            q = bg[aa]
            pssm_scores.append(safe_log2(p / q))

        rows.append([j+1, wt_res] + pssm_scores)

    # Display warning messages
    if verbose:
        if skipped_positions:
            print(f"[Warning] Skipped {len(skipped_positions)} positions where WT is a gap", file=sys.stderr)
        if low_coverage_positions:
            print(f"[Warning] {len(low_coverage_positions)} positions have valid sequences < {min_valid_seqs}:", 
                  file=sys.stderr)
            for pos, count in low_coverage_positions[:10]:  # Only show first 10
                print(f"  Position {pos}: {count} valid sequences", file=sys.stderr)
            if len(low_coverage_positions) > 10:
                print(f"  ... and {len(low_coverage_positions) - 10} more positions", file=sys.stderr)
    
    if not rows:
        raise ValueError("No PSSM rows generated. Please check input parameters.")
    
    return rows


def write_pssm_csv(rows, output_path, precision=6):
    """
    Write PSSM to CSV file.
    
    Args:
        rows (list): PSSM data rows
        output_path (str): Output file path
        precision (int): Floating point precision
    """
    import csv
    
    header = ["pos", "wt"] + list(AA20)
    
    # Ensure output directory exists
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w", newline="", encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            pos, wt, *scores = r
            formatted_scores = [f"{x:.{precision}f}" if x != float('-inf') else '-inf' 
                              for x in scores]
            w.writerow([pos, wt] + formatted_scores)


def main():
    parser = argparse.ArgumentParser(
        description="Generate PSSM from FASTA MSA (excluding first WT sequence from counts)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic usage:
    python msa_to_pssm.py --msa input.fasta --out pssm.csv
  
  Using Robinson background frequencies and custom pseudocount:
    python msa_to_pssm.py --msa input.fasta --out pssm.csv --background robinson --pseudocount 0.5
  
  Ignore WT gap positions and show detailed information:
    python msa_to_pssm.py --msa input.fasta --out pssm.csv --ignore-wt-gaps --verbose
        """
    )
    
    parser.add_argument("--msa", required=True, help="Input MSA file (FASTA format)")
    parser.add_argument("--out", required=True, help="Output CSV file path")
    parser.add_argument("--pseudocount", type=float, default=1.0, 
                       help="Pseudocount value for probability smoothing (default: 1.0)")
    parser.add_argument("--background", choices=["uniform", "robinson"], default="uniform",
                       help="Background frequency mode: uniform (uniform distribution) or robinson (Robinson-Robinson 1991) (default: uniform)")
    parser.add_argument("--ignore-wt-gaps", action="store_true",
                       help="If specified, skip columns where WT position is '-'")
    parser.add_argument("--min-valid-seqs", type=int, default=1,
                       help="Minimum valid sequences per column; warnings issued if below this (default: 1)")
    parser.add_argument("--precision", type=int, default=6,
                       help="Decimal places for output values (default: 6)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Display detailed statistics and warnings")
    
    args = parser.parse_args()
    
    # Parameter validation
    if args.pseudocount <= 0:
        parser.error("--pseudocount must be greater than 0")
    
    if args.min_valid_seqs < 1:
        parser.error("--min-valid-seqs must be >= 1")
    
    if args.precision < 0:
        parser.error("--precision must be >= 0")
    
    try:
        # Read MSA
        if args.verbose:
            print(f"[Info] Reading MSA file: {args.msa}", file=sys.stderr)
        
        names, seqs = read_fasta(args.msa)
        
        if args.verbose:
            print(f"[Info] Successfully read {len(seqs)} sequences of length {len(seqs[0])}", file=sys.stderr)
            print(f"[Info] WT sequence: {names[0]}", file=sys.stderr)
        
        # Analyze MSA
        analyze_msa(seqs, verbose=args.verbose)
        
        # Generate PSSM
        if args.verbose:
            print(f"[Info] Generating PSSM...", file=sys.stderr)
            print(f"[Param] Pseudocount: {args.pseudocount}", file=sys.stderr)
            print(f"[Param] Background frequency: {args.background}", file=sys.stderr)
            print(f"[Param] Ignore WT gaps: {args.ignore_wt_gaps}", file=sys.stderr)
        
        rows = msa_to_pssm(
            seqs,
            pseudocount=args.pseudocount,
            background=args.background,
            ignore_wt_gaps=args.ignore_wt_gaps,
            min_valid_seqs=args.min_valid_seqs,
            verbose=args.verbose
        )
        
        # Write CSV
        if args.verbose:
            print(f"[Info] Writing PSSM to: {args.out}", file=sys.stderr)
        
        write_pssm_csv(rows, args.out, precision=args.precision)
        
        print(f"✓ Successfully generated PSSM with {len(rows)} positions, saved to: {args.out}", file=sys.stderr)
        
    except FileNotFoundError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"✗ Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
