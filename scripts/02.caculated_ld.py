#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    python run_ld_pipeline.py <seq_dir> <output_ld_dir>

Parameters:
    seq_dir        Directory containing *.fa protein MSA files
    output_ld_dir  Directory to store final LD matrix files (*.ld)

Example:
    python run_ld_pipeline.py ./input_fa ./ld_results

Workflow:
    For each *.fa in <seq_dir>:
      1) Convert protein MSA -> DNA FASTA (inline logic of 'protein2dna.py')
      2) Convert DNA FASTA -> VCF (inline logic of 'fasta2vcf.py', using
         first sequence as reference)
      3) Use plink to compute LD matrix
      4) Convert .ld tabs to commas with 'sed'
      5) Remove intermediate .log, .nosex files

The final <prefix>.ld file will be saved in <output_ld_dir>, where <prefix>
matches the *.fa filename (without extension).
"""

import os
import sys
import glob
import subprocess

from Bio import SeqIO


# Optimal codon table (from your script)
OPTIMAL_CODONS = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
    '-': '---', 'X': '---'  # gap or unknown
}

def translate_protein_to_dna(protein_seq):
    """
    Convert a single protein sequence (Bio.Seq object or string)
    to a DNA sequence using the OPTIMAL_CODONS table.
    """
    dna_seq = []
    for aa in protein_seq:
        aa = aa.upper()
        if aa in OPTIMAL_CODONS:
            dna_seq.append(OPTIMAL_CODONS[aa])
        else:
            print(f"Warning: Amino acid {aa} is not in the table. Using '---'.")
            dna_seq.append('---')
    return "".join(dna_seq)

def convert_protein_fasta_to_dna(input_fasta, output_fasta):
    """
    Read each record from 'input_fasta' (protein),
    convert to DNA, and write to 'output_fasta'.
    """
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            dna_sequence = translate_protein_to_dna(str(record.seq))
            out_f.write(f">{record.id}\n{dna_sequence}\n")


def generate_vcf_header(samples):
    header_lines = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)
    ]
    return "\n".join(header_lines)

def classify_indel_and_adjust(ref_base, alt_bases):
    """
    Simplify insertion/deletion handling:
      - If ref_base == '-', treat it as an insertion in others.
      - If alt_bases contain '-', treat them as <DEL>.
    Returns (new_ref, new_alts, is_indel_flag).
    """
    if ref_base == '-':
        # Reference is a gap => 'N' as a placeholder
        new_ref = 'N'
        new_alts = []
        for base in alt_bases:
            if base == '-':
                # same as reference => no real variation
                pass
            else:
                # there's an actual base => treat as <INS>
                if base not in ['A', 'C', 'G', 'T']:
                    new_alts.append("<INS>")
                else:
                    new_alts.append("<INS>")
        if new_alts:
            return new_ref, list(set(new_alts)), True
        else:
            return ref_base, [], False
    else:
        # Reference is not a gap
        new_alts = []
        is_indel = False
        for base in alt_bases:
            if base == '-':
                new_alts.append("<DEL>")
                is_indel = True
            else:
                new_alts.append(base)
        return ref_base, new_alts, is_indel

def generate_vcf_record(chrom, pos, ref, alts, samples, sequences):
    """
    Build one VCF line for the given position. 
    sequences[0] is reference, sequences[1:] are samples.
    """
    # Build a dict: REF -> '0', ALT -> '1','2','3', ...
    alt_dict = {ref: '0'}
    for idx, alt_allele in enumerate(alts):
        alt_dict[alt_allele] = str(idx + 1)

    # Construct genotypes for each sample
    genotypes = []
    for sequence in sequences[1:]:
        b = sequence.seq[pos]
        gt = alt_dict.get(b, '.')
        # Assume diploid => x|x
        genotype = f"{gt}|{gt}"
        genotypes.append(genotype)

    # Basic INFO field
    info_field = f"NS={len(samples)}"
    alt_field = ','.join(alts) if alts else '.'
    format_field = "GT"
    record = f"{chrom}\t{pos+1}\t.\t{ref}\t{alt_field}\t.\t.\t{info_field}\t{format_field}\t" + "\t".join(genotypes)
    return record

def convert_fasta_to_vcf(dna_fasta, output_vcf):
    sequences = list(SeqIO.parse(dna_fasta, "fasta"))
    if not sequences:
        raise ValueError("No sequences found in the FASTA file.")

    ref_sequence = sequences[0].seq
    samples = [seq.id for seq in sequences[1:]]

    with open(output_vcf, 'w') as vcf:
        # Write VCF header
        vcf.write(generate_vcf_header(samples) + "\n")
        # Iterate each position in reference
        for pos in range(len(ref_sequence)):
            ref_base = ref_sequence[pos]
            # Gather alt bases from samples
            base_counts = {}
            for seq_obj in sequences[1:]:
                b = seq_obj.seq[pos]
                base_counts[b] = base_counts.get(b, 0) + 1

            # Filter out any base same as ref (unless ref == '-')
            alt_bases = [b for b, c in base_counts.items() if (b != ref_base and c > 0)]

            new_ref, new_alts, is_indel = classify_indel_and_adjust(ref_base, alt_bases)
            if not new_alts:
                continue  # no variation => skip

            record = generate_vcf_record("chr1", pos, new_ref, new_alts, samples, sequences)
            vcf.write(record + "\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python run_ld_pipeline.py <seq_dir> <output_ld_dir>")
        sys.exit(1)

    seq_dir = sys.argv[1]
    output_ld_dir = sys.argv[2]

    os.makedirs(output_ld_dir, exist_ok=True)

    fa_files = glob.glob(os.path.join(seq_dir, "*.fa"))
    if not fa_files:
        print(f"[ERROR] No *.fa files found in {seq_dir}.")
        sys.exit(2)

    for fa_path in fa_files:
        filename = os.path.basename(fa_path)  # e.g. "mymsa.fa"
        prefix = os.path.splitext(filename)[0]  # e.g. "mymsa"

        print(f"\n[INFO] Processing file: {fa_path} with prefix: {prefix}")

        # Step A: Convert protein MSA -> DNA
        dna_fasta = os.path.join(output_ld_dir, f"{prefix}_dna.fasta")
        convert_protein_fasta_to_dna(fa_path, dna_fasta)
        print(f"  [Done] Protein -> DNA: {dna_fasta}")

        # Step B: Convert DNA FASTA -> VCF
        vcf_file = os.path.join(output_ld_dir, f"{prefix}.vcf")
        convert_fasta_to_vcf(dna_fasta, vcf_file)
        print(f"  [Done] DNA FASTA -> VCF: {vcf_file}")

        # Step C: Run plink to compute LD (r2 matrix)
        plink_out_prefix = os.path.join(output_ld_dir, prefix)
        plink_cmd = [
            "plink", 
            "--vcf", vcf_file,
            "--r2", "square", 
            "--out", plink_out_prefix,
            "--allow-no-sex",
            "--keep-allele-order",
            "--geno", "1",
            "--mind", "1",
            "--maf", "0.00001",
            "--hwe", "0",
            "--snps-only",
            "--allow-extra-chr"
        ]
        subprocess.run(plink_cmd, check=True)
        ld_file = f"{plink_out_prefix}.ld"

        # Step D: Convert tabs to commas in .ld (in-place)
        if os.path.isfile(ld_file):
            sed_cmd = ["sed", "s/\t/,/g", ld_file, "-i"]
            subprocess.run(sed_cmd, check=True)
        else:
            print(f"  [WARNING] LD file not found: {ld_file}. Check plink output.")

        # Step E: Remove .nosex, .log if they exist
        nosex_file = f"{plink_out_prefix}.nosex"
        log_file   = f"{plink_out_prefix}.log"
        dna_fasta = f"{plink_out_prefix}_dna.fasta"
        for tmp in [nosex_file, log_file, dna_fasta]:
            if os.path.isfile(tmp):
                os.remove(tmp)

        print(f"[INFO] Done with {prefix}. LD matrix file => {ld_file}")

    print("\nAll jobs completed. LD files are located in:", output_ld_dir)

if __name__ == "__main__":
    main()
