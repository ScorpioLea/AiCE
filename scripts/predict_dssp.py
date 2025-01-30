#!/usr/bin/env python3

import sys
import subprocess
import os
import shutil

def run_dssp(pdb_file, dssp_file):
    """
    Run mkdssp on the provided PDB file to generate a DSSP file.
    """
    try:
        with open(dssp_file, 'w') as outfile:
            subprocess.run(['mkdssp', pdb_file], check=True, stdout=outfile)
    except subprocess.CalledProcessError as e:
        print(f"Error while running mkdssp: {e}")
        sys.exit(1)

def parse_dssp(dssp_file, output_file):
    """
    Parse the DSSP file and extract secondary structure information.
    Output format: one residue per line: [resnum, one-letter-aa, secondary-structure].
    """
    # Mapping from three-letter amino acid to one-letter
    aa_three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
        'UNK': 'X', 'XAA': 'X', 'XLE': 'J',
    }

    with open(dssp_file, 'r') as dssp, open(output_file, 'w') as out:
        for line in dssp:
            # DSSP lines containing "nohd" typically have secondary structure data
            if line.startswith('nohd'):
                columns = line.strip().split()
                if len(columns) < 5:
                    continue  # Skip incomplete lines

                aa_three = columns[3]
                # If the third column is a digit, it likely isn't a valid residue line
                if aa_three.isdigit():
                    continue

                chain_id = columns[1]
                res_num = columns[2]  # residue number
                ss = columns[4]       # secondary structure symbol

                # Convert three-letter amino acid to one-letter
                aa = aa_three_to_one.get(aa_three.upper(), 'X')

                # Map DSSP secondary structure to a simpler form
                ss_map = {
                    'H': 'H',  # Alpha helix
                    'G': 'H',  # 3-10 helix
                    'I': 'H',  # Pi helix
                    'P': 'H',  # Pi helix (mkdssp might use 'P')
                    'E': 'E',  # Extended strand
                    'B': 'E',  # Beta bridge
                    'T': 'C',  # Turn
                    'S': 'C',  # Bend
                    '.': 'C',  # Coil
                    '?': 'C',  # Unknown
                }
                ss_simple = ss_map.get(ss, 'C')  # Default 'C' (Coil)

                out.write(f"{res_num}\t{aa}\t{ss_simple}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python predict_dssp.py input_structure output.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Check if mkdssp is installed
    if not shutil.which('mkdssp'):
        print("Error: mkdssp is not installed or not in the PATH.")
        sys.exit(1)

    # If input is a .cif file, convert it to .pdb
    file_ext = os.path.splitext(input_file)[1].lower()
    if file_ext == '.cif':
        converted_pdb = input_file + '_temp.pdb'
        try:
            subprocess.run(['cif2pdb', input_file, converted_pdb], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error while running cif2pdb: {e}")
            sys.exit(1)

        # After conversion, check if the first line contains "REMARK cif2pdb"
        # If yes, remove that line
        with open(converted_pdb, 'r') as f_in:
            lines = f_in.readlines()
        if lines and "REMARK cif2pdb" in lines[0]:
            # Overwrite the file without the first line
            with open(converted_pdb, 'w') as f_out:
                f_out.writelines(lines[1:])

        pdb_file = converted_pdb
    else:
        # Otherwise assume it's a .pdb
        pdb_file = input_file

    # Prepare the dssp file name
    dssp_file = pdb_file + '.dssp'

    # Run mkdssp on the pdb file
    run_dssp(pdb_file, dssp_file)

    # Parse the generated DSSP file
    parse_dssp(dssp_file, output_file)

    # Optional: remove the DSSP file
    os.remove(dssp_file)

    # Optional: if we generated a temp pdb from a .cif, you can remove it too
    if file_ext == '.cif':
        os.remove(pdb_file)

    print(f"Secondary structure prediction saved to {output_file}")

if __name__ == '__main__':
    main()
