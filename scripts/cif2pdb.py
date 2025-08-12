#!/usr/bin/env python3
"""
Robust mmCIF -> PDB converter (Python)

Reimplements the C++ `cif2pdb` behavior using BioPython (no manual tokenizing),
including:
  - `-chain`   : filter by author chain ID
  - `-mol`     : molecule type mask (1=protein, 2=RNA, 4=DNA; default 7=all)
  - `-split`   : 0 single file (default), 1 split one PDB per chain
  - `-het`     : 0 ATOM only; 1 include MSE and convert to MET; 2 include all HETATM except HOH; 3 include all HETATM incl. HOH
  - MSE->MET   : when `-het>=1`, convert residue name and SE->SD
  - Chain IDs  : resolve to unique single-letter IDs when not splitting
  - Compressed : reads .cif, .cif.gz, or .cif.bz2

Usage
-----
python cif2pdb.py input.cif output.pdb [-chain A] [-mol 7] [-split 0|1] [-het 1]

Requires: pip install biopython
"""

import argparse
import sys
import os
import gzip
import bz2
from io import StringIO
from typing import Dict, Set

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO, Select

CHAIN_CHARSET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

# ------------------------- helpers -------------------------

def open_text_auto(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    if path.endswith('.bz2'):
        return bz2.open(path, 'rt')
    return open(path, 'rt')


def classify_polymer(resname: str) -> str:
    """Rudimentary classification to emulate the C++ logic.
    Returns one of: 'protein', 'rna', 'dna', or 'other'."""
    rn = (resname or '').strip()
    if len(rn) == 1:
        return 'rna'  # A/U/C/G
    if len(rn) == 2 and rn.startswith('D'):
        return 'dna'  # DA/DG/DC/DT
    if len(rn) >= 4:
        rn = rn[:3]
    if len(rn) == 3:
        return 'protein'
    return 'other'


def map_chain_ids_one_letter(model, keep_chain: str | None) -> Dict[str, str]:
    """Create unique single-letter chain IDs (in-place), return old->new map.
    Mimics the C++ behavior by preferring the last char of the auth asym_id when possible."""
    taken: Set[str] = set()
    mapping: Dict[str, str] = {}

    # First, attempt to keep unique single chars
    for ch in model:
        if keep_chain and ch.id != keep_chain:
            continue
        old = ch.id
        if len(old) == 1 and old not in taken:
            mapping[old] = old
            taken.add(old)

    # Assign remaining
    for ch in list(model):
        if keep_chain and ch.id != keep_chain:
            continue
        old = ch.id
        if old in mapping:
            continue
        cand = old[-1] if old else ' '
        if len(cand) != 1 or cand in taken:
            cand = None
            for c in CHAIN_CHARSET:
                if c not in taken:
                    cand = c
                    break
        if cand is None:
            # drop this chain if truly impossible
            # (extremely unlikely in practice)
            model.detach_child(old)
            continue
        mapping[old] = cand
        taken.add(cand)

    # Apply
    for ch in list(model):
        if keep_chain and ch.id != keep_chain:
            model.detach_child(ch.id)
            continue
        ch.id = mapping[ch.id]

    return mapping


# ------------------------- selection -------------------------
class CIF2PDBSelect(Select):
    def __init__(self, chain_opt: str, mol_mask: int, het_opt: int):
        super().__init__()
        self.chain_opt = chain_opt or ''
        self.mol_mask = mol_mask
        self.het_opt = het_opt
        self.per_res_altloc = {}  # (chain, id) -> chosen altloc

    def accept_chain(self, chain):
        if self.chain_opt:
            return chain.id == self.chain_opt
        return True

    def accept_residue(self, residue):
        # HET flag in Biopython: residue.id[0] == ' ' for ATOM, 'H_' for HETATM, 'W' for water
        hetflag = residue.id[0]
        name = residue.get_resname().strip()

        # Water handling
        if hetflag == 'W' or name in {'HOH', 'WAT', 'H2O'}:
            return self.het_opt == 3  # only include when het=3

        # HETATM handling
        if hetflag != ' ':
            if self.het_opt == 0:
                return False
            if self.het_opt == 1:
                return name == 'MSE'
            if self.het_opt == 2:
                return name not in {'HOH', 'WAT', 'H2O'}
            return True  # het_opt == 3

        # Polymer residue: filter by mol type
        ptype = classify_polymer(name)
        if ptype == 'protein' and not (self.mol_mask & 1):
            return False
        if ptype == 'rna' and not (self.mol_mask & 2):
            return False
        if ptype == 'dna' and not (self.mol_mask & 4):
            return False
        return True

    def accept_atom(self, atom):
        # AltLoc: choose ' ' or the first observed per residue; prefer 'A'
        alt = atom.get_altloc()
        if not alt or alt == ' ':
            return True
        # Per-residue choice
        res = atom.get_parent()
        key = (res.get_parent().id, res.id)
        chosen = self.per_res_altloc.get(key)
        if chosen is None:
            self.per_res_altloc[key] = 'A' if alt == 'A' else alt
            chosen = self.per_res_altloc[key]
        return alt == chosen


# ------------------------- main conversion -------------------------

def convert(infile: str, outfile: str, chain_opt: str, mol_mask: int, split: int, het: int):
    # Parse mmCIF via Biopython
    parser = MMCIFParser(QUIET=True)
    with open_text_auto(infile) as handle:
        structure = parser.get_structure(os.path.basename(infile)[:4] or 'stru', handle)

    # Work on the first model only (C++ behavior)
    model = next(structure.get_models())

    # MSE -> MET if het>=1: mutate in-place before writing
    if het >= 1:
        for chain in model:
            for res in chain:
                if res.get_resname().strip() == 'MSE':
                    res.resname = 'MET'
                    for atom in res:
                        if atom.get_name().strip() == 'SE':
                            atom.set_name('SD')

    # Chain filtering/mapping
    sel_all = CIF2PDBSelect(chain_opt='', mol_mask=mol_mask, het_opt=het)

    if split == 1:
        # Write each chain separately; do not remap IDs
        io = PDBIO()
        io.set_structure(structure)
        for chain in list(model):
            if chain_opt and chain.id != chain_opt:
                continue
            sel = CIF2PDBSelect(chain_opt=chain.id, mol_mask=mol_mask, het_opt=het)
            outpath = '-'
            if outfile != '-':
                prefix, ext = os.path.splitext(outfile)
                outpath = f"{prefix}{chain.id}.pdb"
            if outpath == '-':
                buf = StringIO()
                io.save(buf, select=sel)
                sys.stdout.write(f"REMARK cif2pdb {chain.id} {chain.id}\n")
                sys.stdout.write(buf.getvalue())
                sys.stdout.write("TER\n")
            else:
                with open(outpath, 'w') as fout:
                    buf = StringIO()
                    io.save(buf, select=sel)
                    fout.write(f"REMARK cif2pdb {chain.id} {chain.id}\n")
                    fout.write(buf.getvalue())
                    if not buf.getvalue().endswith('\n'):
                        fout.write('\n')
                    fout.write("TER\n")
    else:
        # Single file: ensure unique one-letter chain IDs
        # First remap chains (in-place)
        mapping = map_chain_ids_one_letter(model, keep_chain=chain_opt or None)

        io = PDBIO()
        io.set_structure(structure)
        sel = CIF2PDBSelect(chain_opt=chain_opt and mapping.get(chain_opt, chain_opt) or '',
                            mol_mask=mol_mask, het_opt=het)
        if outfile == '-':
            buf = StringIO()
            io.save(buf, select=sel)
            # Write REMARK lines for each non-empty chain
            sys.stdout.write(buf.getvalue())
            if not buf.getvalue().endswith('\n'):
                sys.stdout.write('\n')
            # Bio.PDB already writes END; no need to add
        else:
            with open(outfile, 'w') as fout:
                io.save(fout, select=sel)


# ------------------------- CLI -------------------------

def main():
    ap = argparse.ArgumentParser(description='Convert mmCIF to PDB (BioPython)')
    ap.add_argument('infile', help='mmCIF input file (.cif/.cif.gz/.cif.bz2)')
    ap.add_argument('outfile', nargs='?', default='-', help='PDB output (default stdout)')
    ap.add_argument('-chain', default='', help='Auth chain ID to convert')
    ap.add_argument('-mol', type=int, default=7, help='1=protein,2=RNA,4=DNA (bitmask); default 7')
    ap.add_argument('-split', type=int, choices=[0, 1], default=0, help='Split PDB by chain')
    ap.add_argument('-het', type=int, choices=[0, 1, 2, 3], default=1,
                    help='0 ATOM only; 1 include MSE; 2 include all HETATM except HOH; 3 include all HETATM incl. HOH')
    args = ap.parse_args()

    try:
        convert(args.infile, args.outfile, args.chain, args.mol, args.split, args.het)
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        raise


if __name__ == '__main__':
    main()
