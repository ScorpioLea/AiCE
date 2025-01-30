
#!/usr/bin/env python
"""
The scaCore script runs the core calculations for SCA, and stores the output using the python tool pickle.
It now also outputs the sca matrix (db['sca']['Csca']) as a separate tsv file in the same folder.

:Example:
>>> ./scaCore.py PF00071_full.db
"""

from __future__ import division
import sys, time
import os
import numpy as np
import copy
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import scaTools as sca
import pickle
import argparse
from Bio import SeqIO
from scipy.stats import t
from scipy.stats import scoreatpercentile
from scipy.io import savemat

if __name__ == '__main__':

    # 1) Parse inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("database", help='database from running scaProcessMSA')
    parser.add_argument("-n", dest = "norm", default='frob',
                        help="norm type for dimension-reducing the sca matrix. Options: 'spec' or 'frob'. Default: frob")
    parser.add_argument("-t", "--Ntrials", dest ="Ntrials", default=10, type=int,
                        help="number of randomization trials (default: 10)")
    parser.add_argument("-l", dest = "lbda", default=0.03, type=float,
                        help="lambda parameter for pseudo-counting the alignment. Default: 0.03")
    parser.add_argument("-m","--matlab", dest = "matfile",  action = "store_true", default = False,
                        help="write out the results of these calculations to a matlab workspace for further analysis")
    options = parser.parse_args()

    # Check norm
    if options.norm not in ['frob','spec']:
        sys.exit("Error: the option -n must be set to 'frob' or 'spec'.")

    # 2) Load database produced by scaProcessMSA
    db_in = pickle.load(open(options.database, "rb"))
    D_in = db_in['sequence']

    msa_num = D_in['msa_num']
    seqw    = D_in['seqw']
    Nseq    = D_in['Nseq']
    Npos    = D_in['Npos']
    ats     = D_in['ats']
    hd      = D_in['hd']

    # 3) Sequence-level analyses
    print("Computing the sequence projections...")
    Useq, Uica = sca.seqProj(msa_num, seqw, kseq=30, kica=15)
    simMat = sca.seqSim(msa_num)

    # 4) SCA calculations
    print("Computing the SCA conservation and correlation values...")
    Wia, Dia, Di = sca.posWeights(msa_num, seqw, options.lbda)
    Csca, tX, Proj = sca.scaMat(msa_num, seqw, options.norm, options.lbda)

    # 5) Matrix randomizations
    print("Computing matrix randomizations...")
    start = time.time()
    Vrand, Lrand, Crand = sca.randomize(msa_num, options.Ntrials, seqw, options.lbda)
    end = time.time()
    print(f"Randomizations complete, {options.Ntrials} trials, time: {(end - start)/60:.1f} minutes")

    # 6) Organize output data
    #    => We'll write them to the same directory as the input .db file
    input_dir = os.path.dirname(options.database)
    if input_dir == "":
        input_dir = "."
    fn        = os.path.basename(options.database)         # e.g. 'PF00071_full.db'
    fn_noext  = os.path.splitext(fn)[0]                    # e.g. 'PF00071_full'

    # 7) Build sca dictionary
    D = {}
    D['Useq']   = Useq
    D['Uica']   = Uica
    D['simMat'] = simMat
    D['lbda']   = options.lbda
    D['Dia']    = Dia
    D['Di']     = Di
    D['Csca']   = Csca
    D['tX']     = tX
    D['Proj']   = Proj
    D['Ntrials'] = options.Ntrials
    D['Vrand'] = Vrand
    D['Lrand'] = Lrand
    D['Crand'] = Crand

    db = {}
    db['sequence'] = D_in
    db['sca']      = D

    # 8) Write out the new .db file in the same folder
    output_db_path = os.path.join(input_dir, fn_noext + ".db")
    print(f"Calculations complete. Writing updated database to: {output_db_path}")
    with open(output_db_path, "wb") as fdb:
        pickle.dump(db, fdb)

    # If user wants .mat file
    if options.matfile:
        output_mat_path = os.path.join(input_dir, fn_noext + ".mat")
        from scipy.io import savemat
        savemat(output_mat_path, db, appendmat=True, oned_as='column')
        print(f"Also wrote MAT file: {output_mat_path}")

    # 9) Also write the sca matrix as a .tsv file (comma-delimited here)
    output_sca_matrix_path = os.path.join(input_dir, fn_noext + ".sca_matrix.tsv")
    np.savetxt(output_sca_matrix_path, Csca, delimiter=',', comments='')
    print(f"SCA matrix file saved as {output_sca_matrix_path}")

    print("Done.")
