#!/usr/bin/env python
"""
A basic script to filter a fasta file of sequences by size - a useful step to remove partial sequences or sequences that would potentially introduce a large number of gaps in the alignment. This script reads in the alignment, computes the average sequence length, and outputs a new alignment that keeps sequences of length mean +/- tolerance (tolerance default = 50)

:Arguments:
    Input_MSA.fasta (the alignment to be processed)

:Keyword Arguments:
    --tolerance, -t      allowable sequence length variation (in number of amino acids), default: 50
    --output             output file name, default: FilteredAln.fa

:By: Kim Reynolds
:On: 6.5.2015
 
Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds
This program is free software distributed under the BSD 3-clause
license, please see the file LICENSE for details.
"""

import scaTools as sca
import numpy as np
import argparse

if __name__ =='__main__':
        #parse inputs
        parser = argparse.ArgumentParser()
        parser.add_argument("alignment", help='Input Sequence Alignment')
        parser.add_argument("-t","--tolerance", dest = "tol", type = int, default = 50, help="allowable sequence length variation in number of amino acids (alignment will be trimmed to mean +/- tolerance, default = 50)")
        parser.add_argument("--output", dest="outputfile", default = 'FilteredAln.fa', help="specify an outputfile name")

        options = parser.parse_args()

        headers, seqs = sca.readAlg(options.alignment)
        seqLen = np.zeros((len(seqs),1)).astype(int)
        for i,k in enumerate(seqs):
            seqLen[i] = len(k)
        avgLen = seqLen.mean()
        print ("Average sequence length: %i" % avgLen)
        print ("Min: %i, Max %i" % (seqLen.min(), seqLen.max()))
        minsz = avgLen - options.tol;
        maxsz = avgLen + options.tol;
        print ("Keeping sequences in the range: %i - %i" % (minsz, maxsz))

        keepSeqs = list()
        keepHeaders = list()
        for i,k in enumerate(seqLen):
            if (k > minsz) & (k < maxsz):
                keepSeqs.append(seqs[i])
                keepHeaders.append(headers[i])
        
        print ("Keeping %i of %i total sequences" % (len(keepSeqs), len(seqLen)))
        
        f = open(options.outputfile, 'w')
        for i,k in enumerate(keepSeqs):
            f.write('>%s\n' % keepHeaders[i])
            f.write('%s\n' % keepSeqs[i])
        f.close()
    
            
