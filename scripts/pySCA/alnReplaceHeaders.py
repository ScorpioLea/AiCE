#!/usr/bin/env python 
"""
A script that replaces the headers of one fasta file with headers from another. It assumes the sequences in the two fasta files are in identical order. This is useful, for example, when working with Promals3D alignments (which often have the header information truncated). 

:Arguments:
     Headers.fasta      (Alignment that is providing the headers)
     Sequences.fasta    (Alignment that is providing the sequences)  

:Keyword Arguments:
     --output           output file name, default: FixedHeaders.fa

:By: Kim Reynolds
:On: 6.5.2015

Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds
This program is free software distributed under the BSD 3-clause
license, please see the file LICENSE for details.
"""

import scaTools as sca
import argparse 


if __name__ =='__main__':
        #parse inputs
        parser = argparse.ArgumentParser()
        parser.add_argument("alignment_for_headers", help='Alignment that is providing the headers')
        parser.add_argument("alignment_for_seqs", help ='ALignment that is providing the sequences')
        parser.add_argument("--output", dest="outputfile", default = 'FixedHeaders.fa', help="specify an outputfile name")
        options = parser.parse_args()

        print ("WARNING:  This script assumes that the headers of the two input fasta files are in IDENTICAL order. If this is NOT true, the script will give incorrect results");

        headers1, seqs1 = sca.readAlg(options.alignment_for_headers)
        headers2, seqs2 = sca.readAlg(options.alignment_for_seqs)

        if (len(seqs2) != len(headers1)):
            print ("ERROR: The length of the two alignments does not match.")
            exit
                

        f = open(options.outputfile, 'w')
        for i,k in enumerate(headers1):
            f.write('>%s\n' % k)
            f.write('%s\n' % seqs2[i])
        f.close()
    
        
