#!/usr/bin/env python
"""
A script to parse GI numbers from the headers of an alignment with typical Blast formatting.

:Arguments:
    Input_MSA.fasta (the alignment to be processed)

:Keyword Arguments:
    --output             output file name, default: FilteredAln.fa

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
    parser.add_argument("alignment", help='Input Sequence Alignment')
    parser.add_argument("--output", dest="outputfile", default = 'GI_Num', help="specify an outputfile name")
    options = parser.parse_args()

    headers, seqs = sca.readAlg(options.alignment)
    gis = [h.split('_')[1] for h in headers]
    
    f = open(options.outputfile, 'w')
    for gi in gis:
        f.write('%s\n' % gi)
    f.close()
