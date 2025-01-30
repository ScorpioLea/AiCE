"""
pySCA - A SCA toolbox in Python

:By:  

|      Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
|      Kimberly Reynolds (kimberly.reynolds@utsouthwestern.edu)
|      Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)

:On:  August 2014
:version: 6.1 beta

Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds
This program is free software distributed under the BSD 3-clause
license, please see the file LICENSE for details.
"""
from __future__ import division
import os
import subprocess
import copy
import time
import numpy as np
import random as rand
import scipy.sparse
import scipy.sparse.linalg
from scipy.sparse import csr_matrix as sparsify
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import colorsys
import shutil
from Bio.PDB.PDBParser import PDBParser
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import t
from scipy.stats import scoreatpercentile 
from optparse import OptionParser

##########################################################################################
# PATHS
# These have to be changed to be consistent with user-defined paths.

# (this directory should contain the file 'pfamseq.txt' from
# ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/database_files/
path2pfamseq = '~/Documents/Packages/pfamseq.txt'

# the location of your PDB structures
path2structures = 'Inputs/'

# paths to pymol and needle (EMBOSS) applictaions
path2pymol = '/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL'
path2needle = '/usr/local/bin/'

# Also assumes that a folder named 'Outputs' is in the path

##########################################################################################
# CLASSES

class Unit:
    ''' A class for units (sectors, sequence families, etc.) 

        **Attributes:** 
            -  `name`  = string describing the unit (ex: 'firmicutes')
            -  `items` = set of member items (ex: indices for all firmicutes sequences in an alignment)
            -  `col`   = color code associated to the unit (for plotting)
            -  `vect`  = an additional vector describing the member items (ex: a list of sequence weights)

    '''
    def __init__(self):
        self.name = ''
        self.items = set()
        self.col = 0
        self.vect = 0

class Annot:
    ''' A class for annotating sequences 
    
            :Attributes:
                -  `descr` = description (often the sequence header)
                -  `species` = species string
                -  `taxo` = taxonomy string
                -  `seq` = sequence
        '''
    def __init__(self, descr, species, taxo, seq=''):
        self.descr = descr
        self.species = species
        self.taxo = taxo
        self.seq = seq

##########################################################################################
# ALIGNMENT PROCESSING

def readAlg(filename):
    ''' Read in a multiple sequence alignment in fasta format, and return the 
    headers and sequences.

    >>> headers, sequences = readAlg(filename) '''
    filelines = open(filename, 'r').readlines()
    headers = list(); sequences = list(); notfirst = 0
    for line in filelines:
        if line[0] == '>':
            if notfirst > 0: sequences.append(seq.replace('\n','').upper())
            headers.append(line[1:].replace('\n',''))
            seq = ''; notfirst = 1
        elif line != '\n': seq += line
    sequences.append(seq.replace('\n','').upper())
    return headers, sequences

def AnnotPfam(pfam_in, pfam_out, pfam_seq=path2pfamseq):
    ''' Phylogenetic annotation of a Pfam alignment (in fasta format) using information from pfamseq.txt. The output is a fasta file containing phylogenetic annotations in the header (to be parsed with '|' as a delimiter).
    
    Note: the headers for the original alignment take the form >AAA/x-y.  If two entries have same AAA but correspond to different sequences only one of the two sequences will be represented (twice) in the output - this should however not practically be an issue. 
    
    :Arguments:
        -  input PFAM sequence alignment
        -  output file name for the annotated PFAM alignment

    :Keyword Arguments:
        -  `pfam_seq` = path to the file pfamseq.txt

    '''
    start_time=time.time()
    print('Beginning annotation')
    # Reads the pfam headers and sequences:
    headers, sequences = readAlg(pfam_in)
    pfamseq_ids = [h.split('/')[0] for h in headers]
    # Reads the sequence information for those sequences:
    seq_info = dict()
    with open(pfam_seq) as fp:
        for line in fp:
            pf_id = line.split('\t')[1]
            if pf_id in pfamseq_ids:
                seq_info[pf_id] = line
                pfamseq_ids.remove(pf_id)
    end_time=time.time()
    # Writes in output file:
    f = open(pfam_out,'w')
    pfamseq_ids = [h.split('/')[0] for h in headers]
    for i, key in enumerate(pfamseq_ids):
        print('Current step %i, key %s' % (i,key))
        try:
            info = seq_info[key]
        except:
            info = '\t'.join(['unknown']*10 + ['unknown;unknown'])
        #this f.write line works with older pfamseq.txt files (release 4 and before, was 
        #used to annotate the tutorial alignments
        #f.write('>%s|%s|%s|%s\n' % (key, info.split('\t')[6], info.split('\t')[9],\
        #        ','.join([name.strip() for name in info.split('\t')[10].split(';')])))
        #this f.write line works with the new version of pfamseq.txt
        f.write('>%s|%s|%s|%s\n' % (key, info.split('\t')[5], info.split('\t')[8],\
                ','.join([name.strip() for name in info.split('\t')[9].split(';')])))
        f.write('%s\n' % (sequences[i]))
    f.close()
    print('Elapsed time: %.1f min' % ((end_time-start_time)/60))

def clean_al(alg, code='ACDEFGHIKLMNPQRSTVWY', gap='-'):
    ''' Replaces any character that is not a valid amino acid by a gap. 

        :Arguments: amino acid sequence alignment

        :Keyword Arguments:
            -  `code` = list of valid amino acid characters (case sensitive) 
            -  `gap` =  gap character for replacement
    '''
    alg_clean = list()
    for seq in alg:
        seq_clean = ''
        for iaa in range(len(seq)):
            if  code.find(seq[iaa]) != -1: seq_clean += seq[iaa]
            else: seq_clean += gap
        alg_clean.append(seq_clean)
    return alg_clean
    
def MSAsearch(hd, algn, seq, species = None, path2_algprog=path2needle):
    ''' 
    Identify the sequence in the alignment that most closely corresponds to the species of the reference sequence, and return its index.

         **Arguments:**
           -  sequence alignment headers
           -  alignment sequences
           -  selected reference sequence (often from a PDB file)

         **Keyword Arguments:**
           -  `species` =  species of the reference sequence (Used to speed up alignment searching when possible)
           -  `path2_algprog` =  path to an alignment program

         **Example:**
           >>> strseqnum = MSASearch(hd, alg0, pdbseq, 'Homo sapiens')
    '''
    if species is not None:
        species = species.lower()
        key_list = list()
        for (i,h) in enumerate(hd):
            if species in h.lower():
                key_list.append(i)
        hd = [hd[k] for k in key_list]
        algn = [algn[k] for k in key_list]

    try:
        print("Trying MSASearch with ggsearch")
        if not os.path.exists('tmp/'): os.makedirs('tmp/')  
        output_handle = open('tmp/PDB_seq.fasta', 'w')
        SeqIO.write(SeqRecord(Seq(seq), id='PDB sequence'), output_handle, "fasta") 
        output_handle.close()
        f = open("tmp/algn_seq.fasta", "w")
        for i in range(len(algn)):
            f.write(">" + hd[i] + "\n")
            f.write (algn[i] + "\n")
        f.close()
        args = ['ggsearch36','-M 1-'+str(len(algn[0])),'-b','1','-m 8','tmp/PDB_seq.fasta','tmp/algn_seq.fasta']
        output = subprocess.check_output(args)
        i_0 = [i for i in range(len(hd)) if output.split('\t')[1] in hd[i]]  
        if species is not None:
            strseqnum = key_list[i_0[0]]
        else:
            strseqnum = i_0[0]
        shutil.rmtree('tmp')
        return strseqnum
    except:
        try:
            from Bio.Emboss.Applications import NeedleCommandline 
            print("Trying MSASearch with EMBOSS")
            output_handle = open('tmp/PDB_seq.fasta', 'w')
            SeqIO.write(SeqRecord(Seq(seq), id='PDB sequence'), output_handle, "fasta") 
            output_handle.close()
            output_handle = open("tmp/algn_seq.fasta", "w")
            s_records = list()
            for k in range(len(algn)):
                s_records.append(SeqRecord(Seq(algn[k]), id=str(k), description=hd[k]))
            SeqIO.write(s_records, output_handle, "fasta")
            output_handle.close()
            needle_cline = NeedleCommandline(path2_algprog+"needle",\
                asequence="tmp/PDB_seq.fasta",\
                bsequence="tmp/algn_seq.fasta", gapopen=10, gapextend=0.5, outfile="tmp/needle.txt")
            stdout, stderr = needle_cline()
            print(stdout+stderr)
            algres = open('tmp/needle.txt', 'r').readlines()
            score = list()
            for k in algres:
                if (k.find('Identity: ') > 0):
                    score.append(int(k.split()[2].split('/')[0]))
            i_0 = score.index(max(score))
            if species is not None:
                strseqnum = key_list[i_0]
            else:
                strseqnum = i_0
            shutil.rmtree('tmp')
            return strseqnum
        except:
            print("Trying MSASearch with BioPython")
            score = list()
            for k,s in enumerate(algn):
                score.append(pairwise2.align.globalxx(algn, s, one_alignment_only=1, score_only=1))
            i_0 = score.index(max(score)) 
            if species is not None:
                strseqnum = key_list[i_0]
            else:
                strseqnum = i_0 
            print("BP strseqnum is %i" % (strseqnum))
            return strseqnum

def chooseRefSeq(alg):
    ''' This function chooses a default reference sequence if none is given by taking the
    sequence which has the mean pairwise sequence identity closest to that of the entire alignment.
  
    :Example:
       >>> i_ref = chooseRefSeq(msa_num)'''
    
    if (len(alg) > 1000):
        seqw = seqWeights(alg)
        keep_seq =  randSel(seqw, 1000)
    else:
        keep_seq = [k for k in range(len(alg))]
    algNew = [alg[k] for k in keep_seq]
    numAlgNew = lett2num(algNew)
    simMat = seqSim(numAlgNew)
    listS = [simMat[i,j] for i in range(simMat.shape[0]) for j in range(i+1, simMat.shape[1])]
    meanSID = list()
    for k in range(len(simMat)):
        meanSID.append(simMat[k].mean())
    meanDiff = abs(meanSID - np.mean(listS))
    strseqnum = [i for i,k in enumerate(meanDiff) if k==min(meanDiff)]
    ix = keep_seq[strseqnum[0]]
    return strseqnum[0]

def makeATS(sequences, refpos, refseq, iref=0, truncate=False):
    ''' If specified, truncate the alignment to the structure (assumes MSAsearch_ has already been run
    to identify the reference sequence (iref)) and produce a mapping (ats) between alignment positions and the positions in the reference sequence (refpos).
    
    .. _MSAsearch: scaTools.html#scaTools.MSAsearch

     **Arguments:**
       -  sequences
       -  reference positions 
       -  reference sequence
       -  iref, the index of the sequence in the alignment with the highest identity to the reference

    :Keyword Arguments:
       truncate  -- truncate the alignment to positions that align to the reference sequence.

    :Example:
      >>> sequences_trun, ats_new = sca.makeATS(sequences_full, ats_pdb, seq_pdb, i_ref)

    '''
    if truncate == True:
        print("truncating to reference sequence...")
        # Removing gaps:
        pos_ref = [i for i,a in enumerate(refseq) if a != '-']
        seq_ref = ''.join([refseq[i] for i in pos_ref])
        ats_ref = [refpos[i] for i in pos_ref]
        pos_alg = [i for i,a in enumerate(sequences[iref]) if a != '-']
        seq_tr = [''.join([sq[i] for i in pos_alg]) for sq in sequences]
        # Positions to keep in the alignment and pbd sequences
        # (no gap in any of them after co-alignment):
        seqal_ref, seqal_alg, _, _, _ = pairwise2.align.globalms(seq_ref, seq_tr[iref],\
                                         2, -1, -.5, -.1)[0]
        keep_ref, keep_alg = list(), list()
        j_ref, j_alg = 0, 0
        for i in range(len(seqal_ref)):
            if seqal_ref[i] != '-' and seqal_alg[i] != '-':
                keep_ref.append(j_ref)
                keep_alg.append(j_alg)
            if seqal_ref[i] != '-': j_ref += 1
            if seqal_alg[i] != '-': j_alg += 1
        sequences_out = [''.join([sq[i] for i in keep_alg]) for sq in seq_tr]
        ats_out = [ats_ref[i] for i in keep_ref]
    else:
        tmp = sequences[iref].replace('-','.')
        refseq = refseq.replace('-','');
        seqal_ref, seqal_alg, _, _, _ = pairwise2.align.globalms(refseq, tmp,\
                                         2, -1, -.5, -.1)[0]
        print ('Len refseq %i, len refpos %i, Len alg seq %i, len pairalg %i, len gloalg %i' % (len(refseq),len(refpos), len(tmp),len(seqal_alg),len(sequences[0])))
        #print seqal_ref
        #print seqal_alg
        ats_out = list()
        j_ref = 0
        j_pdb = 0
        for i in range(len(seqal_alg)):
            if seqal_alg[i] == '.' and seqal_ref[i] == '-':
                ats_out.insert(j_ref,'-')
                j_ref += 1
            elif seqal_alg[i] != '.' and seqal_alg[i] != '-':
                if seqal_ref[i] != '-':
                    ats_out.insert(j_ref,refpos[j_pdb])
                    j_ref += 1
                    j_pdb += 1
                else:
                    ats_out.insert(j_ref,'-')
                    j_ref += 1
            elif seqal_alg[i] == '.' and seqal_ref[i] != '-':
                ats_out.insert(j_ref, refpos[j_pdb])
                j_ref += 1
                j_pdb += 1
            elif seqal_alg[i] == '-':
                j_pdb += 1
        sequences_out = sequences
    return sequences_out, ats_out

def lett2num(msa_lett, code='ACDEFGHIKLMNPQRSTVWY'):
    ''' Translate an alignment from a representation where the 20 natural amino
    acids are represented by letters to a representation where they are
    represented by the numbers 1,...,20, with any symbol not corresponding to an
    amino acid represented by 0.

    :Example:
       >>> msa_num = lett2num(msa_lett, code='ACDEFGHIKLMNPQRSTVWY') 

    '''
    lett2index = {aa:i+1 for i,aa in enumerate(code)}
    [Nseq, Npos] = [len(msa_lett), len(msa_lett[0])]
    msa_num = np.zeros((Nseq, Npos)).astype(int)
    for s, seq in enumerate(msa_lett):
        for i, lett in enumerate(seq):
            if lett in code:
                 msa_num[s, i] = lett2index[lett]
    return msa_num

def alg2bin(alg, N_aa=20):
    ''' Translate an alignment of size M x L where the amino acids are represented 
    by numbers between 0 and N_aa (obtained using lett2num) to a sparse binary 
    array of size M x (N_aa x L). 
    
    :Example:
      >>> Abin = alg2bin(alg, N_aa=20) '''
    
    [N_seq, N_pos] = alg.shape
    Abin_tens = np.zeros((N_aa, N_pos, N_seq))
    for ia in range(N_aa):
        Abin_tens[ia,:,:] = (alg == ia+1).T
    Abin = sparsify(Abin_tens.reshape(N_aa*N_pos, N_seq, order='F').T)
    return Abin
        
def seqWeights(alg, max_seqid=.8, gaps=1):
    ''' Compute sequence weights for an alignment (format: list of sequences)
    where the weight of a sequence is the inverse of the number of sequences in
    its neighborhood, defined as the sequences with sequence similarity below
    max_seqid. The sum of the weights defines an effective number of sequences.

    **Arguments:**
        -  `alg` = list of sequences

    **Keyword Arguments:**
        -  `max_seqid` 
        -  `gaps` = If gaps == 1 (default), considering gaps as a 21st amino acid, if gaps == 0, not considering them.
  
    :Example:
      >>> seqw = seqWeights(alg)    
  
    '''
    codeaa = 'ACDEFGHIKLMNPQRSTVWY'
    if gaps == 1: codeaa += '-'
    msa_num = lett2num(alg, code=codeaa)
    X2d = alg2bin(msa_num, N_aa=len(codeaa) )
    simMat = (X2d.dot(X2d.T)).todense()/msa_num.shape[1]
    seqw = np.array(1/(simMat>max_seqid).sum(axis=0))
    return seqw

def filterSeq(alg0, sref=0.5, max_fracgaps=.2, min_seqid=.2, max_seqid=.8):
    ''' Take in an alignment (alg0, assumed to be filtered to remove highly gapped positions),
    a reference sequence, the maximum fraction of gaps allowed per sequence (max_fracgaps),
    the minimum and maximum sequence identities to the reference sequence (min_seqid 
    and max_seqid), and return (1) alg, the alignment filtered to remove sequences with 
    more than max_fracgaps (i.e. partial seqs), (2) seqw, a vector of weights for each 
    sequence, (3) seqkeep, the indices of the original alignment (alg0) retained in alg:
    
    **Note:** if sref is set to 0.5, filterSeq calls chooseRefSeq_ to automatically select a
              reference sequence.

    .. _chooseRefSeq: scaTools.html#scaTools.chooseRefSeq
 
    **Example:**
        >>> alg, seqw, seqkeep = filterSeq(alg0, iref, max_fracgaps=.2, min_seqid=.2, max_seqid=.8) 
    
    '''
    if (sref == 0.5): sref = chooseRefSeq(alg0)
    Nseq, Npos = len(alg0), len(alg0[0])    
    # Elimination of sequences with too many gaps:
    seqkeep0 = [s for s in range(Nseq) if alg0[s].count('-')/Npos < max_fracgaps]
    print ("Keeping %i sequences of %i sequences (after filtering for gaps)"% (len(seqkeep0), Nseq))
    # Elimination of sequences too dissimilar to the reference (trimming):
    seqkeep = [s for s in seqkeep0 \
               if sum([alg0[s][i] == alg0[sref][i] for i in range(Npos)])/Npos > min_seqid]
    print ("Keeping %i sequences of %i sequences (after filtering for seq similarity)"% (len(seqkeep), len(seqkeep0)))
    alg = [alg0[s] for s in seqkeep]
    # Sequence weights (smoothing, here effectively treats gaps as a 21st amino acid):
    msa_num = lett2num(alg)
    seqw = seqWeights(alg, max_seqid)
    return alg, seqw, seqkeep
    
def filterPos(alg, seqw=[1], max_fracgaps=.2):
    ''' Truncate the positions of an input alignment to reduce gaps, taking into account sequence weights.

    **Arguments:**
        -  `alg` = An MxL list of sequences 
    
    **Keyword Arguments:**
        -  `seqw` = vector of sequence weights (default is uniform weights)
        -  `max_fracgaps` = maximum fraction of gaps allowed at a position

    **Returns:**
        -  `alg_tr` = the truncated alignment
        -  `selpos` = the index of retained positions (indices start at 0 for the first position)

    :Example:
       >>> alg_tr, selpos = filterPos(alg, seqw, max_fracgaps=.2) 

    '''
    Nseq, Npos = len(alg), len(alg[0])
    if len(seqw) == 1: seqw = np.tile(1, (1, Nseq))
    # Fraction of gaps, taking into account sequence weights:
    gapsMat = np.array([[int(alg[s][i]=='-') for i in range(Npos)] for s in range(Nseq)])
    seqwn = seqw/seqw.sum()
    gapsperpos = seqwn.dot(gapsMat)[0]
    # Selected positions:
    selpos = [i for i in range(Npos) if gapsperpos[i] < max_fracgaps]
    # Truncation:
    alg_tr = [''.join([alg[s][i] for i in selpos]) for s in range(Nseq)]
    return alg_tr, selpos     

def randSel(seqw, Mtot, keepSeq = []):
    ''' Random selection of Mtot sequences, drawn with weights and without replacement.
        The seed for the random number generator is fixed to ensure reproducibility.

        **Arguments:**
            -  `seqw` = the sequence weights
            -  `Mtot` = the total number of sequences for selection

        **Keyword Arguments:**
            -  `keepSeq` = an (optional) list of sequnces to keep. This can be useful if you would like to 
                           retain the reference sequence for example.

        :Example:
      >>> selection = randSel(seqw, Mtot, [iref]) 
    '''
    rand.seed(0)
    return weighted_rand_list(seqw[0], Mtot, keepSeq)

def weighted_rand_list(weights, Nmax, keepList):
    ''' Generate a random list of at most Nmax elements with weights (numpy array) but 
    without replacement. Called by randSel_.

    .. _randSel: scaTools.html#scaTools.randSel

    :Example:
      >>> selection = weighted_rand_list(weights, Nmax, [iref]) 

    '''
    Ntot = min((weights > 0).sum(), Nmax)
    wlist = [w for w in weights]
    selection = list()
    for k in keepList:
        selection.append(k)
        wlist[k] = 0
        Ntot -= 1
    for k in range(Ntot):
            i = weighted_rand_sel(wlist)
            selection.append(i)
            wlist[i] = 0
    return selection

def weighted_rand_sel(weights):
    ''' Generate a random index with probability given by input weights.
    Called by weighted_rand_list_.

    .. _weighted_rand_list: scaTools.html#scaTools.weighted_rand_list

    :Example:
      >>> index = weighted_rand_sel(weights) 

    '''
    rnd = rand.random() * sum(weights)
    for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                    return i

##########################################################################################
# BASIC STATISTICAL FUNCTIONS

def freq(alg, seqw=1, Naa=20, lbda=0, freq0=np.ones(20)/21):
    ''' 
    Compute amino acid frequencies for a given alignment.
    
    **Arguments:**
        -  `alg` = a MxL sequence alignment (converted using lett2num_) 
 
    .. _lett2num: scaTools.html#scaTools.lett2num

    **Keyword Arguments:**
        - `seqw` = a vector of sequence weights (1xM)
        - `Naa` = the number of amino acids
        - `lbda` = lambda parameter for setting the frequency of pseudo-counts (0 for no pseudo counts)
        - `freq0` = expected average frequency of amino acids at all positions

    **Returns:**
        -  `freq1` = the frequencies of amino acids at each position taken independently (Naa*L)
        -  `freq2` = the joint frequencies of amino acids at pairs of positions (freq2, Naa*L * Naa*L) 
        -  `freq0` = the average frequency of amino acids at all positions (Naa)

    :Example:
      >>> freq1, freq2, freq0 = freq(alg, seqw, lbda=lbda) 
   
    '''
    Nseq, Npos = alg.shape
    if type(seqw) == int and seqw == 1: seqw = np.ones((1,Nseq))
    seqwn = seqw/seqw.sum()
    al2d = alg2bin(alg, Naa)
    freq1 = seqwn.dot(np.array(al2d.todense()))[0]
    freq2 = np.array(al2d.T.dot(scipy.sparse.diags(seqwn[0], 0)).dot(al2d).todense())
    # Background:
    block = np.outer(freq0,freq0)
    freq2_bkg = np.zeros((Npos*Naa, Npos*Naa))
    for i in range(Npos): freq2_bkg[Naa*i:Naa*(i+1),Naa*i:Naa*(i+1)] = block
    # Regularizations:
    freq1_reg = (1-lbda)*freq1 + lbda*np.tile(freq0,Npos)
    freq2_reg = (1-lbda)*freq2 + lbda*freq2_bkg
    freq0_reg = freq1_reg.reshape(Npos, Naa).mean(axis=0)
    return freq1_reg, freq2_reg, freq0_reg

def eigenVect(M):
    ''' Return the eigenvectors and eigenvalues, ordered by decreasing values of the 
    eigenvalues, for a real symmetric matrix M. The sign of the eigenvectors is fixed
    so that the mean of its components is non-negative.

    :Example:
       >>> eigenVectors, eigenValues = eigenVect(M) 
    
    '''
    eigenValues,eigenVectors = np.linalg.eigh(M)
    idx = (-eigenValues).argsort()   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    for k in range(eigenVectors.shape[1]):
        if np.sign(np.mean(eigenVectors[:,k])) != 0:
            eigenVectors[:,k] = np.sign(np.mean(eigenVectors[:,k]))*eigenVectors[:,k]
    return eigenVectors, eigenValues

def svdss(X, k=6):
    ''' Singular value decomposition for sparse matrices (top k components). The singular
     values are ordered by decreasing values,  the sign of the singular vectors is fixed,
    and the convention is that X = u.dot(s).dot(v.T):
   
    :Example:
      >>> u, s ,v = svdss(X, k=6)

    '''
    u,s,vt = scipy.sparse.linalg.svds(X, k)
    idx = (-s).argsort()   
    s = s[idx]; u = u[:,idx]
    for j in range(u.shape[1]):
        sign = np.sign(np.mean(u[:,j]))
        u[:,j] = sign*u[:,j]
    v = X.T.dot(u).dot(np.diag(1/s))
    return u, s ,v

def basicICA(x, r, Niter):
    ''' Basic ICA algorithm, based on work by Bell & Sejnowski (infomax). The input data should preferentially be sphered, i.e., x.T.dot(x) = 1 

    **Arguments:**
      -  `x` = LxM input matrix where L = # features and M = # samples
      -  `r` = learning rate / relaxation parameter (e.g. r=.0001)
      -  `Niter` =  number of iterations (e.g. 1000)

    **Returns:**
      -  `w` = unmixing matrix
      -  `change` = record of incremental changes during the iterations.
       
    **Note:** r and Niter should be adjusted to achieve convergence, which should be assessed by visualizing 'change' with plot(range(iter) ,change)

    **Example:**
      >>> [w, change] = basicICA(x, r, Niter)
   
    '''
    [L, M] = x.shape
    w = np.eye(L)
    change = list()
    for _ in range(Niter):
        w_old = np.copy(w)
        u = w.dot(x)
        w += r*(M*np.eye(L)+(1-2*(1./(1+np.exp(-u)))).dot(u.T)).dot(w)
        delta = (w-w_old).ravel()
        change.append(delta.dot(delta.T))
    return [w, change]

def rotICA(V, kmax=6, learnrate=.0001, iterations=10000):
    ''' ICA rotation (using basicICA) with default parameters and normalization of 
    outputs. 
    
    :Example:
       >>> Vica, W = rotICA(V, kmax=6, learnrate=.0001, iterations=10000) 
    '''
    V1 = V[:,:kmax].T
    [W, changes_s] = basicICA(V1, learnrate, iterations)
    Vica = (W.dot(V1)).T
    for n in range(kmax):
        imax = abs(Vica[:,n]).argmax()
        Vica[:,n] = np.sign(Vica[imax,n])*Vica[:,n]/np.linalg.norm(Vica[:,n])
    return Vica, W

##########################################################################################
# SCA FUNCTIONS

def seqSim(alg):
    ''' Take an MxL alignment (converted to numeric representation using lett2num_) 
    and compute a MxM matrix of sequence similarities.
    
    :Example:
      >>> simMat = seqSim(alg)
 
    '''    
    # Get the number of sequences and number of positions:
    [Nseq, Npos] = alg.shape
    # Convert into a M*(20L) (sparse) binary array:
    X2d = alg2bin(alg)
    # Make the product with sparse matrices and convert it back to a dense array:
    simMat = (X2d.dot(X2d.T)).todense()/Npos
    return simMat

def posWeights(alg, seqw=1, lbda=0, N_aa = 20, freq0 = np.array([.073, .025, .050, .061, .042, .072,\
	.023, .053, .064, .089,.023, .043, .052, .040, .052, .073, .056, .063, .013, .033]), tolerance=1e-12):
    ''' Compute single-site measures of conservation, and the sca position weights, :math:`\\frac {\partial {D_i^a}}{\partial {f_i^a}}`

    **Arguments:**
         -  `alg` =  MSA, dimensions MxL, converted to numerical representation with lett2num_
         -  `seqw` = a vector of M sequence weights (default is uniform weighting)
         -  `lbda` = pseudo-counting frequencies, default is no pseudocounts
         -  `freq0` =  background amino acid frequencies :math:`q_i^a`
   
    **Returns:**
         -  `Wia` = positional weights from the derivation of a relative entropy, :math:`\\frac {\partial {D_i^a}}{\partial {f_i^a}}` (Lx20)
         -  `Dia` = the relative entropy per position and amino acid (Lx20)
         -   `Di` = the relative entropy per position (L) 

    :Example:
       >>> Wia, Dia, Di = posWeights(alg, seqw=1,freq0)
    
    '''
    N_seq, N_pos = alg.shape; 
    if type(seqw) == int and seqw == 1: seqw = np.ones((1,N_seq))
    freq1, freq2, _ = freq(alg, Naa=N_aa, seqw=seqw, lbda=lbda, freq0=freq0)
    # Overall fraction of gaps:
    theta = 1 - freq1.sum()/N_pos
    if theta<tolerance: theta = 0
    # Background frequencies with gaps:
    freqg0 = (1-theta)*freq0
    freq0v = np.tile(freq0,N_pos)
    iok = [i for i in range(N_pos*N_aa) if (freq1[i]>0 and freq1[i]<1)]
    # Derivatives of relative entropy per position and amino acid:
    Wia = np.zeros(N_pos*N_aa)
    Wia[iok] = abs(np.log((freq1[iok]*(1-freq0v[iok]))/((1-freq1[iok])*freq0v[iok])))
    # Relative entropies per position and amino acid:
    Dia = np.zeros(N_pos*N_aa)
    Dia[iok] = freq1[iok]*np.log(freq1[iok]/freq0v[iok])\
                   + (1-freq1[iok])*np.log((1-freq1[iok])/(1-freq0v[iok]))
    # Overall relative entropies per positions (taking gaps into account):
    Di = np.zeros(N_pos)
    for i in range(N_pos):
        freq1i = freq1[N_aa*i: N_aa*(i+1)]
        aok = [a for a in range(N_aa) if freq1i[a]>0]
        flogf = freq1i[aok]*np.log(freq1i[aok]/freqg0[aok])
        Di[i] = flogf.sum()
        freqgi = 1 - freq1i.sum()
        if freqgi > tolerance: Di[i] += freqgi*np.log(freqgi/theta)
    return Wia, Dia, Di

def seqProj(msa_num, seqw, kseq=15, kica=6):
    ''' Compute three different projections of the sequences based on eigenvectors of the sequence similarity matrix.
  
    **Arguments:**
       -  `msa_num` = sequence alignment (previously converted to numerical representation using lett2num_)
       -  `seqw` = a vector of sequence weights

    **Keyword Arguments:**
       -  `kseq` = number of eigenvectors to compute
       -  `kica` = number of independent components to compute

    **Returns:**
       -  `Useq[0]/Uica[0]` =  use no weight
       -  `Useq[1]/Uica[1]` =  use sequence weights
       -  `Useq[2]/Uica[2]` =  use sequence weights and positional weights
    
    **Example:**
      >>> Useq, Uica = sca.seqProj(msa_num, seqw, kseq = 30, kica = 15)
 
    '''
    posw, Dia, Di = posWeights(msa_num, seqw)
    Useq = list()
    # 1 - raw:
    X2d = alg2bin(msa_num)
    Useq.append(svdss(X2d,k=kseq)[0])
    # 2 - with sequence weights:
    X2dw = sparsify(np.diag(np.sqrt(seqw[0]))).dot(X2d)
    u, s, v = svdss(X2dw,k=kseq)
    Useq.append(X2d.dot(v).dot(np.diag(1/s)))
    # 3 - with sequence and position weights:
    X2dp = X2d.dot(sparsify(np.diag(posw)))
    X2dpw = sparsify(np.diag(np.sqrt(seqw[0]))).dot(X2dp)
    u, s, v = svdss(X2dpw,k=kseq)
    Useq.append(X2dp.dot(v).dot(np.diag(1/s)))
    # Fixing the sign:
    for U in Useq:
        for j in range(U.shape[1]):
            U[:,j] = np.sign(np.mean(U[:,j]))*U[:,j]
    # Rotation by ICA (default is kica=6):
    Uica = list()
    for U in Useq:
        Uica.append(rotICA(U, kmax=kica)[0]) 
    return Useq, Uica

def scaMat(alg, seqw=1, norm='frob',lbda=0, freq0=np.ones(20)/21,):
    ''' Computes the SCA matrix.
       
     **Arguments:**
        - `alg` =  A MxL multiple sequence alignment, converted to numeric representation with lett2num_
    
     **Keyword Arguments:**
        -  `seqw` =  A vector of sequence weights (default: uniform weights)
        -  `norm` =   The type of matrix norm used for dimension reduction of the
                      SCA correlation tensor to a positional correlation matrix.
                      Use 'spec' for spectral norm and 'frob' for Frobenius
                      norm.  The frobenius norm is the default.
        -  `lbda` =  lambda parameter for setting the frequency of pseudo-counts (0 for no pseudo counts) 
        -  `freq0` = background expectation for amino acid frequencies

     **Returns:**
        -  `Cp` = the LxL SCA positional correlation matrix
        -  `tX` = the projected MxL alignment
        -  `projMat` = the projector

     :Example:
      >>> Csca, tX, projMat = scaMat(alg, seqw, norm='frob', lbda=0.03)

    '''
    N_seq, N_pos = alg.shape; N_aa = 20
    if type(seqw) == int and seqw == 1: seqw = np.ones((1,N_seq)) 
    freq1, freq2, freq0 = freq(alg, Naa=N_aa, seqw=seqw, lbda=lbda, freq0=freq0)
    W_pos = posWeights(alg, seqw, lbda)[0]
    tildeC = np.outer(W_pos,W_pos)*(freq2 - np.outer(freq1,freq1))
    # Positional correlations:
    Cspec = np.zeros((N_pos,N_pos))
    Cfrob = np.zeros((N_pos,N_pos))
    P = np.zeros((N_pos,N_pos,N_aa))
    for i in range(N_pos):
        for j in range(i,N_pos):
            u,s,vt = np.linalg.svd(tildeC[N_aa*i:N_aa*(i+1), N_aa*j:N_aa*(j+1)])
            Cspec[i,j] = s[0]
            Cfrob[i,j] = np.sqrt(sum(s**2))
            P[i,j,:] = np.sign(np.mean(u[:,0]))*u[:,0]
            P[j,i,:] = np.sign(np.mean(u[:,0]))*vt[0,:].T
    Cspec += np.triu(Cspec,1).T
    Cfrob += np.triu(Cfrob,1).T
    # Projector:
    al2d = np.array(alg2bin(alg).todense())
    tX = np.zeros((N_seq,N_pos))
    Proj = W_pos*freq1
    ProjMat = np.zeros((N_pos,N_aa))
    for i in range(N_pos):
        Projati = Proj[N_aa*i:N_aa*(i+1)] 
        if sum(Projati**2) > 0:
            Projati /= np.sqrt(sum(Projati**2))
        ProjMat[i,:] = Projati
        tX[:,i] = al2d[:,N_aa*i:N_aa*(i+1)].dot(Projati.T)
    if norm == 'frob' : Cspec = Cfrob
    return Cspec, tX, Proj

##########################################################################################
# PROJECTIONS OF ANNOATED SEQUENCES

def projUica(msa_ann, msa_num, seqw, kica=6):
    ''' Compute the projection of an alignment (msa_ann) on the kpos ICA components
    of the sequence space of another (msa_num, seqw).This is useful to compare the sequence space of one alignment to another.

    :Example:
      >>> Uica_ann, Uica = projUpica(msa_ann, msa_num_ seqw, kica=6) 
   
    '''
    X2d = alg2bin(msa_num)
    posw, Dia, Di = posWeights(msa_num, seqw)
    X2dw = sparsify(np.diag(np.sqrt(seqw[0]))).dot(X2d)
    u, s, v = svdss(X2dw,k=kica)
    P = v.dot(np.diag(1/s))
    U = X2d.dot(P)
    for j in range(U.shape[1]):
        P[:,j] = np.sign(np.mean(U[:,j]))*P[:,j]
        U[:,j] = np.sign(np.mean(U[:,j]))*U[:,j]
    U, W = rotICA(U, kmax=kica)
    X2d_new = alg2bin(msa_ann)
    U0 = X2d.dot(P)
    U1 = X2d_new.dot(P)
    Ui0 = (W.dot(U0[:,:kica].T)).T
    Ui1 = (W.dot(U1[:,:kica].T)).T
    for n in range(kica):
        imax = abs(Ui0[:,n]).argmax()
        Ui1[:,n] = np.sign(Ui0[imax,n])*Ui1[:,n]/np.linalg.norm(Ui0[:,n])
        Ui0[:,n] = np.sign(Ui0[imax,n])*Ui0[:,n]/np.linalg.norm(Ui0[:,n])
    return Ui1, Ui0

def projAlg(alg, Proj):
    ''' Projection of an alignment (alg) based on a projector (Proj). The input alignment should already be converted to numeric representation using lett2num_.

    :Example:
      >>> tX = projAlg(msa_num, Proj) 

    '''
    N_seq, N_pos = alg.shape; N_aa = 20
    al2d = np.array(alg2bin(alg).todense())
    tX = np.zeros((N_seq,N_pos))
    ProjMat = np.zeros((N_pos,N_aa))
    for i in range(N_pos):
        Projati = Proj[N_aa*i:N_aa*(i+1)] 
        if sum(Projati**2) > 0:
            Projati /= np.sqrt(sum(Projati**2))
        ProjMat[i,:] = Projati
        tX[:,i] = al2d[:,N_aa*i:N_aa*(i+1)].dot(Projati.T)
    return tX

def projUpica(msa_ann, msa_num, seqw, kpos):
    ''' Compute the projection of an alignment (msa_ann) on the kpos ICA components
    of the SCA matrix of another (msa_num, seqw). This is useful to compare the sequence space (as projected by the positional correlations) of one alignment to another. 

    :Example:
      >>> Upica_ann, Upica = projUpica(msa_ann, msa_num_ seqw, kpos) 

    '''
    Csca, tX, Proj = scaMat(msa_num, seqw)
    Vsca, Lsca = eigenVect(Csca)
    Vpica, Wpica = rotICA(Vsca, kmax=kpos)
    Usca = tX.dot(Vsca[:,:kpos]).dot(np.diag(1/np.sqrt(Lsca[:kpos])))
    tX_ann = projAlg(msa_ann, Proj)
    Usca_ann = tX_ann.dot(Vsca[:,:kpos]).dot(np.diag(1/np.sqrt(Lsca[:kpos])))
    Upica = Wpica.dot(Usca.T).T
    Upica_ann = Wpica.dot(Usca_ann.T).T
    for k in range(Upica.shape[1]): 
        Upica_ann[:,k] /= np.sqrt(Upica[:,k].T.dot(Upica[:,k]))
        Upica[:,k] /= np.sqrt(Upica[:,k].T.dot(Upica[:,k]))        
    return Upica_ann, Upica


##########################################################################################
# SECTOR ANALYSIS

def sizeLargestCompo(adjMat):
    ''' Compute the size of the largest component of a graph given its adjacency matrix. 
    Called by numConnected_ (Done by actually listing all the components)

    .. _numConnected: scaTools.html#scaTools.sizeLargestCompo

    :Example:
      >>> s = sizeLargestCompo(adjMat)
    
    '''
    Nnodes = adjMat.shape[0]
    found = list()
    components = list()
    for node in range(Nnodes):
        # If the node has not yet been encountered, 
        # it belongs to a new component:
        if node not in found:
            newcomponent = [node]
            found.append(node)
            i = 0
            # Recursively listing the neighboors:
            while (i < len(newcomponent)):
                newneighbors = [j for j in range(Nnodes)\
                    if (adjMat[newcomponent[i],j] == 1 and j not in newcomponent)]
                newcomponent += newneighbors
                found += newneighbors
                i += 1
            # Adding the new component to the list of all components:
            components.append(newcomponent)
    # Returning only the size of the maximal component:
    return max([len(compo) for compo in components])
    
def numConnected(Vp, k, distmat, eps_list=np.arange(.5,0,-.01), dcontact=5):
    ''' Calculates the number of positions in the largest connected component for groups of positions i
     with :math:`V_p[i,k] > eps` and :math:`V_p[i,k] > V_p[i,kk]`, for :math:`kk != k` and eps in eps_list. Useful for looking evaluating the physical connectivity of different sectors or sub-sectors.
    
    **Arguments**: 
       -  `Vp` = A set of eigenvectors or independent components
       -  `k`  = The eigenvector or independent component to consider
       -  `distmat` = Distance matrix (computed by pdbSeq_)

    .. _pdbSeq:  scaTools.html#scaTools.pdbSeq
    
    **Keyword Arguments:**
        -  `eps_list` = the range of values of eps for which the group is non-empty
        -  `dcontact` = the distance cutoff for defining contacts

    **Example:**
       >>> eps_range, num_co, num_tot = numConnected(Vp, k, distmat, eps_list = np.arange(.5,0,-.01), dcontact=8)
    
    '''
    Npos = distmat.shape[0]
    eps_range = list()
    num_co = list()
    num_tot = list()
    [Npos, kmax] = Vp.shape
    for eps in eps_list:
        Vothers = [0 for i in range(Npos)]
        for kk in range(kmax):
            if kk != k:
                Vothers = [max(Vothers[i], Vp[i,kk]) for i in range(Npos)]
        group = [i for i in range(Npos) if (Vp[i,k]>eps and Vp[i,k]>Vothers[i])]
        eps_range.append(eps)
        num_tot.append(len(group))
        if len(group) > 0:
            adjMat = distmat[np.ix_(group,group)] < dcontact
            num_co.append(sizeLargestCompo(adjMat))
        else:
            num_co.append(0)                  
    return eps_range, num_co, num_tot

def chooseKpos(Lsca,Lrand):
    ''' Given the eigenvalues of the sca matrix (Lsca), and the eigenvalues for the set of randomized matrices (Lrand), return the number of significant eigenmodes.'''
    return Lsca[Lsca>(Lrand[:,1].mean()+(3*Lrand[:,1].std()))].shape[0]

def icList(Vpica, kpos, Csca, p_cut=0.95):
    ''' Produces a list of positions contributing to each independent component (IC) above
    a defined statistical cutoff (p_cut, the cutoff on the CDF of the t-distribution
    fit to the histogram of each IC).  Any position above the cutoff on more than one IC 
    are assigned to one IC based on which group of positions to which it shows a higher
    degree of coevolution. Additionally returns the numeric value of the cutoff for each IC, and the
    pdf fit, which can be used for plotting/evaluation.
    icList, icsize, sortedpos, cutoff, pd  = icList(Vsca,Lsca,Lrand) '''
    #do the PDF/CDF fit, and assign cutoffs
    Npos = len(Vpica)
    cutoff = list()
    scaled_pdf = list()
    all_fits = list()
    for k in range(kpos):
        pd = t.fit(Vpica[:,k])
        all_fits.append(pd)
        iqr = scoreatpercentile(Vpica[:,k],75) - scoreatpercentile(Vpica[:,k],25)
        binwidth=2*iqr*(len(Vpica[:,k])**(-0.33))
        nbins=round((max(Vpica[:,k])-min(Vpica[:,k]))/binwidth)
        h_params = np.histogram(Vpica[:,k], int(nbins))
        x_dist = np.linspace(min(h_params[1]), max(h_params[1]), num=100)
        area_hist=Npos*(h_params[1][2]-h_params[1][1]);
        scaled_pdf.append(area_hist*(t.pdf(x_dist,pd[0],pd[1],pd[2])))
        cd = t.cdf(x_dist,pd[0],pd[1],pd[2])
        tmp = scaled_pdf[k].argmax()
        if abs(max(Vpica[:,k])) > abs(min(Vpica[:,k])):
            tail = cd[tmp:len(cd)]
        else:
            cd = 1 - cd
            tail = cd[0:tmp]
        diff = abs(tail - p_cut);
        x_pos = diff.argmin()
        cutoff.append(x_dist[x_pos+tmp])
    #select the positions with significant contributions to each IC
    ic_init = list()
    for k in range(kpos): ic_init.append([i for i in range(Npos) if Vpica[i,k]> cutoff[k]])
    #construct the sorted, non-redundant iclist    
    sortedpos = list()
    icsize = list()
    ics = list()
    icpos_tmp = list()
    Csca_nodiag = Csca.copy()
    for i in range(Npos): Csca_nodiag[i,i]=0
    for k in range(kpos):
        icpos_tmp = list(ic_init[k])
        for kprime in [kp for kp in range(kpos) if (kp != k)]:
            tmp = [v for v in icpos_tmp if v in ic_init[kprime]] 
            for i in tmp:
                remsec = np.linalg.norm(Csca_nodiag[i,ic_init[k]]) \
                         < np.linalg.norm(Csca_nodiag[i,ic_init[kprime]])
                if remsec: icpos_tmp.remove(i)
        sortedpos += sorted(icpos_tmp, key=lambda i: -Vpica[i,k])
        icsize.append(len(icpos_tmp))  
        s = Unit()
        s.items = sorted(icpos_tmp, key=lambda i: -Vpica[i,k])
        s.col = k/kpos
        s.vect = -Vpica[s.items,k]
        ics.append(s)
    return ics, icsize, sortedpos, cutoff, scaled_pdf, all_fits

def singleBar(x, loc, cols, width=.5):
    ''' Single bar diagram, called by MultiBar_.
       
    .. _MultiBar: scaTools.html#scaTools.MultiBar
    
    :Example:
      >>> singleBar(x, loc, cols, width=.5) 

    '''
    s = 0; y = list()
    for i,v in enumerate(x):
        s += v; y.append(s)
    for i,v in enumerate(x):
        plt.bar(loc, y[len(x)-i-1], width, color=cols[len(x)-i-1])
        
def MultiBar(x, colors = 'wbrgymc', width=.5):
    ''' Multiple bar diagram (plots contributions to each bar from different elements in x as different colors). This can be useful if you'd like to inspect how sector positions are distributed among independent components/eigenmodes. The argument x is a tuple, specifying the number of elements in each bar to be each color. The example below makes a graph with four bars, where the first bar has 99 white elements, 1 blue element and 1 red element.

    :Example: 
     >>> x = [[99, 1, 1], [6, 13, 2], [0, 0, 13], [1, 7, 5]]
     >>> sca.MultiBar(x)
'''
    for i,v in enumerate(x):
        singleBar(v,i,cols=colors)
    plt.xticks(np.arange(len(x)) + width/2., ['S%i' % i for i in range(len(x))],\
               fontsize=14)
    plt.axis([-width/2, len(x)-width/2, 0, max([sum(v) for v in x])+5]);

##########################################################################################
# DIRECT COUPLING ANALYSIS (DCA)

class Pair:
    ''' A class for a pair of positions. 
        :Attributes:
 
           -  `pos` = a pair of amino acid positions (ex: [1,3], supplied as argument p)
           -  `DI`  = the direct information between the two positions (argument x)
           -  `dist` = the physical distance between the positions (argument d)
    '''
    def __init__(self, p, x, d):
        self.pos = p
        self.DI = x
        self.dist = d
        
def directInfo(freq1, freq2, lbda=.5, freq0=np.ones(20)/21, Naa=20):
    ''' Calculate direct information as in the Direct Coupling Analysis (DCA) method proposed by 
    M. Weigt et collaborators (Ref: Marcos et al, PNAS 2011, 108: E1293-E1301).

    :Example:
      >>> DI = directInfo(freq1, freq2, lbda=.5, freq0=np.ones(20)/21, Naa=20)
    '''
    Npos = int(len(freq1)/Naa)
    Cmat_dat = freq2 - np.outer(freq1,freq1)
    # Background:
    block = np.diag(freq0) - np.outer(freq0,freq0)
    Cmat_bkg = np.zeros((Npos*Naa, Npos*Naa))
    for i in range(Npos): Cmat_bkg[Naa*i:Naa*(i+1),Naa*i:Naa*(i+1)] = block
    # Regularizations:
    Cmat = (1-lbda)*Cmat_dat + lbda*Cmat_bkg
    frq = (1-lbda)*freq1 + lbda*np.tile(freq0,Npos)
    # DI at mean-field approx:
    Jmat = -np.linalg.inv(Cmat)
    DI = np.zeros((Npos, Npos))
    for i in range(Npos):
        for j in range(i+1,Npos):
            DI[i,j] = dirInfoFromJ(i, j, Jmat, frq, Naa)
            DI[j,i] = DI[i,j]
    return DI

def dirInfoFromJ(i, j, Jmat, frq, Naa=20, epsilon=1e-4):
    ''' Direct information from the matrix of couplings :math:`J_{ij}` (called by directInfo_). Ref: Marcos et al, PNAS 2011, 108: E1293-E1301

    .. _directInfo: scaTools.html#scaTools.directInfo
    
    **Arguments:**
        - `i` = position 1
        - `j` = position 2
        - `Jmat` = coupling matrix
        - `frq` = frequency
        - `Naa` = number of amino acids

    :Example:
      >>> DI = dirInfoFromJ(i, j, Jmat, frq, Naa=20, epsilon=1e-4)
 
   '''
    W = np.ones((Naa+1,Naa+1))
    W[:Naa,:Naa] = np.exp(Jmat[Naa*i:Naa*(i+1), Naa*j:Naa*(j+1)])
    mui = np.ones(Naa+1)/(Naa+1)
    muj = np.ones(Naa+1)/(Naa+1)
    pi = np.zeros(Naa+1); pi[:Naa] = frq[Naa*i:Naa*(i+1)]; pi[Naa] = 1 - sum(pi)
    pj = np.zeros(Naa+1); pj[:Naa] = frq[Naa*j:Naa*(j+1)]; pj[Naa] = 1 - sum(pj)
    diff = epsilon + 1
    while (diff > epsilon):
        scrai = muj.dot(W.T); scraj = mui.dot(W)
        newi = pi/scrai; newi /= newi.sum()
        newj = pj/scraj; newj /= newj.sum()
        diff = max(abs(newi-mui).max(), abs(newj-muj).max())
        mui = newi; muj = newj
    Pdir = W*np.outer(mui, muj)
    Pdir /= Pdir.sum()
    Pfac = np.outer(pi, pj)
    tiny = 1e-100
    return np.trace(Pdir.T.dot(np.log((Pdir+tiny)/(Pfac+tiny))))

def truncDiag(M, dmax):
    ''' Set to 0 the elements of a matrix M up to a distance dmax from the diagonal.

    :Example:
       >>> Mtr = truncDiag(M, dmax)

    '''
    Mtr = copy.copy(M)
    for i in range(M.shape[0]):
        for j in range(i,min(i+dmax+1,M.shape[1])):
            Mtr[i,j] = 0; Mtr[j,i] = 0
    return Mtr

class Secton:
    ''' A class for sectons. 
  
        **Attributes:**
           -  `pos` = a list of positions
           -  `num` = number of positions in the secton
    '''
    def __init__(self, positions):
        self.pos = positions
        self.num = len(positions)

    def dist(self, distmat):
        ''' returns the distance between the position pair''' 
        return distmat[np.ix_(self.pos, self.pos)]

    def connected(self, distmat, threshold):
        ''' Check the structural connectivity based on the principle that if :math:`M_{ij}` is the adjacency matrix of a graph, :math:`M^n_{ij}` is the number of paths of length :math:`n` between i and j, which must be > 0 for :math:`n` = number of nodes when i and j are in the same connected component. '''
        return (np.linalg.matrix_power(self.dist(distmat)< threshold, self.num) > 0)\
        .sum()/self.num**2 == 1

##########################################################################################
# RANDOMIZATION
def preprocess_frq(frq):
    # 将小于0的值设置为0，将大于1的值设置为1，将NaN值设置为均匀分布的概率
    frq = np.clip(frq, 0, 1)
    nan_mask = np.isnan(frq)
    frq[nan_mask] = 1.0 / frq.shape[1]
    return frq
def randAlg(frq, Mseq):
    ''' Generate a random alignment with Mseq sequences based on the
    frequencies frq[i,a] of amino acids with a = 0,1,...,Naa (0 for gaps).
    
    :Example:
       >>> msa_rand = randAlg(frq, Mseq) 

    '''
    # 预处理频率矩阵 frq
    frq = preprocess_frq(frq)

    Npos = frq.shape[0]
    msa_rand = np.zeros((Mseq, Npos), dtype=int)
    for i in range(Npos):
        Maa = np.random.multinomial(Mseq, frq[i,:])
        col = np.array([], dtype=int)
        for aa, M in enumerate(Maa):
            col = np.append(col, np.tile(aa, M))
        np.random.shuffle(col)
        msa_rand[:, i] = col
    return msa_rand

def randomize(msa_num, Ntrials, seqw=1, norm='frob', lbda=0, Naa=20, kmax=6):
    ''' Randomize the alignment while preserving the frequencies of amino acids at each 
    position and compute the resulting spectrum of the SCA matrix.

    **Arguments:**
        -  `msa_num` = a MxL sequence alignment (converted to numerical representation using lett2num_)
        -  `Ntrials` = number of trials
        -  `seqw` = vector of sequence weights (default is to assume equal weighting)
        -  `norm` = either 'frob' (frobenius) or 'spec' (spectral)
        -  `lbda` = lambda parameter for setting the frequency of pseudo-counts (0 for no pseudo counts)
        -  `Naa` = number of amino acids
        -  `kmax` = number of eigenvectors to keep for each randomized trial
    
    **Returns:**
        -  `Vrand` =  eigenvectors for the :math:`\\tilde {C_{ij}^{ab}}` matrix of the randomized alignment (dimensions: Ntrials*Npos*kmax)
        -  `Lrand` =  eigenvalues for the :math:`\\tilde {C_{ij}^{ab}}` matrix of the randomized alignment  (dimensions: Ntrials*Npos)

    :Example:
       >>> Vrand, Lrand, Crand = randomize(msa_num, 10, seqw, Naa=20, kmax=6)

    '''
    if type(seqw) == int and seqw == 1: seqw = np.ones((1,Nseq))
    Mseq = np.round(seqw.sum()).astype(int)
    Nseq, Npos = msa_num.shape
    Crnd = np.zeros((Npos,Npos))
    # Weighted frequencies, including gaps:
    f1, f2, f0 = freq(msa_num, Naa=20, seqw=seqw, lbda=lbda, freq0=np.ones(20)/21)
    fr1 = np.reshape(f1,(Npos,Naa))
    fr0 = (1-fr1.sum(axis=1)).reshape(Npos,1)
    fr01 = np.concatenate((fr0,fr1),axis=1)
    # Multiple randomizations:
    Vrand = np.zeros((Ntrials,Npos,kmax))
    Lrand = np.zeros((Ntrials,Npos))
    for t in range(Ntrials):
        msa_rand = randAlg(fr01, Mseq)
        Csca = scaMat(msa_rand,norm=norm,lbda=lbda)[0]
        Crnd += Csca
        V, L = eigenVect(Csca)
        Vrand[t,:,:] = V[:,:kmax]
        Lrand[t,:] = L
    Crnd = Crnd/Ntrials
    return Vrand, Lrand, Crnd

##########################################################################################
# DISPLAY

def figWeights(U1, U2, weight):
    ''' A 2d scatter plot with color indicating weight.
        
        :Example:
       >>> figWeights(U1, U2, weight) '''
    seqcol = -np.log(weight)
    seqcol = (seqcol-min(seqcol)+1)/(max(seqcol)-min(seqcol)+1)
    seqorder = sorted(range(len(seqcol)), key=lambda s: seqcol[s])
    for s in seqorder:
        plt.plot(U1[s], U2[s], 'o', color = cm.jet(seqcol[s],1),\
             markeredgecolor = 'none')
       
def figColors():
    ''' Color code for figUnits_.
      
        .. _figUnits: scaTools.html#scaTools.figUnits 

        :Example:
      >>> figColors()

    '''
    plt.rcParams['figure.figsize'] = 6, 6
    for s in np.arange(0,1,.05):
        for a in np.arange(0,1,.01):
            bgr = colorsys.hsv_to_rgb(a,s,1)
            plt.plot(s*np.cos(2*np.pi*a), s*np.sin(2*np.pi*a),'o', markersize=8,\
            markerfacecolor=bgr, markeredgecolor=bgr)
    plt.title(r'Color at angle $\alpha$ encoded with $\alpha/(2\pi)$.', fontsize=16)
    plt.axis([-1.1, 1.1, -1.1, 1.1]);
    
def figUnits(v1, v2, units, marker='o', dotsize=9, notinunits=1):
    ''' 2d scatter plot specified by 'units', which must be a list of elements
    in the class Unit_. See figColors_ for the color code. Admissible color codes are in [0 1] 
    (light/dark gray can also be obtained by using -1/+1).
    For instance: 0->red, 1/3->green, 2/3-> blue.
    
    .. _Unit: scaTools.html#scaTools.Unit
    .. _figColors: scaTools.html#scaTools.figColors

    **Arguments:**
       -  `v1` = xvals
       -  `v2` = yvals
       -  `units` = list of elements in units
    
    **Keyword Arguments**
       -  `marker` = plot marker symbol
       -  `dotsize` = specify marker/dotsize
       -  `notinunits` = if set to 1 : the elements not in a unit are represented in white, if set to 0 
                         these elements are not represented, if set to [w1,w2] : elements with coordinates 
                         w1,w2 are represented in white in the background.

    :Example:
      >>> figUnits(v1, v2, units, marker='o', gradcol=0, dotsize=9, notinunits=1)
    
     ''' 
    Ntot = len(v1)
    # Plot all items in white:
    if notinunits == 1:
        plt.plot(v1,v2,marker, markersize=dotsize, markerfacecolor='w',\
         markeredgecolor='k')
    elif len(notinunits) == 2:
        plt.plot(notinunits[0],notinunits[1],marker, markersize=dotsize,\
        markerfacecolor='w', markeredgecolor='k')
    # Plot items in the units with colors:
    for u in units:
        items_list = list(u.items)
        if u.col>=0 and u.col<1: bgr = colorsys.hsv_to_rgb(u.col,1,1)
        if u.col == 1: bgr = [.3,.3,.3]
        if u.col < 0: bgr = [.7,.7,.7]
        plt.plot(v1[np.ix_(items_list)],v2[np.ix_(items_list)],marker,\
                 markersize=dotsize, markerfacecolor=bgr, markeredgecolor='k')
                                                                       
def figMapping(Csca, tX, kpos, sectors, subfam):
    ''' Function that automates finding the top :math:`k_{pos}` independent components, projection, and plotting.
        Useful to get a representation of the sectors/subfamilies mapping for a given kpos.
        **Arguments:**
            -  `Csca` = the sca matrix
            -  `tX` = the projected alignment
            -  `kpos` = the number of independent components to consider
            -  `sectors` = list of Unit_ elements for each sector
            -  `subfam` = list of Unit_ elements for each sequence family
        
        **Returns:**
            -  `Vpica` = the independent components of Csca

        :Example: 
        >>> Vpica = figMapping(Csca, tX, kpos, sectors, subfam) 

        '''
    Vsca, Lsca = eigenVect(Csca)
    Vpica, Wpica = rotICA(Vsca, kmax=kpos)
    Usca = tX.dot(Vsca[:,:kpos]).dot(np.diag(1/np.sqrt(Lsca[:kpos])))
    Upica = Wpica.dot(Usca.T).T
    for k in range(Upica.shape[1]): Upica[:,k] /= np.sqrt(Upica[:,k].T.dot(Upica[:,k]))
    Usica, Wsica = rotICA(Usca, kmax=kpos)   
    pairs = [[i,i+1] for i in [j for j in range(kpos-1) if (j%2==0)]]
    if (kpos%2==1): pairs.append([kpos-1,0])
    plt.rcParams['figure.figsize'] = (4*len(pairs))+1, 8 
    for n,[k1,k2] in enumerate(pairs):
        plt.subplot(2,len(pairs),n+1)
        figUnits(Vpica[:,k1], Vpica[:,k2], sectors)
        plt.xlabel(r'$V^p_{%i}$' % (k1+1), fontsize=16)
        plt.ylabel(r'$V^p_{%i}$' % (k2+1), fontsize=16)
        plt.subplot(2,len(pairs),n+len(pairs)+1)
        figUnits(Upica[:,k1], Upica[:,k2], subfam, marker='D')
        plt.xlabel(r'$U^p_{%i}$' % (k1+1), fontsize=16)
        plt.ylabel(r'$U^p_{%i}$' % (k2+1), fontsize=16)
    plt.tight_layout()
    return Vpica 

##########################################################################################
# PDB PROCESSING

def pdbSeq(pdbid, chain='A', path2pdb=path2structures, calcDist=1):
    ''' Extract sequence, position labels and matrix of distances from a PDB file.
    
    **Arguments:**
       -  `pdbid` = PDB identifier (four letters/numbers)
       -  `chain` = PDB chain identifier
       -  `path2pdb` = location of the PDB file
       -  `calcDist` = calculate a distance matrix between all pairs of positions, default is 1

    :Example:
       >>> sequence, labels, dist = pdbSeq(pdbid, chain='A', path2pdb=path2structures) '''
    # Table of 3-letter to 1-letter code for amino acids
    aatable = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',\
     'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',\
     'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    # Read PDB structure:
    structure = PDBParser().get_structure(pdbid, path2pdb+pdbid+'.pdb')
    # Fill up sequence and label information
    sequence = ''; labels = list()
    residues = [res for res in structure[0][chain] if res.get_id()[0] == ' ']
    for res in residues:
        labels.append(str(res.get_id()[1])+str(res.get_id()[2]).strip())
        try:
            sequence += aatable[res.get_resname()]
        except:
            sequence += 'X'
    # Distances between residues (minimal distance between atoms, in angstrom): 
    dist = np.zeros((len(residues), len(residues)))
    if (calcDist == 1):
        for n0, res0 in enumerate(residues):
            for n1, res1 in enumerate(residues):
                dist[n0,n1] = min([atom0 - atom1 for atom0 in res0 for atom1 in res1])
        return sequence, labels, dist
    else:
        return sequence, labels

def writePymol(pdb, sectors, ics, ats, outfilename, chain='A',inpath=path2structures, quit=1):
    ''' Write basic a pymol script for displaying sectors and exporting an image.

    :Example: 
      >>> writePymol(pdb, sectors, ics, ats, outfilename, chain='A',inpath=path2structures, quit=1)
 
    '''
    f = open(outfilename, 'w')
    f.write('delete all\n')
    f.write('load %s%s.pdb, main\n' % (inpath, pdb))
    f.write('hide all\n')
    f.write('bg_color white\n')
    f.write('show cartoon, (chain %s)\n' % chain)
    f.write('color white\n\n')
    for k,sec in enumerate(sectors):       
        b,g,r = colorsys.hsv_to_rgb(sec.col,1,1)
        f.write('set_color color%i, [%.3f,%.3f,%.3f]\n' % (k+1,b,g,r))
        f.write('create sector%i, (resi %s) & (chain %s)\n' \
                    % (k+1,','.join([ats[s] for s in sec.items]),chain))
        f.write('color color%i, sector%i\n' % (k+1,k+1))
        f.write('show spheres, sector%i\n' % (k+1))
        f.write('show surface, sector%i\n\n' % (k+1))       
    for k,sec in enumerate(ics):       
        b,g,r = colorsys.hsv_to_rgb(sec.col,1,1)
        f.write('set_color color_ic%i, [%.3f,%.3f,%.3f]\n' % (k+1,b,g,r))
        f.write('create ic_%i, (resi %s) & (chain %s)\n' \
                    % (k+1,','.join([ats[s] for s in sec.items]),chain))
        f.write('color color_ic%i, ic_%i\n' % (k+1,k+1))
        f.write('show spheres, ic_%i\n' % (k+1))
        f.write('show surface, ic_%i\n\n' % (k+1))          
    f.write('zoom\n')
    f.write('set transparency, 0.4\n')
    f.write('ray\n')
    path_list = outfilename.split(os.sep)
    fn = path_list[-1]
    f.write('png %s\n' % fn.replace('.pml',''))
    if quit==1:
        f.write('quit')
    f.close()    
    
def figStruct(pdbid, sectors, ats, chainid='A', outfile = 'Outputs/sectors.pml',\
              quit=1, pymol = path2pymol):
    ''' Make and display an image of the sectors (within a python notebook). By default quit PyMol after running it,
    unless the option 'quit=0' is given. The default name and location of the output can 
    also be changed.

    :Example:
     >>> figStruct(pdbid, sectors, ats, chainid='A', outfile = 'Outputs/sectors.pml',\
              quit=1, pymol = path2pymol) 

    '''
    writePymol(pdbid, sectors, ats, outfile, chain=chainid, quit=1)
    os.system(pymol+' '+outfile)
    img = mpimg.imread(outfile.replace('.pml','.png'))
    plt.imshow(img)
    plt.axis('off')

##########################################################################################
# CYTOSCAPE OUTPUT

def cytoscapeOut(ats, cutoff, Csca, Di, sectors, Vp, outfilename):
    ''' Output tab-delimited text that can be read in by cytoscape. The goal is to enable
    graph representations of the SCA couplings, where residues are nodes, and couplings are edges.
    Within cytoscape, the graph can be color-coded or weighted by Csca, Di, sector definition or Vp. 
    
    :Example:
      >>> cytoscapeOut(ats, cutoff, Csca, Di, sectors, Vp, outfilename)'''
    f = open(outfilename+'.sif', 'w')
    for k in range(len(ats)):
        flag = 0;
        for j in range(k+1,len(ats)):
            if (Csca[k][j] > cutoff):
                f.write(ats[k] + ' aa ' + ats[j] + '\n')
                flag = 1;
        if flag == 0:
            f.write(ats[k] + '\n')
    f.close()
    
    f = open(outfilename+'.eda', 'w')
    f.write('KEY\tSCA\n')
    for k in range(len(ats)):
        for j in range(k+1,len(ats)):
            f.write((ats[k] + ' (aa) ' + ats[j] + '\t  %.4f \n') % Csca[k][j])
    f.close()
    
    s_idx = [0 for k in range(len(ats))]
    for i,j in enumerate(sectors):
        for k in j.items: 
                s_idx[k] = i+1 

    f = open(outfilename+'.noa', 'w')
    f.write('KEY\tCONSERVATION\tSector\tVp1\tVp2\tVp3\n')
    for j,k in enumerate(ats):
        f.write((k + '\t %.4f \t %i \t %.4f \t %.4f \t %.4f \n') %(Di[j], s_idx[j], Vp[j,0], Vp[j,1], Vp[j,2]))
    f.close()
 
