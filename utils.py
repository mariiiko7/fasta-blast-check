import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from glob import glob
from Bio.SearchIO._legacy import NCBIStandalone
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def FASTA(filename):
    '''
    FASTA parser found in http://www.petercollingridge.co.uk/tutorials/bioinformatics/fasta-parser/
    Parser FASTA file and return the list of the names of the sequences and the dictionary of the sequences
    
    Parameters
    -------------------
        filename: path to the fasta file
    
    Returns
    -------------------
        order: list of the names of the sequences (str)
        sequences: dictionary of the sequences (str) with the name of the sequence (str) as a key
    '''
    try:
        f = open(filename)
    except IOError:                     
        print ("The file, %s, does not exist" % filename)
        return
    order = []
    sequences = {}

    for line in f:
        if line.startswith('>'):
            name = line[1:].rstrip('\n')
            #name = name.replace('_', ' ')
            order.append(name)
            sequences[name] = ''
        else:
                sequences[name] += line.rstrip('\n')

    print ("%d sequences found" % len(order))
    return order, sequences

def seq_check(row):
    '''
    Extract sequences from the contigs using start and stop locations, compare the results with the sequences in the annotation file and return whether they match or not
    
    Parameters
    -------------------
        row: pd.Series, whose keys are ['gene_id'(str),'contig_id'(str),'strand'(str),'start'(np.int64),'stop'(np.int64),'function'(str),'comment'(str),'sequence_nt'(str),'sequence_aa'(str),'sequence'(str)]
    
    Returns
    -------------------
        pd.Series, whose keys are ['nt_identity'(np.bool_),'aa_identity'(np.bool_),'res_seq_nt'(str),'res_seq_aa'(str)]
    '''
    nt_identity = 0
    start = row.start
    stop = row.stop
    strand = row.strand
    seq = row.sequence
    seq_nt = row.sequence_nt
    seq_aa = row.sequence_aa

    
    # error handling for wrong strand
    if strand == '+' and (stop < start):
        return False
    elif strand == '-' and (start < stop):
        return False
    
    if strand == '+':
        first_idx = start - 1
        last_idx = stop # not need -1
        res_seq_nt = seq[first_idx:last_idx]
        res_seq_aa = str(Seq(res_seq_nt, generic_dna).translate(table="Bacterial", to_stop=True))
        if res_seq_nt == seq_nt:
            nt_identity = True
        else:
            nt_identity = False
        if res_seq_aa == seq_aa:
            aa_identity = True
        else:
            aa_identity = False
            
    elif strand == '-':
        first_idx = stop - 1
        last_idx = start # not need -1
        res_seq_nt = str(Seq(seq[first_idx:last_idx]).reverse_complement())
        res_seq_aa = str(Seq(res_seq_nt, generic_dna).translate(table="Bacterial", to_stop=True))
        if res_seq_nt == seq_nt:
            nt_identity = True
        else:
            nt_identity = False
        if res_seq_aa == seq_aa:
            aa_identity = True
        else:
            aa_identity = False
#             print(nt_identity)
#             print('seq_ext',res_seq)
#             print('seq_nt',seq_nt)
    else:
        return 'unknown strand value'

    return pd.Series([nt_identity,aa_identity,res_seq_nt,res_seq_aa])

def codon_check(row):
    '''
    Compare the start and stop codons of the genes and the ones in the start_codon_list and stop_codon_list and return whether they match or not
    Find any stop codons in the sequences of the genes and return the location when there are
    
    Parameters
    -------------------
        row: pd.Series, whose keys are ['gene_id'(str),'contig_id'(str),'strand'(str),'start'(np.int64),'stop'(np.int64),'function'(str),'comment'(str),'sequence_nt'(str),'sequence_aa'(str),'sequence'(str),'nt_identity'(np.bool_),'aa_identity'(np.bool_),'res_seq_nt'(str),'res_seq_aa'(str)]
    
    Returns
    -------------------
        pd.Series, whose keys are ['start_codon_identity'(np.bool_),'stop_codon_identity'(np.bool_),'data'(list)]
    '''
    seq = row.sequence_nt
    data = []
    start_codon_list = ['ATG','GTG','TTG','CTG','ATT']
    stop_codon_list = ['TAG','TGA','TAA'] 
    
    if seq[:3] in start_codon_list:
        start_codon_identity = True
    else:
        start_codon_identity = False
    if seq[-3:] in stop_codon_list:
        stop_codon_identity = True
    else:
        stop_codon_identity = False

    codon_check = seq
    i = 0
    while (len(codon_check) > 3):
        if codon_check[:3] in stop_codon_list:
            position = i*3+1
            data.append(position)
        i += 1
        codon_check = codon_check[3:]  

    return pd.Series([start_codon_identity,stop_codon_identity,data])    
