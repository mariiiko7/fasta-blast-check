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

def print_all_hsps(alignment):
    '''
    print all data in result.alignment
    '''
    hsps = alignment.hsps
    index = 0
    hsp = hsps[index]
    print('title>>>' + str(alignment))
    print('bits>>>' + str(hsp.bits))
    print('expect>>>' + str(hsp.expect))
    print('frame>>>' + str(hsp.frame))
    print('gaps>>>' + str(hsp.gaps))
    print('identities>>>' + str(hsp.identities))
    print('match>>>' + str(hsp.match))
    print('num_alignments>>>' + str(hsp.num_alignments))
    print('positives>>>' + str(hsp.positives))
    print('query>>>' + str(hsp.query))
    print('query_start>>>' + str(hsp.query_start))
    print('sbjct>>>' + str(hsp.sbjct))
    print('sbjct_start>>>' + str(hsp.sbjct_start))
    print('score>>>' + str(hsp.score))
    print('strand>>>' + str(hsp.strand))
    

def get_top_alignment(result):
    '''
    Extract top alignment from result.alignments and return their rank, identities, start_difference, query_start, sbjct_start, start_difference, alignment_length and alignment
    
    Parameters
    -------------------
        result: NCBIStandalone.BlastParser.parse(file) (refer: http://www.dalkescientific.com/writings/NBN/blast_parsing.html)
    
    Returns
    -------------------
        target_data_list: a list of data, whose keys are ['rank'(int), 'identities'(float), 'query_start'(int), 'sbjct_start'(int),'start_difference'(int), 'alignment_length'(int), alignment'(result.alignment)]
    '''
# alignments are sorted by rank in default
# it should have one or zero element. in order to apply the same function as n31 to this return, make return value a list
    target_data_list = []
    data = {}
    if len(result.alignments) > 0:
        alignment = result.alignments[0]
        data['identities'] = alignment.hsps[0].identities[0]/alignment.hsps[0].identities[1]
        data['query_start'] = alignment.hsps[0].query_start
        data['sbjct_start'] = alignment.hsps[0].sbjct_start
        data['start_difference'] = alignment.hsps[0].query_start - alignment.hsps[0].sbjct_start
        data['alignment_length'] = len (alignment.hsps[0].query)
        data['alignment'] = alignment
        target_data_list.append(data)
    return target_data_list

def make_data_frame_from_list(result_list):
    '''
    Extract information of top alignment and convert it to data frame which contains g_id, gene_name, identities, query_start, sbjct_start, start_difference and alignment_length
    
    Parameters
    -------------------
        result_list: NCBIStandalone.BlastParser.parse(file) (refer: http://www.dalkescientific.com/writings/NBN/blast_parsing.html)
    
    Returns
    -------------------
        df: dataframe, whose values are ['g_id'(str),  'identities'(float), 'query_start'(int), 'sbjct_start'(int),'start_difference'(int), 'alignment_length'(int)]
    '''
    final_outputs = []
    for result in result_list:
        data = {}
        data['g_id'] = result['g_id']
        data['gene_name'], data['identities'], data['query_start'], data['sbjct_start'], data['start_difference'], data['alignment_length'], = wrap_get_attributes(result['outputs'])
        final_outputs.append(data)
    df = pd.DataFrame(final_outputs)
    df = df.sort_values('g_id')    
    return df
    
def wrap_get_attributes(outputs):
    '''
    Extract information of top alignment which contains gene_name, identities, query_start, sbjct_start, start_difference and alignment_length
    
    Parameters
    -------------------
        outputs: NCBIStandalone.BlastParser.parse(file) (refer: http://www.dalkescientific.com/writings/NBN/blast_parsing.html)
    Returns
    -------------------
        'gene_name'(str), 'identities'(float), 'query_start'(int), 'sbjct_start'(int),'start_difference'(int), 'alignment_length'(int)
    '''
    if len(outputs) < 1:
        gene_name = np.nan
        identities = np.nan 
        query_start = np.nan
        sbjct_start = np.nan
        start_difference = np.nan
        alignment_length = np.nan
    else:
        for output in outputs:
            gene_name = output['alignment'].title
            identities = output['identities']
            query_start = output['query_start']
            sbjct_start = output['sbjct_start']
            start_difference = output['start_difference']
            alignment_length = output['alignment_length']
    return gene_name, identities, query_start, sbjct_start, start_difference, alignment_length

