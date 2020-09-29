"""
Description: Functions for mutating exising parameters derived from Genbank file 
to seq (CDS) features for machine learning downstream.
Author: Koen van den Berg
"""

# Imports
import statistics

# Definitions
def p_collectNc(record):
    '''Computes the Nc for a given CDS sequence
    Forumla Faa: Faa = (n*P-1)/(n-1)
                 P = sum(p_i^2)
                 p_i = n_i/n
    Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
    Fn is the average homozygosity (Faa) for the aa 
    having a degeneracy of n. 
    --------------------------------------------------
    record: 
        dictionary, contains data for each rowid, uses:
        - cds -- string of dna/rna
    --------------------------------------------------
    returns:
        float, Nc value
    '''
    CDS = record['dna']
    codon_dict = {
          "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0, "CTT": 0,
          "CTC": 0, "CTA": 0, "CTG": 0, "ATT": 0, "ATC": 0,
          "ATA": 0, "ATG": 0, "GTT": 0, "GTC": 0, "GTA": 0,
          "GTG": 0, "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
          "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0, "AAT": 0,
          "AAC": 0, "AAA": 0, "AAG": 0, "GAT": 0, "GAC": 0,
          "GAA": 0, "GAG": 0, "TCT": 0, "TCC": 0, "TCA": 0,
          "TCG": 0, "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
          "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0, "GCT": 0,
          "GCC": 0, "GCA": 0, "GCG": 0, "TGT": 0, "TGC": 0,
          "TGA": 0, "TGG": 0, "CGT": 0, "CGC": 0, "CGA": 0,
          "CGG": 0, "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
          "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}

    # Counting occurences of codons and filling codon_dict
    codon_list = [CDS[i:i+3] for i in range(0, len(CDS), 3)]
    for codon in codon_list[:-1]: # skip last codon
        if codon in codon_dict:
            codon_dict[codon] += 1
    synonymous_codon_dict = {
           "CYS": ["TGT", "TGC"],
           "ASP": ["GAT", "GAC"],
           "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
           "GLN": ["CAA", "CAG"],
           "MET": ["ATG"],
           "ASN": ["AAC", "AAT"],
           "PRO": ["CCT", "CCG", "CCA", "CCC"],
           "LYS": ["AAG", "AAA"],
           #"STP": ["TAG", "TGA", "TAA"],
           "THR": ["ACC", "ACA", "ACG", "ACT"],
           "PHE": ["TTT", "TTC"],
           "ALA": ["GCA", "GCC", "GCG", "GCT"],
           "GLY": ["GGT", "GGG", "GGA", "GGC"],
           "ILE": ["ATC", "ATA", "ATT"],
           "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
           "HIS": ["CAT", "CAC"],
           "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
           "TRP": ["TGG"],
           "VAL": ["GTA", "GTC", "GTG", "GTT"],
           "GLU": ["GAG", "GAA"],
           "TYR": ["TAT", "TAC"]}

    # Computing Faa's and adding them to degeneracy groups
    F_dict = {'F1':[], 'F2':[], 'F3':[], 'F4':[], 'F6':[]}
    for aa in synonymous_codon_dict:
        P = 0
        n = 0
        local_codon_list = synonymous_codon_dict[aa]
        for codon in local_codon_list:
            n += codon_dict[codon]
        try:
            for codon in local_codon_list:
                ni = codon_dict[codon]
                pi = ni/n
                P += pi**2
            Faa = ((n * P)-1)/(n-1)
        except(ZeroDivisionError):
            Faa = 0
        Fn = len(synonymous_codon_dict[aa]) # degeneracy group n
        if Faa > 0:
            F_dict[f'F{Fn}'].append(Faa)

    # If no codons are found for a certain group (e.g. F4), then let
    # it be 1. This will thus not impact the Nc in the forumla
    # specified below. F3 will be the average of F2 and F4:
    F_dict = {k:v if len(v)>0 and k!='F3' else [1] for k,v in F_dict.items()}
    if len(F_dict['F3']) == 0:
        F3 = (statistics.mean(F_dict['F2']) + statistics.mean(F_dict['F4']) )/ 2
        F_dict['F3'].append(F3)
    Nc = 2 + 9/statistics.mean(F_dict['F2']) + \
         1/statistics.mean(F_dict['F3']) + \
         5/statistics.mean(F_dict['F4']) + \
         3/statistics.mean(F_dict['F6'])

    # Nc cannot be higher than 61 biologically speaking
    Nc = Nc if Nc < 61 else 61
    return(Nc)
