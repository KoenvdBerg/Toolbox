"""
Description: Functions for loading in the RNA-seq and proteomics data
Author: Koen van den Berg
"""

# Imports
from pathlib import Path
import pandas as pd
import sys
import math

# Definitions
def collectData(record, ref, replicates, dt):
    '''collects the data from the reference, protein 
    or mRNA
    --------------------------------------------------
    record: 
        dictionary, contains data for each rowid, uses:
        - locus_tag -- string
        - old_locus_tag -- string
        - prot_id -- string
    ref:
        dict, mrna/protein data from csv that resides in dict
    replicates:
        dict, contains group information per sample
    dt:
        string, type of data that is loaded in (mrna/prot)
    --------------------------------------------------
    returns:
        iterator object containing k=mrnakey, v=value
    '''
    locus = record['locus_tag']
    old_locus = record['old_locus_tag']
    prot_id = record['protein_id']
    if locus in ref:
        values = ref[locus]
    elif old_locus in ref:
        values = ref[old_locus]
    elif prot_id in ref:
        values = ref[prot_id]
    else:
        l_ref = list(ref.values())[1] # taking some value to extract column names
        ref_keys = l_ref.keys()
        values = {k:'NA' for k in ref_keys}
    for k, v in values.items():
        # extract replicate group:
        try:
            group = replicates[k]
        except(KeyError):
            sys.stdout.write(f'Error: data column encountered in input protein/mRNA file that is not\
            present in the biological_replicates file.\nColumn: {k}\nTo solve this\
            error, put the column in the biological_replicates file and mention to\
            which biological group it belongs, or remove the column from the input\
            data file.\n\n')
            exit()
        try:
            v = float(v)
            if math.isnan(v): # NAN values being stubborn and thus removed here
                v = 'NA'
            yield(f'{dt}={group}|{k}', v) # gensym replaced return(v)
        except(ValueError, TypeError):
            # ValueError: if v is NA
            # TypeError: if v is NoneType
            v = 'NA'
            yield(f'{dt}={group}|{k}', v) # gensym replaced return(v)


def collectDummy(record, dt):
    """Collects NAs for absent mrna or protein data
    --------------------------------------------------
    record:
        dictionary, contains data for each rowid, uses:
        - NONE
    dt:
        string, data type
    --------------------------------------------------
    yields:
        name, NA
    """
    yield(f'{dt}=DUMMY', 'NA')


def dataToDict(data_file, sheetname = None, output_format = 'index'):
    """parses the .csv or .xlsx to dictionary
    --------------------------------------------------
    data_file:
        string, absolute path to the .csv file or .xlsx file
    sheetname:
        string, name of the sheetname. None default
    output_format:
        string, name of pandas output dictionary
    --------------------------------------------------
    returns:
        ret = {rowkey1: {col1:v, col2:v, col3:v, ...}
               rowkey2: {...}}
    """
    file_type = Path(data_file).suffix
    if file_type == '.xlsx':
        if sheetname:
            data = pd.read_excel(data_file, sheet_name = sheetname, index_col=0)
        else:
            if output_format == 'records':
                data = pd.read_excel(data_file)
            else:
                data = pd.read_excel(data_file, index_col=0)
        ret = data.to_dict(output_format) # {row1: {col1: 0.1, col2: 0.6}, {row2:{...}}, ...}
    elif file_type == '.csv':
        if output_format == 'records':
            data = pd.read_csv(data_file, sep = ',')
        else:
            data = pd.read_csv(data_file, index_col=0, sep = ',')
        ret = data.to_dict(output_format) # {row1: {col1: 0.1, col2: 0.6}, {row2:{...}}, ...}
    return(ret)
