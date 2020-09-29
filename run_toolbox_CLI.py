#!/usr/bin/env python3

'''
MODIFIED VERSION Translational-Toolbox
--------------------------------------------------
Author: Koen van den Berg
University: VU & WUR
Department of Systems Bioinformatics
Date: Jan - Jul 2020
--------------------------------------------------
NOTE: This work is only suitable for review purposes. Big chunks of
the full version of the software are left out and this program will
thus not run.

'''

# Imports
from genbank.genbank import GenbankData
from collect_features import load_data, compute_features
import argparse

# Definitions
def get_arguments_cl():
    '''parses the command line arguments to a tuple
    --------------------------------------------------
    None
    --------------------------------------------------
    returns: 
        tuple with arguments
    '''
    # Main parser
    parser = argparse.ArgumentParser(description="Collects relevant mRNA\
    sequence parameters from annotated whole genome genbank files")
    subparsers = parser.add_subparsers(dest='mode') #title='subcommands', dest='mode')

    # collect_features
    parameter_parser = subparsers.add_parser('collect_features',
                                             help='Collect transcript parameters that are\
                                             relevant for protein modelling.',
                                             add_help=False)
    required = parameter_parser.add_argument_group('Required Arguments')
    optional = parameter_parser.add_argument_group('Optional Arguments')
    required.add_argument("-G", "--genome_file",
                            help='bacterial genome genbank file that is fully annotated. Examples are\
                            e.coli.gff, l.lactis.gff, etc... Genbank v<> or higher',
                            required=True)
    required.add_argument("-R", "--RNA_file",
                          help='Absolute path to comma seperated file (.csv) or Excel spreadsheet\
                          (.xlsx) containing a first column with a row ID similar to:\
                          {locus_tag, old_locus_tag, protein_id}. The remaining columns contain\
                          the mrna counts/fractions for all biological replicates in the\
                          condition. If NO mRNA data is available, specify "None" here and use\
                          the "-m" flag.',
                          required=True)
    required.add_argument("-P", "--protein_file",
                          help='Absolute path to comma seperated file (.csv) or Excel spreadsheet\
                          (.xlsx) containing a first column with a row ID similar to:\
                          {locus_tag, old_locus_tag, protein_id}. The remaining columns contain\
                          the protein counts/fractions for all biological replicates in the\
                          condition. If NO protein data available, specify "None" here and use\
                          the "-p" flag.',
                          required=True)
    required.add_argument("-B", "--biological_replicates",
                          help='NONE',
                          required=True)
    required.add_argument("-T", "--tAI_reference",
                          help='NONE',
                          required=True)
    required.add_argument("-O", "--output_dir",
                          help='NONE',
                          required=True)
    optional.add_argument("-h", "--help",
                          action="help",
                          help='show this help message and exit')
    optional.add_argument("-sc", "--scan_tRNA",
                          help='NONE',
                          required = False,
                          action = 'store_true')
    optional.add_argument("-p", "--no_protein_data",
                          help='Use this if there is NO protein data available. Specify "None" at\
                          the -P flag. Fills the output with a DUMMY variable that only\
                          contains NAs. Does not work together with the -r2 or -c options. \
                          Will work with expr as alias: the predict_proteins module, but \
                          no models can be created using a result with dummy protein values \
                          (default = false)',
                          required = False,
                          action = 'store_true')
    optional.add_argument("-m", "--no_mRNA_data",
                          help='Use this if there is NO mRNA data available. Specify "None" at the -R\
                          flag. Fills the output with a DUMMY variable that only contains\
                          NAs. Does not work together with the -r2 or -c options. Result cannot\
                          be applied to modelling/predicting (default = false)',
                          required = False,
                          action = 'store_true')
    optional.add_argument("-t", "--transform_data",
                          help='NONE',\
                          required = False,
                          action = 'store_true')
    optional.add_argument("-pl", "--plasmids",
                          help='NONE',
                          nargs='+',
                          required = False)
    optional.add_argument("-r2", "--R2_upper_limit",
                          help='NONE',
                          required = False,
                          action = 'store_true')
    optional.add_argument("-c", "--continue_models",
                          help='NONE',
                          required = False,
                          action = 'store_true')
    optional.add_argument("-d", "--developer_mode",
                          help='NONE',
                          required = False,
                          action = 'store_true')
    return(parser.parse_args())

##################################################
# MAIN collect parameters
##################################################
def main_collect_features(args):
    """Collects the model parameters for organism
    --------------------------------------------------
    args:
        namespace object that contains user arguments
    --------------------------------------------------
    """
    ##############################
    # PARSING GENBANK FILES
    ##############################
    # Here the list that holds the parameters is initialized
    CDS_record_list = GenbankData()
    CDS_record_list.parseData(args.genome_file)

    ##############################
    # LOADING RNA AND PROTEIN
    ##############################
    # Loading biological replicate group information
    #replicate_groups = load_data.dataToDict(args.biological_replicates)
    #replicate_groups = {k:list(v.values())[0] for k,v in replicate_groups.items()}
    
    # Loading RNA-seq data
    if args.no_mRNA_data:
        # Fill with NAs
        CDS_record_list.addNewDatatoRecord('NA mRNA values',
                       load_data.collectDummy, 'mrna', setting='yield', loading=True)
    else:
        # Normal procedure
        mrna_reference = load_data.dataToDict(args.RNA_file)
        CDS_record_list.addNewDatatoRecord('mrna samples', load_data.collectData,
                       mrna_reference, replicate_groups, 'mrna',
                       setting='yield', loading=True)

    # Loading Proteomics data
    if args.no_protein_data:
        # Fill with NAs
        CDS_record_list.addNewDatatoRecord('NA protein values',
                       load_data.collectDummy, 'prot', setting='yield', loading=True)
    else:
        # Normal procedure
        protein_reference = load_data.dataToDict(args.protein_file)
        CDS_record_list.addNewDatatoRecord('protein samples',
                       load_data.collectData, protein_reference,
                       replicate_groups, 'prot', setting='yield',
                       loading=True)

    ##############################
    # SEQUENCE PARAMETERS
    ##############################
    # Computing and adding Nc
    Nc = CDS_record_list.addNewDatatoRecord('Nc',
                                            compute_features.p_collectNc,
                                            setting='return',
                                            loading=True)

    print(CDS_record_list.record_list[:4])

    # Performing more mutations to the existing JSON data format in
    # the full version, but due to confidentiality that cannot be
    # disclosed.

if __name__ == '__main__':
    args = get_arguments_cl()

    if args.mode == 'collect_features':
        main_collect_features(args)

