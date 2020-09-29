"""
Description: Class of genbank data
Author: Koen van den Berg
"""

# Imports
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import sys

class GenbankData:
    """Describes the Genbank parse manner of data and corresponding data
    structures.
    """

    # Initialization Method
    def __init__(self):
        """Initialization of Genbank Class datastructures
        """
        self.record = {}
        self.record_list = []

    # Other Methods
    def parseData(self, infile):
        """Parses the genbank record into self.record_list
        --------------------------------------------------
        infile:
            string, name of the infile
        --------------------------------------------------
        returns:
            None
        """
        i = 0
        gbkcontents = SeqIO.parse(infile, 'genbank')
        for record in gbkcontents:
            self.record['organism'] = record.annotations['organism']
            all_locs = [feature.location for feature in record.features if feature.type == 'CDS']
            for feature in record.features:
                if feature.type == 'CDS':
                    # Parse straightforward items:
                    f_rec = self._collectFeature(feature)
                    # Parse the coding sequence (cds):
                    self.record['dna'] = str(self.record['location'].extract(record.seq))

                    # Parse the 5'UTR sequence:
                    fiveprime_UTR = self._collectUpstream(all_locs, i)
                    self.record['fiveprime_UTR'] = fiveprime_UTR.extract(record.seq)
                    i += 1

                    self.record_list.append(self.record)

    def _collectFeature(self, feature):
        """Parses the easy features into self.record
        --------------------------------------------------
        feature:
            Biopython object, containing qualifiers
        --------------------------------------------------
        returns:
            None
        """
        d = feature.qualifiers
        self.record['locus_tag'] = d['locus_tag'][0] if 'locus_tag' in d else ''
        self.record['old_locus_tag'] = d['old_locus_tag'][0] if 'old_locus_tag' in d else ''
        self.record['gene'] = d['gene'][0] if 'gene' in d else ''
        self.record['product'] = d['product'][0] if 'product' in d else ''
        self.record['protein_id'] = d['protein_id'][0] if 'protein_id' in d else ''
        self.record['translation'] = d['translation'][0] if 'translation' in d else ''
        loc = feature.location
        strand = feature.location.strand
        self.record['location'] = loc
        self.record['strand'] = strand
                    
    def _collectUpstream(self, mem, idx):
        """finds the upstream region for a location for both plus and minus
        strands
        --------------------------------------------------
        mem:
            list, containing all feature locations from the genbank data
        idx:
            int, index of feature CDS
        --------------------------------------------------
        returns:
            ret = FeatureLocation object of upstream region
        """
        n = 30 # max size upstream region
        # PLUS STRANDS
        if self.record['strand'] == 1:
            for i in range(1,len(mem)):
                loc = mem[idx-i] # search backward for previous_end
                if loc.strand == 1:
                    previous_end = loc.end
                    break
            current_start = self.record['location'].start
            igr = current_start - previous_end
            if igr > n: # Take only the upstream region upto n
                ret = FeatureLocation(current_start-n, current_start, strand=1)
            elif igr < 0: # If genes are too close dont return region
                ret = FeatureLocation(0,0)
            else: # Return intergenetic region between 0 and n
                ret = FeatureLocation(current_start-igr, current_start, strand=1)
        # MINUS STRANDS
        if self.record['strand'] == -1:
            try: # search forward for next start
                for i in range(1,len(mem)):
                    loc = mem[idx+i]
                    if loc.strand == -1:
                        next_start = loc.start
                        break
            except(IndexError): # this is for the last location
                next_start = 1e20 # set to way too large number will be
                                  # clipped of by n anyways
            current_end = self.record['location'].end
            igr = next_start - current_end
            if igr > n:
                ret = FeatureLocation(current_end, current_end+n, strand=-1)
            elif igr < 0:
                ret = FeatureLocation(0,0)
            else:
                ret = FeatureLocation(current_end, current_end+igr, strand=-1)
        return(ret)

    
    def addNewDatatoRecord(self, description, func, *args, **kwargs):
        """Fills the record list with new columns i.e. each record gains a
        new entry
        --------------------------------------------------
        record_list:
            list, contains the records in dict format
        description:
            string, will show in terminal and as column header 
            in end result
        func:
            function, function that computes the parameter
        *args:
            arguments for the function
        **kwargs:
            yield/return: setting for a yielding/returning function
            loading_off: setting for turning of progress bar
        --------------------------------------------------
        returns:
            description = string
        """
        if kwargs['loading'] == True:
            sys.stdout.write(f'{description:50} :: ...working...\r')
        # Fill record list:
        for i, record in enumerate(self.record_list, start=1):
            # Functions that yield:
            if kwargs['setting'] == 'yield':
                for p1, p2 in func(record, *args):
                    record[p1] = p2
            # Functions that return:
            elif kwargs['setting'] == 'return':
                record[description] = func(record, *args)
        if kwargs['loading'] == True:
            sys.stdout.flush()
            sys.stdout.write(f'{description:50} :: DONE            \n')
        return(description)


        
