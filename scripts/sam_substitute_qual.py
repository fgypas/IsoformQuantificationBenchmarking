## Foivos Gypas, Alexander Kanitz
## 08-02-2015

## Substitutes QUAL field entries of a SAM file with the values of a pre-computed cPickle 'SEQ -> QUAL' dictionary

## Imports
import sys
from Bio.Seq import Seq
import cPickle as pickle

## Inputs
dict_path = sys.argv[1] # e.g. /import/bc2/home/zavolan/GROUP/ISO_BENCH/resources/seq.qual.dictionary.hsa1.cPickle

## Functions

def rev_complement(inp):
    seq =Seq(inp)
    return str(seq.reverse_complement())

## Main script

## Load dictionary
with open(dict_path,'rb') as fp:
    pickle_dict = pickle.load(fp)


# Variables initialization
prev_name = 'tpt'
removed_quality = ''

for line in sys.stdin:
    
    if(line[0] == '@'):
        sys.stdout.write(line)
    else:
    
        sp_line = str(line).strip().split("\t")
        current_name = sp_line[0]
        current_strand = int(sp_line[1])
        current_sequence = sp_line[9]
        
        if(prev_name == 'tpt'): # initialize
            prev_name = current_name
            prev_strand = current_strand
            if(current_strand == 0 or current_strand == 256):
                removed_quality = pickle_dict[current_sequence].pop(0)
                sp_line[10] = str(removed_quality)
            else:
                removed_quality = pickle_dict[rev_complement(current_sequence)].pop(0)[::-1]
                sp_line[10] = str(removed_quality)
        
        elif(prev_name != current_name and prev_name != 'tpt'):
            
            if(current_strand == 0 or current_strand == 256):
                removed_quality = pickle_dict[current_sequence].pop(0) # remove first quality score from list
                sp_line[10] = str(removed_quality)
            else:
                removed_quality = pickle_dict[rev_complement(current_sequence)].pop(0)[::-1]
                sp_line[10] = str(removed_quality)
                
        elif(prev_name == current_name and prev_name != 'tpt'):
            
            if(current_strand == 0 or current_strand == 256):
                if(prev_strand == 0 or prev_strand == 256):
                    sp_line[10] = str(removed_quality)
                else:
                    removed_quality = removed_quality[::-1]
                    sp_line[10] = str(removed_quality)
            else:
                if(prev_strand == 0 or prev_strand == 256):
                    removed_quality = removed_quality[::-1]
                    sp_line[10] = str(removed_quality)
                else:
                    sp_line[10] = removed_quality
        
        sys.stdout.write("\t".join(sp_line)+"\n")
            
        prev_name = current_name
        prev_strand = current_strand
        prev_sequence = current_sequence
