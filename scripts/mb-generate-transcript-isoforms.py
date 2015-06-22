#!/usr/bin/python

from __future__ import with_statement #Workaround to support generators in Python 2.5
from collections import defaultdict

import os
import sys
import optparse
import gzip

class FastaImploder:
    fasta_path = ""
    output_fasta_path = ""
    output_map_path = ""
    skips = 0
    max_sequences = -1
    debug_return = None
    def __init__(self):
        parser = optparse.OptionParser("usage: %prog -e exon_file -o fasta_output -m mapping_output -n number_of_skips(optional)")
        parser.add_option("-e", "--exonfasta", dest="exonfasta", type="string", help="Path fasta file containing the exons", default=None)
        parser.add_option("-o", "--outputfasta", dest="outputfasta", type="string", help="Path to file that will contain the fasta isoform output.", default=None)
        parser.add_option("-m", "--outputmap", dest="outputmap", type="string", help="Path to file that will contain the mapping output.", default=None)
        parser.add_option("-n", "--num-skipping-events", dest="numberOfSkippingEvents", default=0, type="int", help="Number of skips to generate")


        (options, args) = parser.parse_args(sys.argv)
        if None in options.__dict__.values():
            parser.error("Incorrect number of arguments. Did you specify input and output files?")
            sys.exit(-1)
        else:
            self.fasta_path = options.exonfasta
            self.output_fasta_path = options.outputfasta
            self.output_map_path = options.outputmap
            self.skips = options.numberOfSkippingEvents

        #Create a file to store in the permutations

        self.output_fasta_path = os.path.abspath(self.output_fasta_path)
        self.output_map_path = os.path.abspath(self.output_map_path)
        open(self.output_fasta_path, 'w').close() #NOTE: This will clear any content already in the file
        open(self.output_map_path, 'w').close() #NOTE: This will clear any content already in the file

                
    # Returns a tuple for each lines from file at filepath, splitted on 'delimiter' token (\t by defaults)
    # This generator is ment to return a one time use iterator. It is important to note that
    # an instance of a generator can only yield each tuple once. However, this avoids the need
    # to store content in memory as the generator iterates until it is consumed and discarded.
    def readfile(self, filepath, delimiter):
        if filepath[-2:] == "gz":
            g = gzip.open(filepath, 'r')
            for line in g:
                yield tuple(line.split(delimited))
        else:
            with open(filepath, 'r') as f:
                for line in f:
                    yield tuple(line.split(delimiter))

    # Similar to readfile, but doesnt break it down in a tuple (not required as we want the full string anyways)
    def read_fasta(self, filepath):
        if filepath[-2:] == "gz":
            g = gzip.open(filepath, 'r')
            for line in g:
                yield line
        else:
            with open(filepath, 'r') as f:
                for line in f:
                    yield line

    # Help method to get the line numbers in a file
    def file_len(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
            return i + 1
        
    # Create two generators to map Fasta to Bed. Each file has its own generator.
    def implode_fasta(self):
         fasta_lines = {}
         fasta = self.read_fasta(self.fasta_path)

         root_exon = ""
         root_exon_list = []
         root_exon_strand = {}
         counter = 0
         map = {}
         for header in fasta:
               # The fasta header is implicitely removed by header, and the actually sequence is returned by  as fasta.next() is called
               exon_name = header.split("|")[0]
               exon_name = exon_name[1:] # Get rid of that pesky '>'
               part_number = exon_name.split("_")[-1] # Get the part number

               if part_number == "0":
                   root_exon = exon_name
                   root_exon_list.append(root_exon)
                   root_exon_strand[root_exon] = header.split("|")[2]
                   fasta_lines[root_exon] = []

               fasta_sequence = fasta.next()
               fasta_sequence = fasta_sequence.strip().upper() #Trim whitespaces at the begining and tail and append to next line
               if root_exon == '':
                   print fasta_sequence
                   print part_number
                   print header
                   raw_input()
               
               fasta_lines[root_exon].append(fasta_sequence)

         sequence_counter = 0
         with open(self.output_fasta_path, 'a') as f:
             for root in root_exon_list:
                 # Implode sequences
                 exploded_transcript= []
                 exon_numbers = []
                 exon_number = 0
                 for fasta_line in fasta_lines[root]:
                     exploded_transcript.append(fasta_line)
                     exon_numbers.append(exon_number)
                     exon_number = exon_number + 1
                 if exon_number < 2:
                     continue

                 # Handle negative strand
                 if root_exon_strand[root] == "-": #Negative strand
                     exploded_transcript.reverse()
                     exon_numbers.reverse()
                 elif root_exon_strand[root] == "+": #Positive strand
                     pass # Do nothing
                 else: #This should not happen...
                     print "COULD NOT GET STRAND; Expecting '-' or '+', got",root_exon_strand[root]
                     sys.exit(-1)


                 parts_string = '|'.join(str(v) for v in exon_numbers)
                 imploded_transcript = "".join(exploded_transcript)

                 root = root[:-1] # Remove the root exon's number after the _ 

                 map[sequence_counter] = ">"+ root + "|"  + parts_string
                 f.write(">"+str(sequence_counter)+"\n")
                 f.write(imploded_transcript)
                 f.write("\n")

                 sequence_counter += 1
                
                 #Insert the set of ordered strings with 1 deletion
                 if len(exon_numbers) > 2  and self.skips >= 1:               
                     for i in range(len(exon_numbers)-1):
                         if i == 0:
                             continue
                         parts_list = exon_numbers[:]
                         parts_list.pop(i)
                         parts_list = '|'.join(str(x) for x in parts_list)
                         
                         modified_exploded_transcript = []
                         modified_exploded_transcript = exploded_transcript[:]
                         modified_exploded_transcript.pop(i)
                         
                         imploded_transcript = "".join(modified_exploded_transcript)
                         map[sequence_counter] = ">"+ root + "|"  + parts_list
                         
                         f.write(">"+str(sequence_counter)+"\n")
                         f.write(imploded_transcript)
                         f.write("\n")
                         
                         sequence_counter += 1

                 #Insert the set of ordered strings with 2 deletion
                 if len(exon_numbers) > 3 and self.skips >= 2:
                     for i in range(len(exon_numbers)-2):
                         if i == 0:
                             continue
                         parts_list = exon_numbers
                         parts_list.pop(i)
                         parts_list.pop(i)
                         parts_list = '|'.join(str(x) for x in parts_list)
                         
                         modified_exploded_transcript = []
                         modified_exploded_transcript = exploded_transcript[:]
                         modified_exploded_transcript.pop(i)
                         modified_exploded_transcript.pop(i)
                         
                         
                         imploded_transcript = "".join(modified_exploded_transcript)
                         map[sequence_counter] = ">"+ root + "|"  + parts_list
                         
                         f.write(">"+str(sequence_counter)+"\n")
                         f.write(imploded_transcript)
                         f.write("\n")
                         
                         sequence_counter += 1
                         
         f.close()

         with open(self.output_map_path, 'a') as file:
            for i in range(len(map)):
                if map[i][0] == ">":
                    real_name = map[i][1:];
                else:
                    real_name = map[i]
                file.write(str(i)+"\t"+real_name)
                file.write("\n")

#EOFastaImploder

# Main Execution
exploder = FastaImploder()
#exploder.validate_fasta_vs_bed()
exploder.implode_fasta()
