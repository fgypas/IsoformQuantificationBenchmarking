#!/usr/bin/python

from __future__ import with_statement # Workaround to support Python 2.5

from sets import Set
import gzip
import sys
import os
import optparse
import multiprocessing


#import affinity

IF1_NAME = 0
IF1_CHR = 1
IF1_STRAND = 2
IF1_START = 3
IF1_STOP = 4

EXON_STRAND = 0
EXON_START = 1
EXON_STOP = 2
EXON_RANGE = 3
EXON_CHR = 4

MATCH_START = 1
MATCH_STOP = 2
MATCH_NAME = 3
MATCH_SCORE = 4

input_if1_path = ""
input_matches = ""
input_map_path = ""
complete_coordinates = None
max_processors = -1;

parser = optparse.OptionParser("usage: %prog -i if1_file -d map_file")
parser.add_option("-i", "--if1-input", dest="if1", type="string", help="Path to the IF1 file corresponding to the matches bed file)", default=None)
parser.add_option("-c", "--complete_coordinates", dest="complete_coordinates", type="int", help="If True, the full list of start and end coordinates for each traversed exons will be returned. Otherwise, only the start and end coordinates on the transcript will be  returned.", default=0)       
parser.add_option("-d", "--dictionary", dest="dictionary", type="string", help="Path the the dictonary (hash) that maps the exon names to the matches bed file's numeric values. If this parameter is not provided the same path as the if1 file shall be assumed with the extension replaces to '.map'.", default="None") # $ is used as an invalid token for a path.
parser.add_option("-p", "--processors", dest="processors", type="int", help="Number of processors to use", default=1)

(options, args) = parser.parse_args(sys.argv)
if None in options.__dict__.values():
    parser.error("Incorrect number of arguments. Did you specify the if1 (-i) and matches (-m) for the input files?")
    sys.exit(-1)
else:
    input_if1_path = os.path.abspath(options.if1)
    input_map_path = os.path.abspath(options.dictionary)
    complete_coordinate = options.complete_coordinates
    max_processors = options.processors

#affinity.set_process_affinity_mask(0,2**max_processors-1)

# Returns a tuple for each lines from file at filepath, splitted on 'delimiter' token (\t for bed)
# This generator is ment to return a one time use iterator. It is important to note that
# an instance of a generator can only yield each tuple once. However, this avoids the need
# to store content in memory as the generator iterates until it is consumed and discarded.
def readfile(filepath, delimiter):
    if filepath[-2:] == "gz":
        g = gzip.open(filepath, 'r')
        for line in g:
            yield tuple(line.split(delimiter))

    else:            
        with open(filepath, 'r') as f:
            for line in f:
                yield tuple(line.split(delimiter))

# Call's the generator and appended yielded content to a file. Wil break after max_lines argument is provided.
transcript_map = {}
exon_dict = {}
def get_coordinates():
    # Create map
    map_as_tuples = readfile(input_map_path, "\t")

    for muple in map_as_tuples:
        mkey = muple[0]
        mvalue = muple[1]
        if mvalue[0] == '>':
            mvalue = mvalue[1:] # If sequence starts by > then pop it out of the string
        transcript_map[mkey] = mvalue

    lines_as_tuples = readfile(input_if1_path, "\t")

    for data in lines_as_tuples:
        exon = data[IF1_NAME]
        exon_name = exon.split("_")[0] # ie. ENST00000237247.6_0 exon_name = ENST00000237247.6
        exon_part = exon.split("_")[1] #                         exon_part = 0
        start = int(data[IF1_START])
        stop = int(data[IF1_STOP].strip())
        strand = data[IF1_STRAND]
        chr = data[IF1_CHR]
        exon_dict[(exon_part, exon_name)] = (strand, start, stop, str(int(stop)-int(start)-1), chr)



    # After ths we have a exon_dict fully populated as above. We need now to parse the matches.
    count = 0
    matches_as_lines  = []
    pool = multiprocessing.Pool(max_processors)
    chunky_bacon = 100000
    for line in sys.stdin.readlines():
        matches_as_lines.append(line)
        if count == chunky_bacon:
            res = pool.map(process_line, matches_as_lines, max_processors)            
            for r in res:
                print r
            res = None
            count = 0
            matches_as_lines = []
        else:
            count +=1
    if count > 0:
        res = pool.map(process_line, matches_as_lines, max_processors)            
        for r in res:
            print r
    
    #matches_as_tuples = readfile(input_matches_path, "\t")
    
    ## That's all folks


def process_line(data_line):
    #Split on \t
    data = data_line.split("\t")

    if data[6].strip() == "-":
        return
    remainder = []
    if len(data) > 6:
        for i in range(6,len(data)):
            remainder.append(data[i])
       
    transcript_id = data[0]

    transcript_name = transcript_map[transcript_id]

    transcript_parts = transcript_name.split("|")

    exon_name = transcript_parts.pop(0)

    exon_name = exon_name.split("_")[0]

    start_point = int(data[MATCH_START])
    stop_point = int(data[MATCH_STOP])
    test_no = data[3].split("|")[0]

    match_range = stop_point - start_point
    match_name = data[MATCH_NAME]
    match_score = data[MATCH_SCORE]

    transcript_coordinates = {}
    transcript_pointer = 0

    flag = False
    strand = exon_dict[("0", exon_name)][EXON_STRAND]
    nucleotides_per_exons = []


    for t in transcript_parts:
        index = t.strip() # Pesky \n ...
        current_exon = exon_dict[(index, exon_name)]

        range_list = range(current_exon[EXON_START], current_exon[EXON_STOP]+1)
        virtual_start = 0
        virtual_stop = 0
        current_exon_nuc_count = 0
        if strand == "+":
            virtual_start = start_point
            virtual_stop = stop_point -1
        elif strand ==  "-":
            virtual_start = start_point
            virtual_stop = stop_point -1
            #virtual_start = stop_point - 1
            #virtual_stop = start_point
            range_list.reverse()
        for i in range_list:
            if transcript_pointer >= virtual_start and transcript_pointer <= virtual_stop:
                   current_exon_nuc_count += 1

            transcript_coordinates[transcript_pointer] = i
            transcript_pointer += 1
        if current_exon_nuc_count > 0:
            nucleotides_per_exons.append(str(current_exon_nuc_count))

    nuc_per_exons = "|".join(nucleotides_per_exons)
    if strand == "-":
        text_to_append = exon_dict[("0", exon_name)][EXON_CHR],"\t",str(transcript_coordinates[stop_point-1]-1),"\t",str(transcript_coordinates[start_point]),"\t",match_name,"\t",match_score,"\t",strand
    elif strand == "+":
        text_to_append = exon_dict[("0", exon_name)][EXON_CHR],"\t",str(transcript_coordinates[start_point]-1),"\t",str(transcript_coordinates[stop_point-1]),"\t",match_name,"\t",match_score,"\t",strand

    text_to_append = "".join(text_to_append)
    remainder_string = ""
    if len(remainder) > 0:
        remainder_string += "\t" + "\t".join(remainder)
        remainder_string = remainder_string.rstrip()

    text_to_append += remainder_string + "\t" + nuc_per_exons
    return text_to_append
                            
#EOPositionFromBed

# Execution
get_coordinates()
