#!/usr/bin/python

from __future__ import with_statement # Workaround to support Python 2.5

import sys
import os
import optparse
import re
import itertools
import multiprocessing
import time

max_processors = multiprocessing.cpu_count()
#def __init__(self):

parser = optparse.OptionParser("usage: %prog -p [processors, optional]")
parser.add_option("-p", "--processors", dest="processors", type="int", help="Number of processors to use", default=1)

(options, args) = parser.parse_args(sys.argv)
if None in options.__dict__.values():
    parser.error("Incorrect number of arguments. Did you specify the input sam file?")
    sys.exit(-1)
else:
    if options.processors >= 1:
        max_processors = options.processors       

    

# Returns a tuple for each lines from file at filepath, splitted on 'delimiter' token (\t for bed)
# This generator is ment to return a one time use iterator. It is important to note that
# an instance of a generator can only yield each tuple once. However, this avoids the need
# to store content in memory as the generator iterates until it is consumed and discarded.
def readfile(filepath, delimiter):
    with open(filepath, 'r') as f:
        for line in f:
            yield tuple(line.split(delimiter))

# Similar to readfile (see above) but instead it returns the entire line and not a tuple
def readlines(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            yield line               


def align():
    ts = time.time()
    ### Multiprocessing ### 
    pool = multiprocessing.Pool(max_processors)
    ### THIS WAS MULTIPROCESSING ###

    #sam_tuples = readfile(input_sam, "\t")
    count = 0
    data_lines = []
    batch_size = 100000
    for sam_line in sys.stdin.readlines():
        if sam_line[0] == "@":
            continue
        sam_tuple = tuple(sam_line.split("\t"))
        data_lines.append(sam_tuple)
        if count == batch_size:
            res = pool.map(generate_bed_line, data_lines, max_processors)
            for r in res:
                print r
            res = None
            count = 0
            data_lines = []
        else:
            count += 1
    if count > 0:
        res = pool.map(generate_bed_line, data_lines, max_processors)
        for r in res:
            print r
            
#    for firstsam in sam_tuples:
#        sam_list = []
#        sam_list.append(firstsam)


def f(x):
    return x

def generate_bed_line(sam):

    if sam[:1] == "@":
        return
    
    # CONSTANTS
    BED_CHR=0
    BED_START=1
    BED_STOP=2
    BED_READ=3
    BED_SCORE=4
    BED_STRAND=5
    BED_MID=6

    SAM_NUMBER=0
    SAM_STRAND=1
    SAM_CHR=2
    SAM_START=3
    SAM_SCORE = 4
    SAM_MID = 5
    SAM_STAR = 6 # What's this??
    SAM__ = 7
    SAM__ = 8
    SAM_READ=9
    SAM_STAR = 10
    SAM__ = 42
    SAM_INDEX_SPACE = 11
    SAM_MISMATCHES = 12

    SAM_POSITIVE = "0"
    SAM_NEGATIVE = "16" #what the hell were they thinking? Who wrote that spec?

    #Now reading in the sam_tuples sequentially.
    #Need to convert this in some kind of bed format

    MID_length = get_match_length(sam[SAM_MID])
    strand = "?"

    if sam[SAM_STRAND] == SAM_POSITIVE:
        strand = "+"
    elif sam[SAM_STRAND] == SAM_NEGATIVE:
        strand = "-"

    bed_line = ""
    bed_line += sam[SAM_CHR]+"\t"
    bed_line += str(int(sam[SAM_START])-1)+"\t"
    bed_line += str(int(sam[SAM_START])-1+MID_length)+"\t"
    bed_line += sam[SAM_NUMBER]+"\t"
    bed_line += sam[SAM_SCORE]+"\t"
    bed_line += strand
    #At this point, the bed file is almost complete, just needs the MMID
    #Make a tuple with it!

    bed = bed_line.split("\t")             


    g_length = MID_length

    g = "".ljust(g_length, "?")
    g = list(g)

    #print sam
    r = sam[SAM_READ]
    r = list(r)

    G = []
    R = []

    current_MID = sam[SAM_MID]
    exploded_MID = re.split(r'(\d+|[A-Za-z])', current_MID)
    for item in exploded_MID:
        if item == '':
            exploded_MID.remove(item)

    ### Do alignment read vs genome
    for i in range(len(exploded_MID)/2):
        number = int(exploded_MID[i*2])
        type = exploded_MID[i*2 +1]

        if type == "M":
            for _ in itertools.repeat(None, number):
                G.append(g.pop(0))
                R.append(r.pop(0))                        
        elif type == "I":
            for _ in itertools.repeat(None, number):
                G.append("-")
                R.append(r.pop(0))
        elif type == "D":
            for _ in itertools.repeat(None, number):
                G.append(g.pop(0))
                R.append("-")


    mismatch = sam[SAM_MISMATCHES]
    mismatch = mismatch.split(":")[-1]
    mismatch = mismatch.replace("^", "")

    mismatch_splitted = re.split(r'(\d+|[A-Za-z])', mismatch)
    for item in mismatch_splitted:
        if item == '':
            mismatch_splitted.remove(item)

    pointer = 0

    for c in mismatch_splitted:

        if c == '':
            continue # Skip c and this loop when this happens

        # Skip the next G[pointer] that is not a minus strand
        if c.isdigit():
            value = int(c)
            while value > 0:
                if G[pointer] == '-':
                    pointer += 1
                else:
                    pointer += 1
                    value -= 1
        else:
            while G[pointer] == '-':
                pointer += 1 
            G.pop(pointer)
            G.insert(pointer, c)
            pointer+=1


    match_count = 0
    MMID = ""


    g_index = 0
    r_index = 0

    for index in range(len(G)):
        cg = G[index]
        cr = R[index]
        if cg == "?":
            match_count += 1
        else:
            if match_count > 0:
                MMID += str(match_count)
                match_count = 0

            if cg == "-":
                MMID += "I"+cr

            elif cr == "-":
                MMID += "D"+cg

            else:
                MMID += "M"+cg+cr

        index += 1

    if match_count > 0:
        MMID += str(match_count)
        match_count = 0

    index_space = sam[SAM_INDEX_SPACE]
    index_space = index_space.split(":")[-1]

    new_bed = list(bed)
    new_bed.pop(4)
    new_bed.insert(4, index_space)
    new_bed.append(MMID)
    new_file_line = "\t".join(new_bed)
    return new_file_line

# Parses a string containing Matches (M) and InDels (I and D) and returns the number of positions represented by that string
def get_match_length( match_string):
    number_string = ""
    sequence_length = 0

    match_splitted = re.split(r'(\d+)', match_string)
    match_splitted.pop(0) # re produces a weird '' as elem 0  ... 

    for i in range(len(match_splitted)/2):
        value = match_splitted[i*2]
        type = match_splitted[i*2 + 1]

        if type == "M" or type == "D":
            sequence_length += int(value)
        else:
            pass

    return sequence_length

                
align()
