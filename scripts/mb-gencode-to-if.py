#!/usr/bin/python

from __future__ import with_statement #Workaround to support Python 2.5 (on bc2)'

"""Script which takes in a bed file and converts it to Internal Format number 1

   Output: NAME CHROM STRAND START END

"""

import sys
import gzip

class InternalFormat1Maker:
    input_path = ""
    output_path = ""
    max_lines = -1
    def __init__(self):
        if len(sys.argv) < 2:
            print "Usage: %prog <pathToBedFile>"
            sys.exit()
        if len(sys.argv) > 1:
            self.input_path = sys.argv[1]

# Returns a tuple for each lineos from file at filepath, splitted on 'delimiter' token (\t by defaults)
# This generator is ment to return a one time use iterator. It is important to note that
# an instance of a generator can only yield each tuple once. However, this avoids the need
# to store content in memory as the generator iterates until it is consumed and discarded.
    def readfile(self, filepath, delimiter):
        if filepath[-2:] == "gz":
            g = gzip.open(filepath, 'r')
            for line in g:
                yield tuple(line.split(delimiter))                          
        else:
            with open(filepath, 'r') as f:
                for line in f:
                    yield tuple(line.split(delimiter))

# Call's the generator and appended yielded content to a file. Wil break after max_lines argument is provided.
    def convertToIF1(self):
        lines_as_tuples = self.readfile(self.input_path, "\t")
        i = 0

        # Get column number for each required field
        header_tuple = lines_as_tuples.next()

        counter = 0
        NAME_FIELD = None
        CHROM_FIELD = None
        STRAND_FIELD = None
        EXONSTARTS_FIELD = None
        EXONENDS_FIELD = None
        
        for field in header_tuple:
            field_part = field.split(".")
            field_part = field_part[-1] # Get last element
            if field_part == "name" :
                NAME_FIELD = counter
            elif field_part == "chrom" :
                CHROM_FIELD = counter
            elif field_part == "strand" :
                STRAND_FIELD = counter
            elif field_part == "exonStarts" :
                EXONSTARTS_FIELD = counter
            elif field_part == "exonEnds" :
                EXONENDS_FIELD = counter
            counter += 1

        # Start needs to be incremented to be inclusive.
        list_of_tuples = []
        for data in lines_as_tuples:
            name = data[NAME_FIELD]
            chrom = data[CHROM_FIELD]
            strand = data[STRAND_FIELD]
            exonStarts = data[EXONSTARTS_FIELD].split(",")
            exonEnds = data[EXONENDS_FIELD].split(",")

            for i in range(len(exonStarts)):
                if exonStarts[i] == '':
                    continue #This is done because the end of exonsStart and exonsEnd has a comma ( , ) at the end, and this will make add an empty value to the array
                a_tuple = (name+"_"+str(i),chrom,strand,str(int(exonStarts[i])+1),exonEnds[i])
                list_of_tuples.append(a_tuple)
                        

        for tuple in list_of_tuples:
            newLine = "\t".join(tuple)
            print newLine
                


#EOInternalFormat1Maker

# Execution

fileMaker = InternalFormat1Maker()
fileMaker.convertToIF1()
