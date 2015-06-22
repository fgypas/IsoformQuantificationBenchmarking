#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 21-MAY-2015                                       ###
#########################################################

import sys

input = sys.argv[1]
output = sys.argv[2]

f = open(output,'w')
with open(input) as fp:
    for line in fp:
        if(line[0]==">"):
            splited = line.strip().split()
            if(len(splited)>=3):
                gene = splited[3].replace("gene:","")
                transcript = splited[0].replace(">","")
                f.write(">"+gene+"$$"+transcript+"\n")
        else:
            f.write(line)
f.close()
