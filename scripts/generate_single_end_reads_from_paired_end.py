########################################
## Author: Foivos Gypas
## email: foivos.gypas@unibas.ch
## Date: 26 MAY 2015
########################################

from Bio import SeqIO
import sys

handle = open(sys.argv[1], "rU")
count = 1
for record in SeqIO.parse(handle, "fasta") :
	if(record.description[-3] is 'S'):
		print('>s'+str(count))
		print(record.seq.upper())
		count+=1
handle.close()
