## Alexander Kanitz, Biozentrum, University of Basel
## 30-SEP-2014
## Based on code provided by Simon Anders (http://seqanswers.com/forums/showthread.php?t=12070)

# Usage: downsample_fastq.py <input file: FASTQ> <fraction: FLOAT> <output file: FASTQ>

import sys, random, itertools, HTSeq

in_file = iter( HTSeq.FastaReader( sys.argv[1] ) )
fraction = float( sys.argv[2] )

out_file = open( sys.argv[3], "w" )

for read in in_file:
   if random.random() < fraction:
      read.write_to_fasta_file( out_file )
      
out_file.close()
