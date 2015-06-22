##################################################
## Authors: Foivos Gypas, Alexander Kanitz
## Affiliation: Biozentrum, University of Basel
## Created: 02-JUN-2014
## Last modified: 15-DEC-2014
##################################################

import sys
import argparse
import numpy as np


def main():

	parser = argparse.ArgumentParser("Calculate mean and std from SAM file (supplied as STDIN).")

	parser.add_argument('--mean', required = True, dest = 'mean_out', help = 'Text file stating the mean read length')
        parser.add_argument('--sd', required = True, dest = 'sd_out', help = 'Text file stating the standard deviation of the read lengths')
	parser.add_argument('--multimappers', action='store_true', help = 'Count multimapped reads once')

	args = parser.parse_args()

	if(args.multimappers):
		seq_lengths, seqs = get_seq_lengths_multimappers()
	else:
		seq_lengths, seqs = get_seq_lengths()
	mean, sd = calc_mean_sd(np.array(seq_lengths))

	mean_file = open(args.mean_out, 'w')
	mean_file.write(str(int(round(mean))))
	mean_file.close()

	sd_file = open(args.sd_out, 'w')
	sd_file.write(str(int(round(sd))))
	sd_file.close()

def calc_mean_sd(seq_lengths):

	return (np.mean(seq_lengths), np.std(seq_lengths))

def get_seq_lengths_multimappers():

	print("Calulating read lengths (mulimappers mode ON)")

	seqs = 0
	seq_lengths = []
	dict = {}

	for line in sys.stdin:
		if line[:1] == "@":
			continue
		split_line =line.strip().split('\t')
		read_name = split_line[0]
		seq = split_line[9]
		if (read_name in dict):
			continue
		else:	
			dict[read_name] = read_name
			if valid_seq(seq):
				seqs += 1
				seq_lengths.append(len(seq))
			else:	
				print("Invalid sequence detected: ")
				print(seq)
	return (seq_lengths, seqs)
				
	

def get_seq_lengths():
	
	print("Calulating read lengths (mulimappers mode OFF)")

	seqs = 0
	seq_lengths = []

	for line in sys.stdin:
		if line[:1] == "@":
			continue
		split_line =line.strip().split('\t')
		seq = split_line[9]
		if valid_seq(seq):
			seqs += 1
			seq_lengths.append(len(seq))
		else:	
			print("Invalid sequence detected: ")
			print(seq)
	
	return (seq_lengths, seqs)

def valid_seq(seq):
	for letter in seq.lower():
		if letter not in "atcgn":
			return(False)
	return(True)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(-1)
