#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 21-MAY-2015                                       ###
#########################################################

import argparse

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument('--in', required = True, dest = 'inp_file', help = 'Insert the output of RSEQ tool.')
	parser.add_argument('--out', required = True, dest= 'out_file', help = 'Output file.')
	args = parser.parse_args()

	f_out = open(args.out_file, 'w')
	with open(args.inp_file) as f_in:
		for line in f_in:
			cols = line.strip().split('\t')
			tr = cols[3]
			tr_v = cols[4]
			
			if(',' in tr):
				tr_sp = tr.split(',')
				tr_v_sp = tr_v.split(',')
				if(len(tr_sp)==len(tr_v_sp)):
					for i in range(0,len(tr_sp)):
						f_out.write(tr_sp[i]+"\t"+tr_v_sp[i]+"\n")
				else:
					sys.stderr.write("Error exists!!!\n")
					sys.exit(-1)
			else:
				f_out.write(tr+"\t"+tr_v+"\n")
	f_out.close()			

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt!\n")
		sys.exit(-1)	
