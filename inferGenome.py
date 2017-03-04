import argparse
import numpy as np

# Function for parsing the arguments
def parse_args():

	parser = argparse.ArgumentParser(description='input parameters')
	parser.add_argument('--sham', dest='sham', help='Input sham file')
	parser.add_argument('--out', dest='out', help='output tsv file')
	parser.add_argument('--p01', dest='p01', help='0 to 1 error rate')
	parser.add_argument('--p10', dest='p10', help='1 to 0 error rate')
	parser.add_argument('--prior1', dest='prior1', help='prior probability that each base is a 1')

	args = parser.parse_args()
	return args
	
# Function for reading sham file
def read_sham_file(sham):
    with open (sham, 'r') as shamfile:
        data = shamfile.readlines()
    locations = []
    reads = []
    for d in data:
        rd = d.split('\t')
        locations.append(int(rd[0]))
        reads.append(rd[1].rstrip())
    return locations, reads
	
# Function for writing the result
def write_result(out, posterior_prob, max_len):
	with open (out, 'w') as out_file:
		for i in range(max_len):
			out_file.write("{}\t{}\n".format(i, format(posterior_prob[i], '.3f')))	
	
def main():

	# Parsing the arguments
	args = parse_args()
	sham = args.sham
	out = args.out
	p01 = float(args.p01)
	p10 = float(args.p10)
	prior1 = float(args.prior1)
	
	# Reading the input file
	locations, reads = read_sham_file(sham)
	
	# Calculating the max size of the predicted string to pre-alocate the arrays
	max_len = max([locations[i]+len(reads[i]) for i in range(len(reads))])
	
	# Number of zeros at every location
	zeros = np.zeros(max_len, int)
	# Number of ones at every location
	ones = np.zeros(max_len, int)
	
	# Loop over reads to calculate the number of ones and zeros at every location
	for i in range(len(reads)):
		rd_ones = np.array([int(n) for n in reads[i]])
		lc = locations[i] 
		ones[lc:lc+rd_ones.shape[0]] += rd_ones 
		zeros[lc:lc+rd_ones.shape[0]] += (1-rd_ones)
		
	# Pre-alocate the posterior probabilities
	posterior_prob = np.zeros(max_len, float)
	
	# Calculate the posterior probabilities at every location
	for i in range(max_len):
		numerator = (p10**zeros[i]) * ((1-p10)**ones[i]) * prior1
		denominator = numerator + (p01**ones[i]) * ((1-p01)**zeros[i]) * (1-prior1)
		posterior_prob[i] = numerator / denominator
		
	# Writing the result to the output file
	write_result(out, posterior_prob, max_len)
			
if __name__ == "__main__":
    main()		