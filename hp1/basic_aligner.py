import sys
import argparse
import time
import zipfile


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need

"""

def generate_kmers(k, text):
    kmers = []
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])
    return kmers

def index_genome(genome):
	k = 10
	indices = {}
	for i in range(len(genome)-k+1):
		kmer = genome[i:i+k]
		if kmer in indices:
			indices[kmer].append(i)
		else:
			indices[kmer] = [i]
	return indices

def lookup(kmers, indices):
    possible = {}
    for kmer in kmers:
        if kmer in indices:
            possible[kmer] = indices[kmer]
    return possible

def calc_diff(seq1, seq2):
    return sum(1 for s1, s2 in zip(seq1, seq2) if s1 != s2)

def find_snps(input_reads, reference):
    indices = index_genome(reference)
    snps = []
    for pairs in input_reads:
        kmers = []
        for read in pairs:
            kmers = generate_kmers(10, read)
            lookups = lookup(kmers, indices)
            for kmer in kmers:
                if kmer in lookups:
                    index = lookups[kmer]
                    j = read.find(kmer)
                    for ind in index:
                        difference = calc_diff(reference[ind-j:ind-j+len(read)], read)
                        if difference < 3:
                            for i in range(len(read)-1):
                                if read[i] != reference[ind+i-j]:
                                    snps.append([reference[ind+i-j], read[i], ind+i-j])
    return snps

def correct_snps(snps):
    freq = {}
    for snp in snps:
        if tuple(snp) not in freq:
            freq[tuple(snp)] = 1
        else:
            freq[tuple(snp)] += 1
    corrected_snps = []
    for key, value in freq.items():
        if value > 100:
            corrected_snps.append(list(key))
    return corrected_snps
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    """
        TODO: Call functions to do the actual read alignment here
        
    """
    
    snps = correct_snps(find_snps(input_reads, reference))

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
