import sys
import argparse
import numpy as np
import time
import zipfile
from collections import defaultdict
import pickle
from multiprocessing import Pool
import os.path
from os import path
start_time = time.time()

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count > 3000000:
                    break
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
def generate_kmers(k, read):
    return [read[i:i+k] for i in range(0, len(read), k)]

def kmers_dict(k, input_reads):
    dkmers = {}
    for read in input_reads:
        dkmers[read] = generate_kmers(k, read)
    return dkmers

def index_genome(k, genome):
    indices = defaultdict(list)
    for i in range(len(genome)-k):
        sub = genome[i:i+k]
        indices[sub].append(i)
    return indices

def calc_diff(seq1, seq2):
    return sum(1 for s1, s2 in zip(seq1, seq2) if s1 != s2)

def correct_snps(snps):
    print('correcting!')
    freq = {}
    for snp in snps:
        if tuple(snp) not in freq:
            freq[tuple(snp)] = 1
        else:
            freq[tuple(snp)] += 1
    corrected_snps = []
    for key, value in freq.items():
        if value > 26:
            corrected_snps.append(list(key))
    return corrected_snps

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

def is_noncontig_substr(needle, haystack, ref_index):
    i = 0
    l = len(needle)
    inserted = ""
    first_index = None
    for j in range(len(haystack)):
        if i < l and haystack[j] == needle[i]:
                i += 1
        else:
            if first_index == None:
                first_index = ref_index + j
            inserted += haystack[j]
    if i == l:
            return [inserted, first_index]
    return None

def find_snps(input_reads, reference, indices, dkmers):
    print('finding snps!')
    snps = {}
    corrected_snps= []
    for read in input_reads:
        kmers = dkmers[read]
        for kmer in kmers:
            if kmer in indices:
                index = indices[kmer]
                j = read.find(kmer)
                for ind in index:
                    if(ind-j < 0 or ind-j+len(read) > len(reference)):
                        continue
                    difference = calc_diff(reference[ind-j:ind-j+len(read)], read)
                    if difference < 4:
                        for i in range(len(read)):
                            if read[i] != reference[ind+i-j]:
                                if (reference[ind+i-j], read[i], ind+i-j) not in snps:
                                    snps[(reference[ind+i-j], read[i], ind+i-j)] = 1
                                else:
                                    snps[(reference[ind+i-j], read[i], ind+i-j)] += 1
                                    if snps[(reference[ind+i-j], read[i], ind+i-j)] == 26:
                                        corrected_snps.append([reference[ind+i-j], read[i], ind+i-j])
    return corrected_snps
    
def find_indels(k, input_reads, reference, indices, dkmers):
    print('finding indels!')
    ins = []
    dels = []
    for read in input_reads:
        kmers = dkmers[read]
        if kmers[0] and kmers[2] in indices:
            index1 = indices[kmers[0]]
            index3 = indices[kmers[2]]
            for ind1 in index1:
                if ind1+32 not in index3:
                    for ind3 in index3:
                        if 2*k < ind3-ind1 < 2*k+4:
                            deleted_candidate = is_noncontig_substr(kmers[1], reference[ind1+k:ind3], ind1+k)
                            if deleted_candidate != None and deleted_candidate not in dels:
                                dels.append(deleted_candidate)
                        elif 0 < ind3-ind1 < 2*k:
                            inserted_candidate = is_noncontig_substr(reference[ind1+k:ind3], kmers[1], ind1+k)
                            if inserted_candidate != None and inserted_candidate not in ins:
                                ins.append(inserted_candidate)
        if kmers[0] and kmers[1] in indices:
            index1 = indices[kmers[0]]
            index2 = indices[kmers[1]]
            for ind1 in index1:
                if ind1+16 not in index2:
                        for ind2 in index2:
                            if k < ind2-ind1 < k+4:
                                deleted_candidate = is_noncontig_substr("", reference[ind1+k:ind2], ind1+k)
                                if deleted_candidate != None and deleted_candidate not in dels:
                                    dels.append(deleted_candidate)
        if kmers[1] and kmers[2] in indices:
            index2 = indices[kmers[1]]
            index3 = indices[kmers[2]]
            for ind2 in index2:
                if ind2+16 not in index3:
                        for ind3 in index3:
                            if k < ind3-ind2 < k+4:
                                deleted_candidate = is_noncontig_substr("", reference[ind2+k:ind3], ind2+k)
                                if deleted_candidate != None and deleted_candidate not in dels:
                                    dels.append(deleted_candidate)
    return ins, dels

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
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
    reads = []
    for pairs in input_reads:
        for read in pairs:
            reads.append(read)

    dkmers = kmers_dict(10, reads)
    indices = index_genome(10, reference)

    snps = find_snps(reads, reference, indices, dkmers)

    dkmers = kmers_dict(16, reads)
    indices = index_genome(16, reference)

    indels = find_indels(16, reads, reference, indices, dkmers)
    insertions = indels[0]
    deletions = indels[1]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

print("--- %s seconds ---" % (time.time() - start_time))