from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


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


"""
    TODO: Use this space to implement any additional functions you might need

"""

def generate_kmers(k, text):
    kmers = []
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])
    return kmers

def correct_kmers(kmers):
    freq = {}
    for kmer in kmers:
        if kmer not in freq:
            freq[kmer] = 1
        else:
            freq[kmer] += 1
    corrected_kmers = []
    for key, value in freq.items():
        if value > 3:
            corrected_kmers.append(key)
    return corrected_kmers

def debruijn_graph(kmers):
    graph = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix not in graph:
            graph[prefix] = [suffix]
        else:
            graph[prefix].append(suffix)
    return graph

def degrees(graph):
    indegree = {}
    outdegree = {}
    for node, next_nodes in graph.items():
        outdegree[node] = len(next_nodes)
        if node not in indegree:  
            indegree[node] = 0
        for next_node in next_nodes:
            if next_node not in indegree:
                indegree[next_node] = 1
            else:
                indegree[next_node] += 1
            if next_node not in outdegree:
                outdegree[next_node] = 0
    return indegree, outdegree

def remove_tips(graph):
    outdegree = degrees(graph)[1]
    while 0 in outdegree.values():
        for onode in list(outdegree.keys()):
            if outdegree[onode] == 0:
                node_to_remove = onode
                for gnode, next_nodes in graph.items(): 
                    if node_to_remove in next_nodes:
                        for i in range(next_nodes.count(node_to_remove)):
                            next_nodes.remove(node_to_remove)
                        outdegree[gnode] -= 1
                outdegree.pop(onode)
    return graph 
                                   
def generate_contigs(graph, node, indegree, outdegree):
    contigs = []
    for next_node in graph[node]:
        new_path = [node, next_node]
        ins, outs = indegree[next_node], outdegree[next_node]
        while ins == 1 and outs == 1:
            node = next_node
            next_node = graph[node][0]
            new_path.append(next_node)
            ins, outs = indegree[next_node], outdegree[next_node]
        contigs.append(new_path)
    return contigs

def debruijn_to_paths(graph):
    paths = []
    indegree, outdegree = degrees(graph)
    for node in outdegree:
        ind, outd = indegree[node], outdegree[node]
        if outd > 0 and not (outd == 1 and ind == 1):
            paths.extend(generate_contigs(graph, node, indegree, outdegree))
    return paths

def assemble_contigs(paths):
    sequence = paths[0]
    for kmer in paths[1:]:
        sequence += kmer[-1]
    return sequence

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    """
            TODO: Call functions to do the actual assembly here

    """
    
    kmers = []
    k = 28
    for pair in input_reads:
        kmers += generate_kmers(k, pair[0][::-1])
        kmers += generate_kmers(k, pair[1])

    graph = remove_tips(debruijn_graph(correct_kmers(kmers)))
    paths = debruijn_to_paths(graph)
    contigs = map(assemble_contigs, paths)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
