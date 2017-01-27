#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

####################################################################
##This script takes an alignment (fasta)  and outputs the most   
##likely lineage of M.tb.
####################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_arguments(): 
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Lineage Assignment")
    input_file = parser.add_mutually_exclusive_group(required=True)
    input_file.add_argument('-i', '--inputFile',
        help ='alignment in fasta format', 
        type=is_file)
    return parser.parse_args()

args = get_arguments()
alnIN = args.inputFile

def make_dict():
    d = {}
    """for every diagnostic snp, record alt allele and lineage"""
    with open("/home/oneill/scripts/MtbLineageClassifier/Coll2014_LinSpeSNPs_final.csv", 'r') as infile:
        for i, line in enumerate(infile):
            line = line.strip().split("\t")
            if i > 0:
                lin = line[0]
                pos = int(line[1])
                allele = line[3].split("/")[1]
                if pos in d:
                    print("Position {0} found more than once, deleting".format(pos))
                    del d[pos]
                else:
                    d[pos] = [allele, lin]
    return d

def assign_lin(alnIN, d):
    assignments = {}
    for seq_record in SeqIO.parse(alnIN, "fasta"):
        assignments[seq_record.id] = []
        print("Processing sample {0}".format(seq_record.id))
        for i in range(len(seq_record.seq)):
            pos = i + 1
            if pos in d:
                if seq_record.seq[i] == d[pos][0]:
                    assignments[seq_record.id].append(d[pos][1])
    return assignments

def analyze_assignments(assignments):
    for sample in assignments:
        print("\n" + sample)
        Lins = set(assignments[sample])
        if len(Lins) > 0:
            for i in Lins:
                v = assignments[sample].count(i)
                print("{0} snps supporting {1}".format(str(v), i))
            lin = max(assignments[sample], key = assignments[sample].count)
            print("Most probable assignment for {0} is {1}".format(sample, lin))
        else:
            print("No diagnostic snps identified in {0}".format(sample))


d = make_dict() 
assignments = assign_lin(alnIN, d)
analyze_assignments(assignments)
