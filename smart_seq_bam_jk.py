#!/usr/bin/python
from __future__ import print_function
import sys
import math
import argparse
import ntpath
import os
import re
import fileinput
from itertools import zip_longest

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)

def joinStrings(stringList):
    return ''.join(str(string) for string in stringList)

parser = argparse.ArgumentParser(description='Hola\n')

parser.add_argument('-f', "--file", dest="file", metavar='file', type=str, default = None,
                    help='sam-file')
parser.add_argument('-l', "--limit", dest="lim", metavar='lim', type=int,
                    default=0, help='stop after this number of lines (for debug)')
parser.add_argument('-nb', "--nbins", dest="nbins", metavar='nbins', type=int,
                    help='number of bins', default=1)

args = parser.parse_args()


file = args.file
lim = args.lim
nbins = args.nbins
BAM_FREVERSE = 0x10
STRAND = ["+", "-"]

regex = re.compile("(\d+)(\w)")
n_counts = {}

n = 0

if(file==None):
    f = fileinput.input()
else:
    f = open(file).readlines()

for line in f:
    n += 1
    if(n>lim & lim>0):
        break
    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()
       
    cell_id = temp[0].split(".")[0]
    flag = int(temp[1])
    ref = temp[2]
    pos = int(temp[3])
    #qual = int(temp[4])
    cigar = temp[5]
    if cigar == '*':
        continue    
    #match = temp[9]

    strand = (flag & BAM_FREVERSE)>0
    offset = 0
    for increment, operation in grouper(list(filter(None, regex.split(cigar))), 2):
        increment = int(increment)
        if operation == 'M':
            pos += increment
            offset += increment
        elif operation == 'I':
            offset += increment
        elif operation == 'D':
            pos += increment
        elif operation == 'N':
           
            curr_key = joinStrings([ref, "_", pos - 1, "_", pos + increment, "_", STRAND[strand]], "\t", cell_id)

            if curr_key not in n_counts:
                n_counts[curr_key] = 1
            else:
                n_counts[curr_key] = n_counts[curr_key] + 1
                            
            pos += increment
        elif operation == 'S':
            offset += increment

for curr_key in n_counts:
    print(curr_jk, n_counts[curr_key], sep = "\t", end="\n")

