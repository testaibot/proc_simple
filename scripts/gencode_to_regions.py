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

def joinStrings(separator, stringList):
    return separator.join(str(string) for string in stringList)


def process_annotation_line(line):
    if "transcript_id" not in line: return
    
    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()
    
    print(temp[0], temp[3], temp[4], temp[6], sep = '\t', end = '\t')
    temp9 = temp[8].split(';')
    for it in temp9:
        it = it.strip()
        it = it.split(' ')
        if(it[0] in ["gene_name", "transcript_name"]):
            it[1] = it[1].strip('\"')
            print('\t', it[1], sep = '', end = '')
    print()

parser = argparse.ArgumentParser(description='Hola\n')
parser.add_argument('-f', "--file", dest="file", metavar='file', type=str, default = None, help='annotation file')
args = parser.parse_args()

file = args.file

if(file==None):
    f = fileinput.input()
else:
    f = open(file).readlines()


for line in f:
    process_annotation_line(line)

