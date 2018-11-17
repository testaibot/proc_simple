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

regex = re.compile("(\d+)(\w)")
def regex_split(cigar):
    return list(filter(None, regex.split(cigar)))
    

def process_sam_line(line):
    temp = line.split('\t')
    
    cigar = temp[5]
    if cigar == '*':
        return   

    temp[-1] = temp[-1].rstrip()
       
    cell_id = temp[0].split(".")[0]
    flag = int(temp[1])
    ref = temp[2]
    pos = int(temp[3])
 
    match = temp[9]

    strand = (flag & BAM_FREVERSE)>0
    offset = 0
    match_start = 0
    match_intervals = []
    event_list = []

    cigar_split = regex_split(cigar)
    for increment, operation in grouper(cigar_split, 2):
        increment = int(increment)
        if operation == 'M':
            pos += increment
            offset += increment
        elif operation == 'I':
            offset += increment
        elif operation == 'D':
            pos += increment
        elif operation == 'N':
            bin = nbins-1 if offset>=nbins else offset

            event = [ref, pos - 1, pos + increment, STRAND[strand], cell_id]
            event_list.append(event)

            match_intervals.append(match[match_start:offset])
            match_start = offset

            pos += increment
        elif operation == 'S':
            offset += increment

    match_intervals.append(match[match_start:offset])

    for i in range(len(event_list)):
        event_list[i].append(match_intervals[i])
        event_list[i].append(match_intervals[i+1])
        event = joinStrings("\t", event_list[i])
        yield event

parser = argparse.ArgumentParser(description='Hola\n')

parser.add_argument('-f', "--file", dest="file", metavar='file', type=str, default = None,
                    help='sam-file')
parser.add_argument('-nb', "--nbins", dest="nbins", metavar='nbins', type=int,
                    help='number of bins', default=2147483647)

args = parser.parse_args()


file = args.file
nbins = args.nbins
BAM_FREVERSE = 0x10
STRAND = ["+", "-"]

if(file==None):
    f = fileinput.input()
else:
    f = open(file).readlines()

for line in f:
    for event in process_sam_line(line):
        print(event)
