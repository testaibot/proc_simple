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


known_event2 = {
                ('M','S'):2, 
                ('S','M'):1
               }

def is_event2(event2):
    if event2 in known_event2:
        return True
    return False

known_event3 = {
                ('M','N','M'):2, 
                ('M','I','M'):2,
                ('M','D','M'):2
               }

def get_pos_offset(increment, operation):
    pos = 0
    offset = 0

    if operation == 'M':
        pos += increment
        offset += increment
    elif operation == 'I':
        offset += increment
    elif operation == 'D':
        pos += increment
    elif operation == 'N':
        pos += increment
    elif operation == 'S':
        offset += increment

    return(pos, offset)

def po_inc(old, new):
    old = (old[0] + new[0], old[1] + new[1])
    return(old)

def event_offsets(old_offset, operations_range, increments_range):  
    offset_list = [old_offset]
    for i in range(len(operations_range)):
        old_offset = po_inc(old_offset, get_pos_offset(increments_range[i], operations_range[i]))
        offset_list.append(old_offset)

    return offset_list


def is_event3(event3):
    if event3 in known_event3:
        return True
    return False

def process_cigar_operations(cigar_increments, cigar_operations):
    offset = (0, 0)
    i = 0
    inf_loop_breaker = 0
    while i < len(cigar_operations)-1:
        if i < len(cigar_operations)-1:
            operations_range = cigar_operations[i:i+2]
            increments_range = cigar_increments[i:i+2]
            if(is_event2(operations_range)):
                ev_offsets = event_offsets(offset, operations_range, increments_range)
                yield (operations_range, ev_offsets)
                offset = ev_offsets[len(ev_offsets)-1]
                i+=known_event2[operations_range]
        
        if i < len(cigar_operations)-2:
            operations_range = cigar_operations[i:i+3]
            increments_range = cigar_increments[i:i+3]
            if(is_event3(operations_range)):
                ev_offsets = event_offsets(offset, operations_range, increments_range)
                yield (operations_range, ev_offsets)
                offset = ev_offsets[len(ev_offsets)-1]
                i+=known_event3[operations_range]
        
        if inf_loop_breaker > 100:
            print(cigar_operations)
            raise Exception("Unknown cigar_template")
        else: inf_loop_breaker+=1

def bulk_1(line):
    spl = {} 

    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()
       
    spl["cell_id"] = temp[0].split(":")[0]
    spl["flag"] = int(temp[1])
    spl["ref"] = temp[2]
    spl["pos"] = int(temp[3])
    spl["cigar"] = temp[5]
    if spl["cigar"] == '*':
        return 
    spl["match"] = temp[9]
    return spl 

def smart_seq2(line):
    spl = {}   

    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()
       
    spl["cell_id"] = temp[0].split(".")[0]
    spl["flag"] = int(temp[1])
    spl["ref"] = temp[2]
    spl["pos"] = int(temp[3])
    spl["cigar"] = temp[5]
    if spl["cigar"] == '*':
        return 
    spl["match"] = temp[9]
    return spl

def process_sam_line(line):

    spl = smart_seq2(line)
    if spl == None: return

    cell_id = spl["cell_id"]
    flag = spl["flag"]
    ref  = spl["ref"]
    pos = spl["pos"]
    cigar = spl["cigar"]
    match  = spl["match"]

    strand = (flag & BAM_FREVERSE)>0
    offset = 0
    match_start = 0
    match_intervals = []
    event_list = []

    cigar_split = regex_split(cigar)

    operations = tuple(cigar_split[i] for i in range(1, len(cigar_split), 2))
    increments = tuple(int(cigar_split[i]) for i in range(0, len(cigar_split), 2))

    if len(operations) == 1: return

    
    for event, offsets_range in process_cigar_operations(increments, operations):
                
        if event == ('M', 'N', 'M'):
            if(offsets_range[1][1] - offsets_range[0][1] < 15): continue
            if(offsets_range[3][1] - offsets_range[2][1] < 15): continue

            if(not stat_only): 
                print(
                    ref, 
                    pos + offsets_range[1][0] - 1, 
                    pos + offsets_range[2][0], 
                    STRAND[strand], 
                    joinStrings(",", event),
                    cell_id, 
                    match[offsets_range[0][1]:offsets_range[1][1]], 
                    match[offsets_range[2][1]:offsets_range[3][1]],
                    sep = "\t"
                 )
        if event == ('M', 'S'):
            if(offsets_range[1][1] - offsets_range[0][1] < 15): continue
            if(offsets_range[2][1] - offsets_range[1][1] < 15): continue

            ev_stats_send(match[offsets_range[0][1]:offsets_range[1][1]])
            ev_stats_send(match[offsets_range[1][1]:offsets_range[2][1]])
            if(not stat_only): 
                print(
                    ref,
                    pos + offsets_range[1][0] - 1, 
                    pos + offsets_range[1][0] - 1, 
                    STRAND[strand], 
                    joinStrings(",", event),
                    cell_id, 
                    match[offsets_range[0][1]:offsets_range[1][1]], 
                    match[offsets_range[1][1]:offsets_range[2][1]],
                    sep = "\t"
                 )

        if event == ('S', 'M'):
            if(offsets_range[1][1] - offsets_range[0][1] < 15): continue
            if(offsets_range[2][1] - offsets_range[1][1] < 15): continue

            ev_stats_send(match[offsets_range[0][1]:offsets_range[1][1]])
            ev_stats_send(match[offsets_range[1][1]:offsets_range[2][1]])
            
            if(not stat_only):
                print(
                    ref,
                    pos + offsets_range[1][0] - 1, 
                    pos + offsets_range[1][0] - 1, 
                    STRAND[strand], 
                    joinStrings(",", event),
                    cell_id, 
                    match[offsets_range[0][1]:offsets_range[1][1]], 
                    match[offsets_range[1][1]:offsets_range[2][1]],
                    sep = "\t"
                 )

ev_stats = {}
def ev_stats_send(item):
    if item not in ev_stats:
        ev_stats[item] = 0
    ev_stats[item]+=1

parser = argparse.ArgumentParser(description='Hola\n')
parser.add_argument('-f', "--file", dest="file", metavar='file', type=str, default = None,
                    help='sam-file')
parser.add_argument('-nb', "--nbins", dest="nbins", metavar='nbins', type=int,
                    help='number of bins', default=2147483647)

parser.add_argument('-s', "--stat", dest="stat_only", metavar='stat_only', type=int,
                    help='print contamination stat only', default=False)

args = parser.parse_args()

file = args.file
nbins = args.nbins
stat_only = args.stat_only
BAM_FREVERSE = 0x10
STRAND = ["+", "-"]

if(file==None):
    f = fileinput.input()
else:
    f = open(file).readlines()


if stat_only:
    from tqdm import tqdm
    f = tqdm(f)

for line in f:
    process_sam_line(line)

if stat_only:
    for ev in sorted(ev_stats, key=ev_stats.get):
        print("ev_stat:", ev, ev_stats[ev])
