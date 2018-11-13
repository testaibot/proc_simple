#!/usr/bin/python

from __future__ import print_function
import sys
import argparse
import ntpath
import os
import fileinput
from itertools import zip_longest


known_events = {
                #(1,0):"jk",
                (2,1,0):"A5", 
                (1,2,0):"A3", 
                (2,1,2,0):"cass",
                #(2,1,2,1,2,0):"mut_ex",
                #(3,2,3,2,3,0):"mut_ex_cass"
               }

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)

def joinStrings(separator, stringList):
    return separator.join(str(string) for string in stringList)


def pos_inc(chr, k, d, sign):
    if chr not in d:
        d[chr] = {}
        
    if k not in d[chr]:
        d[chr][k] = 0
        
    d[chr][k] += sign

def pos_append(cell_id, chrom, key, umi, d):
    if cell_id not in d:
        d[cell_id] = {}
        
    if chrom not in d[cell_id]:
        d[cell_id][chrom] = {}
        
    if key not in d[cell_id][chrom]:
        d[cell_id][chrom][key] = {}
        
    if umi not in d[cell_id][chrom][key]:
        d[cell_id][chrom][key][umi] = 0
        
    d[cell_id][chrom][key][umi] += 1
    
def event_append(chr, key, val, d):
    val = tuple(val)
    if chr not in d:
        d[chr] = {}
        
    if key not in d[chr]:
        d[chr][key] = {}
        
    if val not in d[chr][key]:
        d[chr][key][val] = 0
        
    d[chr][key][val] += 1
  

def classify_event(ev):
    ev = tuple(event)
    if ev in known_events:
        return known_events[ev]
    
    return joinStrings(",",ev)


def get_strand_counts(cell_id_list, chr, pos, l_count, umi = True):
    count_strand_cl = [0,0]
    count_strand = [0,0]
    for cell_id in cell_id_list:
        if chr in l_count[cell_id]:
            if (pos, "+") in l_count[cell_id][chr]:
                count_strand[0] = len(l_count[cell_id][chr][(pos, "+")]) if umi else sum(l_count[cell_id][chr][(pos, "+")].values())

            if (pos, "-") in l_count[cell_id][chr]:
                count_strand[1] = len(l_count[cell_id][chr][(pos, "-")]) if umi else sum(l_count[cell_id][chr][(pos, "-")].values())
            count_strand_cl[0] += count_strand[0]
            count_strand_cl[1] += count_strand[1]
    return count_strand

def get_ev_counts(cell_id_list, coord, l_tot, r_tot):
    ev_counts = []
    for pos in coord:
        l_counts = get_strand_counts(cell_id_list, chr, pos, l_tot)
        r_counts = get_strand_counts(cell_id_list, chr, pos, r_tot)
        l_strand = "=" if l_counts[0] == l_counts[1] else "+" if l_counts[0] > l_counts[1] else "-"
        r_strand = "=" if r_counts[0] == r_counts[1] else "+" if r_counts[0] > r_counts[1] else "-"
                    
        l_counts = sum(l_counts)
        r_counts = sum(r_counts)
        if r_counts == l_counts:
            return(None)
            #raise Exception('Error: expected simple linear element. Input: composite')
        elif r_counts > l_counts:
            counts = r_counts
            strand = r_strand
        else:
            counts = l_counts
            strand = l_strand
        ev_counts.append([strand, counts])
    return(ev_counts)

def get_ev_region(ev_coords):
    return joinStrings(":", [joinStrings(",", [ev_coords[0], ev_coords[len(ev_coords)-1]]), 
                                            joinStrings(",", ev_coords)])

def get_ev_strand(ev_counts):
    ev_strand = 0
    for ev_count in ev_counts:
        sign = 0 if ev_count[0] == "=" else 1 if ev_count[0] == "+" else -1
        ev_strand = ev_count[1] * sign
    ev_strand = None if ev_count[0] == 0 else 1 if ev_count[0] == "+" else -1
                
    if ev_strand == 0: return None
    elif ev_strand > 0: ev_strand = "+"
    elif ev_strand < 0: ev_strand = "-"
    
    return ev_strand

parser = argparse.ArgumentParser(description='Hola\n')

parser.add_argument('-f', "--file", dest="file", metavar='file', type=str, default = None,
help='sam-file')
parser.add_argument('-fl', "--filter", dest="filter", metavar='filter', type=int, default = 10,
help='sam-file')

args = parser.parse_args()

file = args.file
filter_thresh = args.filter

if(file==None):
    f = fileinput.input()
else:
    f = open(file).readlines()

jk_uniq = {}
l_tot = {}
r_tot = {}
lr_uniq = {}
known_cell_id = {}

for line in f:
    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()

    chr = temp[0]
    l = int(temp[1])
    r = int(temp[2])
    strand = temp[3]
    cell_id = temp[4]
    if cell_id not in known_cell_id:
        known_cell_id[cell_id] = 0
    known_cell_id[cell_id] += 1
    umi_l = temp[5]
    umi_r = temp[6]

    jk_id = (chr, l, r)
    
    if jk_id not in jk_uniq:
        jk_uniq[jk_id] = True
        
        pos_inc(chr, l, lr_uniq, 1)
        pos_inc(chr, r, lr_uniq, -1)

    pos_append(cell_id, chr, (l, strand), umi_l, l_tot)
    pos_append(cell_id, chr, (r, strand), umi_r, r_tot)

    
events = {}
for chr in lr_uniq:   
    lr_uniq[chr][0] = 0
    lr_sorted_keys = sorted(list(lr_uniq[chr].keys()))
    event = []
    pos = []
    for i in range(1, len(lr_sorted_keys)):
        lr_uniq[chr][lr_sorted_keys[i]] = lr_uniq[chr][lr_sorted_keys[i]] + lr_uniq[chr][lr_sorted_keys[i - 1]]
        event.append(lr_uniq[chr][lr_sorted_keys[i]])
        pos.append(lr_sorted_keys[i])
        
        if event[len(event)-1] == 0:
            event = classify_event(event)
            event_append(chr, event, pos, events)
            event = []
            pos = []
    del(lr_uniq[chr][0])


cell_clusters = {}
for cell_id in known_cell_id:
    cell_clusters[cell_id] = {cell_id:1}    

for chr in events:
    for ev_class in events[chr]:
        for ev_coords in events[chr][ev_class]:
            for cluster in cell_clusters:
                ev_known = False
                ev_strand = "none"
                ev_balanced = "none"
                ev_region = "none"
                psi = "none"

                if ev_class in ["cass"]:
                    ev_counts = get_ev_counts(cell_clusters[cluster], ev_coords, l_tot, r_tot)
                    if ev_counts == None:
                        continue
                    
                    ev_region = get_ev_region(ev_coords)

                    ev_strand = get_ev_strand(ev_counts)
                    if ev_strand == None:
                        continue

                    a = ev_counts[0][1]
                    b = ev_counts[1][1]
                    c = ev_counts[2][1]
                    d = ev_counts[3][1]

                    ev_balance = abs(a-b)-abs(d-c)
                    ev_shift = round(2*b/(b+c),5) 

                    psi = round((b+c)/(a+d),5)
                    
                    ev_counts = int((ev_counts[0][1]+ev_counts[1][1]+ev_counts[2][1]+ev_counts[3][1])/2)
                    if ev_counts < filter_thresh: 
                        continue

                    ev_known = True

                if ev_class in ["A5"]:
                    ev_counts = get_ev_counts(cell_clusters[cluster], ev_coords, l_tot, r_tot)
                    if ev_counts == None:
                        continue
                        
                    ev_region = get_ev_region(ev_coords)

                    ev_strand = get_ev_strand(ev_counts)
                    if ev_strand == None:
                        continue

                    a = ev_counts[0][1]
                    b = ev_counts[1][1]
                    c = ev_counts[2][1]

                    ev_balance = a-b-c 
                    ev_shift = "none"

                    psi = round(b/(b+c),5)
                    
                    ev_counts = int((ev_counts[0][1]+ev_counts[1][1]+ev_counts[2][1])/2)
                    if ev_counts < filter_thresh:
                        continue

                    ev_known = True

                if ev_class in ["A3"]:
                    ev_counts = get_ev_counts(cell_clusters[cluster], ev_coords, l_tot, r_tot)
                    if ev_counts == None:
                        continue

                    ev_region = get_ev_region(ev_coords)

                    ev_strand = get_ev_strand(ev_counts)
                    if ev_strand == None:
                        continue

                    a = ev_counts[0][1]
                    b = ev_counts[1][1]
                    c = ev_counts[2][1]

                    ev_balance = c-b-a
                    ev_shift = "none"

                    psi = round(a/(a+b),5)
                    
                    ev_counts = int((ev_counts[0][1]+ev_counts[1][1]+ev_counts[2][1])/2)
                    if ev_counts < filter_thresh:
                        continue
                    
                    ev_known = True

                if ev_known:
                    print(
                            cluster,
                            chr, 
                            ev_coords[0],
                            ev_coords[len(event)-1],
                            ev_strand,
                            ev_class,
                            ev_counts,
                            ev_balance,
                            #ev_shift,
                            psi,
                            ev_region
                    )


