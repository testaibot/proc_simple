#!/usr/bin/python

from __future__ import print_function
import sys
import argparse
import ntpath
import os
import fileinput
from itertools import zip_longest


known_events = {(1,0):"jk",
                (2,1,0):"A5", 
                (1,2,0):"A3", 
                (2,1,2,0):"cass",
                (2,1,2,1,2,0):"mut_ex",
                (3,2,3,2,3,0):"mut_ex_cass"
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

def pos_append(chr, k, d, val):
    val = tuple(val)
    if chr not in d:
        d[chr] = {}
        
    if k not in d[chr]:
        d[chr][k] = {}
        
    if val not in d[chr][k]:
        d[chr][k][val] = 0
        
    d[chr][k][val] += 1
  

def classify_event(ev):
    ev = tuple(event)
    if ev in known_events:
        return known_events[ev]
    
    return joinStrings(",",ev)


def get_strand_counts(chr, pos, l_count, umi = True):
    count_strand = [0,0]
    if (pos, "+") in l_count[chr]:
        count_strand[0] = len(l_count[chr][(pos, "+")]) if umi else sum(l_count[chr][(pos, "+")].values())
    
    if (pos, "-") in l_count[chr]:
        count_strand[1] = len(l_count[chr][(pos, "-")]) if umi else sum(l_count[chr][(pos, "-")].values())
        
    return count_strand

def get_ev_counts(coord, l_tot, r_tot):
    ev_counts = []
    for pos in coord:
        l_counts = get_strand_counts(chr, pos, l_tot)
        r_counts = get_strand_counts(chr, pos, r_tot)
        l_strand = "=" if l_counts[0] == l_counts[1] else "+" if l_counts[0] > l_counts[1] else "-"
        r_strand = "=" if r_counts[0] == r_counts[1] else "+" if r_counts[0] > r_counts[1] else "-"
                    
        l_counts = sum(l_counts)
        r_counts = sum(r_counts)
        if r_counts == l_counts:
            raise Exception('Error: expected simple linear element. Input: composite')
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

args = parser.parse_args()

file = args.file

if(file==None):
    f = fileinput.input()
else:
    f = open(file).readlines()

jk_uniq = {}
l_tot = {}
r_tot = {}
lr_uniq = {}
lr_tot = {}

for line in f:
    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()

    chr = temp[0]
    l = int(temp[1])
    r = int(temp[2])
    strand = temp[3]
    cell_id = temp[4]
    m_l = temp[5]
    m_r = temp[6]

    jk_id = (chr, l, r)
    
    if jk_id not in jk_uniq:
        jk_uniq[jk_id] = True
        
        pos_inc(chr, l, lr_uniq, 1)
        pos_inc(chr, r, lr_uniq, -1)

    pos_append(chr, (l, strand), l_tot, (cell_id, m_l))
    pos_append(chr, (r, strand), r_tot, (cell_id, m_r))

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
            pos_append(chr, event, events, pos)
            event = []
            pos = []
    del(lr_uniq[chr][0])

for chr in events:
    for ev_class in events[chr]:
        for ev_coords in events[chr][ev_class]:
            ev_known = False
            
            ev_strand = "none"
            ev_balanced = "none"
            ev_region = "none"
            psi = "none"
            
            if ev_class in ["cass"]:
                ev_known = True
                ev_counts = get_ev_counts(ev_coords, l_tot, r_tot)
                
                ev_region = get_ev_region(ev_coords)
                
                ev_strand = get_ev_strand(ev_counts)
                if ev_strand == None: ev_known = False
                    
                a = ev_counts[0][1]
                b = ev_counts[1][1]
                c = ev_counts[2][1]
                d = ev_counts[3][1]
                
                ev_balance = abs(a-b)-abs(d-c)
                ev_shift = round(2*b/(b+c),5) 

                psi = round((b+c)/(a+d),5)
                if a+d < 10: ev_known = False
                
                ev_counts = int((ev_counts[0][1]+ev_counts[1][1]+ev_counts[2][1]+ev_counts[3][1])/2)
                
            if ev_class in ["A5"]:
                ev_known = True
                ev_counts = get_ev_counts(ev_coords, l_tot, r_tot)
                
                ev_region = get_ev_region(ev_coords)
                
                ev_strand = get_ev_strand(ev_counts)
                if ev_strand == None: ev_known = False
                    
                a = ev_counts[0][1]
                b = ev_counts[1][1]
                c = ev_counts[2][1]
                
                ev_balance = a-b-c 
                ev_shift = "none"

                psi = round(b/(b+c),5)
                if a < 10: ev_known = False
                
                ev_counts = int((ev_counts[0][1]+ev_counts[1][1]+ev_counts[2][1])/2)
                
            if ev_class in ["A3"]:
                ev_known = True
                ev_counts = get_ev_counts(ev_coords, l_tot, r_tot)
                
                ev_region = get_ev_region(ev_coords)
                
                ev_strand = get_ev_strand(ev_counts)
                if ev_strand == None: ev_known = False
                    
                a = ev_counts[0][1]
                b = ev_counts[1][1]
                c = ev_counts[2][1]
                
                ev_balance = c-b-a
                ev_shift = "none"

                psi = round(a/(a+b),5)
                if c < 10: ev_known = False
                
                ev_counts = int((ev_counts[0][1]+ev_counts[1][1]+ev_counts[2][1])/2)
                
            if ev_known:
                print(
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


