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

parser = argparse.ArgumentParser(description='Hola\n')
parser.add_argument('-f', "--file", dest="file", metavar='file', type=str, default = None, help='annotation file')
args = parser.parse_args()

file = args.file

if(file==None):
    f = fileinput.input()
else:
    f = open(file).readlines()


gene_regions = {}
for line in f:
    
    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()
    
    chr = temp[0]
    start = int(temp[1])
    end = int(temp[2])
    strand = temp[3]
    strand_int = 1 if strand == "+" else -1
    gene = temp[4]
    trans = temp[5]
    
    if chr not in gene_regions:
        gene_regions[chr] = {}
    
    if gene not in gene_regions[chr]:
        gene_regions[chr][gene] = [start, end, strand_int]
    else:
        if gene_regions[chr][gene][0] < start:
            gene_regions[chr][gene][0] = start
        if gene_regions[chr][gene][1] < end:
            gene_regions[chr][gene][1] = end
        gene_regions[chr][gene][2] += strand_int

gene_regions1 = {}
for chr in gene_regions:
    if chr not in gene_regions1:
        gene_regions1[chr] = {}
    
    for gene in gene_regions[chr]:
        start = gene_regions[chr][gene][0]
        end = gene_regions[chr][gene][1]
        n_frag = gene_regions[chr][gene][2]
        if start not in gene_regions1[chr]:
            gene_regions1[chr][start] = {}

        if end not in gene_regions1[chr]:
            gene_regions1[chr][end] = {}

        if n_frag not in gene_regions1[chr][start]:
            gene_regions1[chr][start][n_frag] = {}
        gene_regions1[chr][start][n_frag][gene] = True

        if n_frag not in gene_regions1[chr][end]:
            gene_regions1[chr][end][n_frag] = {}
        gene_regions1[chr][end][n_frag][gene] = True
                     
genes_nonoverhang = {}
for chr in gene_regions1:
    sorted_keys = sorted(list(gene_regions1[chr].keys()))
    prev_key = None
    prev_val = None
    
    if chr not in genes_nonoverhang:
        genes_nonoverhang[chr] = {}
    
    for pos in sorted_keys:
        #if len(gene_regions1[chr][pos]) > 1: continue
        max_val_sign = None
        max_val = 0
        for key in gene_regions1[chr][pos]:
            key = int(key)
            if key == 0: continue
            sign = key/abs(key)
            val = abs(key)
            if val > max_val:
                max_val = val
                max_val_sign = sign
        if max_val == 0: continue
                    
        best_val = int(max_val*max_val_sign)
        
        curr_key = list(gene_regions1[chr][pos][best_val].keys())

        if prev_key == None:
            prev_key = curr_key[0]
            prev_val = best_val
        else:
            if prev_key in curr_key:
                genes_nonoverhang[chr][prev_key] = True
                prev_key = None
                prev_val = None
            else:
                #if best_val > prev_val:
                prev_key = curr_key[0]
                #    prev_val = best_val
                #else: pass
            
sum_first_pass = 0
for chr in genes_nonoverhang:
    sum_first_pass += len(genes_nonoverhang[chr])


gene_regions1 = {}
for chr in gene_regions:
    if chr not in gene_regions1:
        gene_regions1[chr] = {}
        
    for gene in gene_regions[chr]:
        if chr in genes_nonoverhang and gene in genes_nonoverhang[chr]:
            start = gene_regions[chr][gene][0]
            end = gene_regions[chr][gene][1]
            n_frag = gene_regions[chr][gene][2]
            if start not in gene_regions1[chr]:
                gene_regions1[chr][start] = {}

            if end not in gene_regions1[chr]:
                gene_regions1[chr][end] = {}

            if n_frag not in gene_regions1[chr][start]:
                gene_regions1[chr][start][n_frag] = {}
            gene_regions1[chr][start][n_frag][gene] = True

            if n_frag not in gene_regions1[chr][end]:
                gene_regions1[chr][end][n_frag] = {}
            gene_regions1[chr][end][n_frag][gene] = True
        
genes_nonoverhang = {}
for chr in gene_regions1:
    sorted_keys = sorted(list(gene_regions1[chr].keys()))
    prev_key = None
    prev_val = None
    
    if chr not in genes_nonoverhang:
        genes_nonoverhang[chr] = {}
    
    for pos in sorted_keys:
        #if len(gene_regions1[chr][pos]) > 1: continue
        max_val_sign = None
        max_val = 0
        for key in gene_regions1[chr][pos]:
            key = int(key)
            if key == 0: continue
            sign = key/abs(key)
            val = abs(key)
            if val > max_val:
                max_val = val
                max_val_sign = sign
        if max_val == 0: continue
                    
        best_val = int(max_val*max_val_sign)
        
        curr_key = list(gene_regions1[chr][pos][best_val].keys())

        if prev_key == None:
            prev_key = curr_key[0]
            prev_val = best_val
        else:
            if prev_key in curr_key:
                genes_nonoverhang[chr][prev_key] = True
                prev_key = None
                prev_val = None
            else:
                #if best_val > prev_val:
                prev_key = curr_key[0]
                #    prev_val = best_val
                #else: pass
        
sum_second_pass = 0
for chr in genes_nonoverhang:
    sum_second_pass += len(genes_nonoverhang[chr])


if sum_second_pass - sum_first_pass == 0:
    for chr in gene_regions1:
        sorted_keys = sorted(list(gene_regions1[chr].keys()))
        
        for pos_start, pos_end in grouper(sorted_keys, 2):
            gene_start = list(list(gene_regions1[chr][pos_start].values())[0].keys())[0]
            gene_end = list(list(gene_regions1[chr][pos_end].values())[0].keys())[0]
            if gene_start == gene_end:
                print(chr, pos_start, pos_end, gene_start, sep = "\t")
            else: 
                raise Exception("Unknown start-end match error")
else:
    raise Exception("Unknown error second pass did not converge")
