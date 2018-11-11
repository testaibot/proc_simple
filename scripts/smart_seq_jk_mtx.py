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


args = parser.parse_args()


file = "out_mtx"
if not os.path.exists(file):
    os.mkdir(file)

mtx = []
jk_names = {}
jk_idex = 1
cell_names = {}
cell_idex = 1

f = fileinput.input()
f_mtx = open(joinStrings([file, "/", 'mtx.mtx']), 'w+')

lines_count = 0
for line in f:
    lines_count +=1
    temp = line.split('\t')
    temp[-1] = temp[-1].rstrip()
       
    jk_id = temp[0]
    cell_id = temp[1]
    reads = temp[2]
    
    if jk_id not in jk_names:
        jk_names[jk_id] = jk_idex
        jk_idex += 1

    if cell_id not in cell_names:
        cell_names[cell_id] = cell_idex
        cell_idex += 1

    f_mtx.write(joinStrings([jk_names[jk_id],"\t", cell_names[cell_id],"\t", reads, "\n"]))


f_jk = open(joinStrings([file, "/", 'jk_names.txt']), 'w+')
f_cell = open(joinStrings([file, "/", 'cell_names.txt']), 'w+')
for item in jk_names:
    f_jk.write(item)
    f_jk.write("\n")

for item in cell_names:
    f_cell.write(item)
    f_cell.write("\n")

f_header = open(joinStrings([file, "/", 'header.txt']), 'w+')
f_header.write("%%MatrixMarket matrix coordinate real general")
f_header.write(joinStrings(["\n", jk_idex-1,"\t", cell_idex-1,"\t", lines_count, "\n"]))


