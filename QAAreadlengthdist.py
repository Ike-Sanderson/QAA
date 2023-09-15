#!/usr/bin/env python
# Author: Ike Sanderson <ikes@uoregon.edu>

#version 0.1
#This file inputs space separated text files produced by htseq and 
#read them to output a histogram with base pair length on the X and 
#read count on the Y.

import argparse
import gzip

parser = argparse.ArgumentParser(description="A program to open a space separated text file, output readcount PNG")
parser.add_argument("-f", "--fn", type=str, help="filename.txt", required=True)
parser.add_argument("-o", type=str, help="outfile.png name", required=True)
args = parser.parse_args()

# Global variables
R1 = args.fn1
outfile = args.o

