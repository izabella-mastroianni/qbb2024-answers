#!/usr/bin/env python3

import sys

import numpy

# Q1

# get gene tissue file name (gene_tissue.tsv)
filename = sys.argv[1]
# open file - did before for loop (which would autoclose it) so have to explicitly close it later
fs = open(filename, mode ='r')
# create dict to hold samples for gene-tissue pairs
relevant_samples = {}
# step through file
for line in fs: 
    # split line into fields: remove end of line character and tab
    fields = line.rstrip("\n").split("\t")
    # create key from gene and tissue in () to make tuple so we make it an immutable
    key = (fields[0], fields[2])
    # initalize dict from key with list to hold samples, use [] to make it a list
    relevant_samples[key] = []
fs.close()

# print(relevant_samples)

# Q2

# get metadata file name (DS)
filename = sys.argv[2]
# open file - did before for loop (which would autoclose it) so have to explicitly close it later
fs = open(filename, mode ='r')
# skip line
fs.readline()
# create dict to hold samples for tissue name
tissue_samples = {}
# step through file
for line in fs: 
    # split line into fields: remove end of line character and tab
    fields = line.rstrip("\n").split("\t")
    # create key from gene and tissue 
    key = fields[6]
    value = fields[0]
    # initalize dict from key with list to hold samples, use [] to make it a list
    tissue_samples.setdefault(key, [])
    # append to empty list specified above
    tissue_samples[key].append(value)
fs.close()

# print(tissue_samples)

# Q3 
# get metadata file name (tmp gene expression file)
filename = sys.argv[3]
# open file 
fs = open(filename, mode ='r')
# skip 2 lines
fs.readline()
fs.readline()
header = fs.readline().rstrip("\n").rsplit("\t")
# header is a list that contains all the sample names in the RNAseq file
header = header[2:]

# Q4

tissue_columns = {}
# key = tissue and value = samples set in the following line with items function 
for tissue, samples in tissue_samples.items():
    tissue_columns.setdefault(tissue, [])
    for sample in samples: 
        if sample in header:
            position = header.index(sample)
            tissue_columns[tissue].append(position)

#need to sort tissue_columns by tissue or somehow count samples per tissue
sorted_tissue_columns = tissue_columns.sort_values(by='tissue')
print(sorted_tissue_columns)

# Q5

#The tissue type with the largest number of samples is ___ and the fewest is ___.

# Q6



