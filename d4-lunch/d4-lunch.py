#!/usr/bin/env python3

import sys

import numpy

# head -n 500 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct > test_data.gct
# ./d4-lunch.py gene_tissue.tsv GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct

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
    key = (fields[0])
    # initalize dict from key with list to hold samples, use [] to make it a list
    relevant_samples[key] = fields[2]
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
for tissue,samples in tissue_samples.items():
    tissue_columns.setdefault(tissue, [])
    for sample in samples: 
        if sample in header:
            # what position is the sample we are looking at right now
            position = header.index(sample)
            tissue_columns[tissue].append(position)

# print(tissue_columns)



# Q5
# need to count samples per tissue, determine which has the most (increasing variable value starting at 0)
maxValue = 0
maxValueKey = " "
for tissue,sample in tissue_columns.items():
    if len(sample) > maxValue: 
        maxValue = len(sample)
        maxValueKey = tissue
        
# print(maxValueKey)

# need to count samples per tissue to find which has the least, start mnValue at very high number so all are smaller and it can work down to least
minValue = 20000000
minValueKey = " "
for tissue,sample in tissue_columns.items():
    if len(sample) < minValue: 
        minValue = len(sample)
        minValueKey = tissue
        
# print(minValueKey)


#The tissue type with the largest number of samples is muscle - skeletal and the fewest is Cells - Leukemia cell line (CML).

# Q6 (got moved to advanced/optional, not finishing but did start)

#relevant_tissues to find the genes and tissue_columns to find the columns, then pull out expression values from tmp

#read in file we need
f = open("test_data.gct", "r")
for l in f: 
    l = l.strip().split("\t")
    # print(l[0])

    geneName = l[0]

# the keys of the relvant_samples dictionary are the gene names
    if geneName in relevant_samples.keys():
        myTissue = relevant_samples[geneName] # find the tissue your gene is expressed in
        print(tissue_columns[myTissue])


