#!/usr/bin/env python3

#steps/pseudocode for question 1:
    #unzip file (gunzip)
    #make copy of file (cp)
    #open file to check what data looks like (vim then :s unwrap or less -S)
    #open file (sys.argv[1]) - says to open what was written second in the command line
    #skip first 2 lines (fs.readline() twice - saying read that line and don't do anything with it)
    #split column headers (separated by tabs and skip first two entries)
    #create way to store gene names 
    #create way to store expression values
    #for each line
        #split line
        #save field 0 into gene IDs
        #save field 1 into gene names
        #save fields 2+ into expression values 

#Q1
import sys

import numpy

#open file
fs = open(sys.argv[1], mode='r')
#skip first 2 lines
fs.readline()
fs.readline()
#split column header by tabs and strip new line character and remove first two lines
line = fs.readline()
fields = line.strip("\n").split("\t")
tissues = fields[2:]
#create way to hold gene names, gene IDs, and expression values
gene_names = []
gene_IDs = []
expression = []
#for each line
for line in fs: 
    #split line
    fields = line.strip("\n").split("\t")
    #save field 0 into gene IDs
    gene_IDs.append(fields[0])
    #save field 1 into gene names
    gene_names.append(fields[1])
    #save fields 2+ into expression values 
    expression.append(fields[2:])
fs.close()

#Q2
tissues = numpy.array(tissues)
gene_IDs = numpy.array(gene_IDs)
gene_names = numpy.array(gene_names)
expression = numpy.array(expression, dtype = float)
#we had to specify the expression data type because it is reading it as a string since that is how it is defined in the orginal file
# print(gene_IDs)
# print(gene_names)
# print(expression)

#Q3 - told to skip

#Q4
expression10 = expression[:10, :]
mean_expression10 = numpy.mean(expression10, axis=1)
#print(mean_expression10)

#Q5

mean_expression = numpy.mean(expression)
median_expression = numpy.median(expression)
# print(mean_expression)
# print(median_expression)

#the mean gives you a better global view of the data because the median is so skewed by all of the 0 expression values

#Q6
#add 1 to all expression values, then continue using log transformed data
expression1 = expression + 1
log2_expression1 = numpy.log2(expression1)
mean_log2_expression1 = numpy.mean(log2_expression1)
median_log2_expression1 = numpy.median(log2_expression1)
# print(mean_log2_expression1)
# print(median_log2_expression1)

#now the median and mean are much closer to each other in value, with a difference of ~1. The median is almost the same between the transformed and non-transformed data, but the means are significantly different


#Q7
#numpy.sort returns a sorted array, it doesn't sort an array in place
cp_log2_expression1 = numpy.copy(log2_expression1)
sort_cp_log2_expression1 = numpy.sort(cp_log2_expression1, axis=1)
#print(numpy.sort(log2_expression1,axis=1)[:5,:])
#print(sort_cp_log2_expression1[0:5,:])
diff_array = sort_cp_log2_expression1[:,-1] - sort_cp_log2_expression1[:,-2]
#print(diff_array)

#Q8
where_diff_array = numpy.where(diff_array >=10)
#print(where_diff_array)
print(numpy.shape(where_diff_array))

#33 genes meet this requirement 
