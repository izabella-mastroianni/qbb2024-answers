#!/usr/bin/env python3

import numpy


## Exercise 3 ##

# make a list to ensure allele frequencies don't repeat/only get unique ones 
allele_frequencies = []

# list for read depth
read_depth = []

# specify VCF file and output AF.txt file
vcf_file = open('biallelic.vcf')
af_output_file = open('AF.txt', 'w')

# add a header to the output file for af
af_output_file.write('Allele_Frequency\n')



for line in vcf_file:
    if line.startswith('#'):
        # look for header with allele frequency info
        if line.startswith('#INFO'):
            print(line.strip())
        continue

    fields = line.rstrip('\n').split('\t')
    # INFO field is the 8th column, python counting starts at 0
    info_field = fields[7]

    # since AF is allele frequency and what we want, find that field and pull out the values
    for entry in info_field.split(';'):
        if entry.startswith('AF='):
            af_value = entry.split('=')[1]
            # don't want af of 0, so exclude
            if af_value != '0':
                if af_value not in allele_frequencies:
                    allele_frequencies.append(af_value)

for af in allele_frequencies:
    af_output_file.write(f"{af}\n")   



# close files
vcf_file.close()
af_output_file.close()

