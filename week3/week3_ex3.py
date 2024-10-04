#!/usr/bin/env python3

import numpy


## Exercise 3 ##


# specify VCF file and output AF.txt file
vcf_file = open('biallelic.vcf')
af_output_file = open('AF.txt', 'w')
dp_output_file = open('DP.txt', 'w')

# add a header to the output file for af, only contains allele frequencies
af_output_file.write('Allele_Frequency\n')

# add header for depth of reads file, only contains read depths
dp_output_file.write('Read_Depth\n')




for line in vcf_file:
    if line.startswith('#'):
        # look for header with allele frequency info
        if line.startswith('#INFO'):
            print(line.strip())
            # skip header lines
        continue
    # split lines into fields
    fields = line.rstrip('\n').split('\t')

# For AF
    # INFO field is the 8th column, python counting starts at 0
    info_field = fields[7]
    # set af_value to none outside loop
    af_value = None

    # since AF is allele frequency and what we want, find that field and pull out the values
    for entry in info_field.split(';'):
        if entry.startswith('AF='):
            af_value = entry.split('=')[1]
            # stop when it finds an AF
            break 
    # when find AF, write it into output file
    if af_value is not None:
        af_output_file.write(af_value + '\n')

# For DP
        dp_value = None

    for sample in fields[9:]:
        sample_data = sample.split(':')
        dp_value = sample_data[2]
        if dp_value is not None: 
            dp_output_file.write(dp_value + '\n')




# close files
vcf_file.close()
af_output_file.close()
dp_output_file.close()

