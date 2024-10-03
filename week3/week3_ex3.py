#! usr/bin/env python3

import numpy


## Exercise 3 ##

# specify VCF file and output AF.txt file
vcf_file = open('biallelic.vcf')
output_file = open('AF.txt', 'w')

# add a header to the output file
output_file.write('Allele_Frequency\n')


for line in vcf_file:
    if line.startswith('#'):
        continue
    fields = line.rstrip('\n').split('\t')
    # INFO field is the 8th column, python counting starts at 0
    info_field = fields[7]


    # grab what you need from `fields`


# close files
vcf_file.close()
output_file.close()

