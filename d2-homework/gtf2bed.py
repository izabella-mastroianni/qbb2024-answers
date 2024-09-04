#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )


for line in my_file:

    if line.startswith("##"):
        continue 

    splitline = line.split("\t")
    splitline2 = splitline[8].replace('"','').split(";")
 
    print(splitline[0], splitline[3], splitline[4], splitline2[2].lstrip('gene_name '))


my_file.close()
