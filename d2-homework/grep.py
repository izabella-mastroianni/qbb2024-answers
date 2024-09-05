#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )
pattern = sys.argv[1]


for line in my_file:
    line = line.rstrip("\n")
    if pattern in line:
        print( line )
  
my_file.close()