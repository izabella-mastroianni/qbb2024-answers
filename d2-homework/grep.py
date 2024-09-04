#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )


for line in my_file:
    line = line.rstrip("\n")
    print( line )
  
my_file.close()