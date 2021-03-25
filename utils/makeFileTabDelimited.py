#! /usr/bin/env python
import sys
import pandas as pd
handle = open(sys.argv[1],"r")
readlines = [filter(None,line.replace("\n","").split(" ")) for line in handle.readlines()]
print('\n'.join(map('\t'.join,readlines)))

