#make it take sys.argv
import sys
import os

#assign files to variables:
unique_erik = sys.argv[1]
unique_me = sys.argv[2]
intersect = sys.argv[3]
unique_erik_unmap = sys.argv[4]
unique_me_unmap = sys.argv [5]
intersect_unmap = sys.argv[6]

#output file
outfile = sys.argv[7]

#a function that count lines in the files, unless the first charatcher is a '#'
def count_lines(file):
    " count lines and return a interger "
    