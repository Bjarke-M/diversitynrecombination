#make it take sys.argv
import sys
import os

# Assign files to variables
unique_erik = sys.argv[1]
unique_me = sys.argv[2]
intersect = sys.argv[3]
unique_erik_unmap = sys.argv[4]
unique_me_unmap = sys.argv[5]
intersect_unmap = sys.argv[6]

# Output file
outfile = sys.argv[7]

def get_sample_name(file):
    # Split the file path by '/' to get the filename
    file_name = file.split('/')[-1]
    # Split the filename by '.' to remove the extension
    file_name = file_name.split('.')[0]
    return file_name

# A function that counts lines in the files unless the first character is a '#'
def count_lines(file):
    open_file = open(file, 'r')
    count = 0
    for line in open_file:
        if line[0] != '#':
            count += 1
    return count

def write_to_file(sample_name, count_unique_erik, count_unique_me, count_intersect, count_unique_erik_unmap, count_unique_me_unmap, count_intersect_unmap, outfile):
    open_outfile = open(outfile, 'a')
    open_outfile.write(sample_name + ',' + str(count_unique_erik) + ',' + str(count_unique_me) + ',' + str(count_intersect) + ',' + str(count_unique_erik_unmap) + ',' + str(count_unique_me_unmap) + ',' + str(count_intersect_unmap) + '\n')
    open_outfile.close()
    return outfile

write_to_file(get_sample_name(unique_erik), count_lines(unique_erik), count_lines(unique_me), count_lines(intersect), count_lines(unique_erik_unmap), count_lines(unique_me_unmap), count_lines(intersect_unmap), outfile)
