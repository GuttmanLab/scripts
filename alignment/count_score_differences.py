import argparse
import pysam
from argparse import RawTextHelpFormatter
from collections import Counter
from itertools import izip

parser = argparse.ArgumentParser(formatter_class = RawTextHelpFormatter, 
    description =
    "This program accepts two BAM files each containing the same reads, but aligned\n" +
    "to different reference genomes. It keeps a tally of the number of times it\n" +
    "sees each difference, e.g., a difference of 0 was seen 10 times, a difference\n" +
    "of 1 was seen 2 times, etc. These differences are then printed to STDOUT as a\n" +
    "tab-delimited text file, with the difference magnitude in the first column and\n" +
    "number of occurrences in the second:\n\n" +
    "-40\t10\n-39\t20\n-38\t1\n-37\t2\n\n" + 
    "and so on. This program examines tags from the Bowtie2 aligner, and will fail\n" +
    "if these are not present. Rows with zero occurrences are skipped.")

parser.add_argument('--input1', action = 'store', metavar = 'FILE',
                    help = 'Input BAM file 1')
parser.add_argument('--input2', action = 'store', metavar = 'FILE',
                    help = 'Input BAM file 2')
args = parser.parse_args()

count = Counter()

with pysam.AlignmentFile(args.input1, "rb") as input_file1,\
     pysam.AlignmentFile(args.input2, "rb") as input_file2:

    for read1, read2 in izip(input_file1.fetch(until_eof = True),
                             input_file2.fetch(until_eof = True)):

        if read1.query_name != read2.query_name:
            exit('Reads in input files do not match: ' +
                 read1.query_name + ', ' + read2.query_name)

        # Paired-end fragments
        if read1.is_paired and read2.is_paired:

            if read1.is_proper_pair and read2.is_proper_pair:
                score1 = int(read1.get_tag("AS")) + int(read1.get_tag("YS"))
                score2 = int(read2.get_tag("AS")) + int(read2.get_tag("YS"))
                diff = score1 - score2
                count[diff] += 1

        # Single-read fragments
        else:

            if not read1.is_unmapped and not read2.is_unmapped:
                score1 = int(read1.get_tag("AS"))
                score2 = int(read2.get_tag("AS")) 
                diff = score1 - score2
                count[diff] += 1

for (diff, num) in sorted(count.iteritems()):
    print str(diff) + "\t" + str(num) 
