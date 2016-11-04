import argparse
import pysam
from itertools import izip

parser = argparse.ArgumentParser(description =
    "This program accepts two BAM files each containing the same reads, but " +
    "aligned to different reference genomes. It compares the two alignemnt " +
    "scores of each read and assigns the read to the reference to which it " +
    "aligns better. If the read aligns equally well (see help menu) to both " +
    "references, the read will be classified as ambiguous. This program " +
    "expects the input files to have the same reads, whether they align or " +
    "not, in the same order. It also makes decisions based on the AS and YS " +
    "tags from the Bowtie2 aligner, and will fail if these are not present.")

parser.add_argument('--input1', action = 'store', metavar = 'FILE',
                    help = 'Input BAM file 1')
parser.add_argument('--input2', action = 'store', metavar = 'FILE',
                    help = 'Input BAM file 2')
parser.add_argument('--output1', action = 'store', metavar = 'FILE',
                    help = 'Output BAM file for records from input file 1')
parser.add_argument('--output2', action = 'store', metavar = 'FILE',
                    help = 'Output BAM file for records from input file 2')
parser.add_argument('--ambiguous1', action = 'store', metavar = 'FILE',
                    help = 'Output BAM file for ambiguous input-file-1 records')
parser.add_argument('--ambiguous2', action = 'store', metavar = 'FILE',
                    help = 'Output BAM file for ambiguous input-file-2 records')
parser.add_argument('--threshold_max', action = 'store', metavar = 'Y',
                    type = int, default = 0,
                    help = ('If the difference in scores falls within the ' +
                            'range [X, Y] (inclusive), classify the reads as' +
                            'ambiguous. (default = 0)'))
parser.add_argument('--threshold_min', action = 'store', metavar = 'X',
                    type = int, default = 0,
                    help = ('If the difference in scores falls within the ' +
                            'range [X, Y] (inclusive), classify the reads as' +
                            'ambiguous. (default = 0)'))
args = parser.parse_args()

if int(args.threshold_max) < int(args.threshold_min):
    exit('Invalid threshold parameters. Max must not be less than min.\n' +
         'Max: ' + args.threshold_max + ', min: ' + args.threshold_min)

tmax = int(args.threshold_max)
tmin = int(args.threshold_min)

with pysam.AlignmentFile(args.input1, "rb") as input_file1,\
     pysam.AlignmentFile(args.input2, "rb") as input_file2,\
     pysam.AlignmentFile(args.output1, "wb",
                         template = input_file1) as output_file1,\
     pysam.AlignmentFile(args.output2, "wb",
                         template = input_file2) as output_file2,\
     pysam.AlignmentFile(args.ambiguous1, "wb",
                         template = input_file1) as ambiguous_file1,\
     pysam.AlignmentFile(args.ambiguous2, "wb",
                         template = input_file2) as ambiguous_file2:

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

                if diff <= tmax and diff >= tmin:
                    ambiguous_file1.write(read1)
                    ambiguous_file2.write(read2)
                elif diff > tmax:
                    output_file1.write(read1)
                else:
                    output_file2.write(read2)

            elif read1.is_proper_pair:
                # Read 2 not paired properly
                output_file1.write(read1)

            elif read2.is_proper_pair:
                # Read 1 not paired properly
                output_file2.write(read2)

            else:
                ambiguous_file1.write(read1)
                ambiguous_file2.write(read2)

        # Single-read fragments
        else:

            if read1.is_unmapped and read2.is_unmapped:
                ambiguous_file1.write(read1)
                ambiguous_file2.write(read2)

            elif read1.is_unmapped:
                # Read 2 is mapped
                output_file2.write(read2)

            elif read2.is_unmapped:
                # Read 1 is mapped
                output_file1.write(read1)

            else:
                score1 = int(read1.get_tag("AS"))
                score2 = int(read2.get_tag("AS")) 

                diff = score1 - score2
            
                if diff <= tmax and diff >= tmin:
                    ambiguous_file1.write(read1)
                    ambiguous_file2.write(read2)
                elif diff > tmax:
                    output_file1.write(read1)
                else:
                    output_file2.write(read2)
