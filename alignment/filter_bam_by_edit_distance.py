#!/usr/bin/python
import argparse
import pysam

def main():
    args = parse_arguments()
    filter_reads(args)

def parse_arguments():

    parser = argparse.ArgumentParser(description =
            "This program removes reads from a BAM file according to the " +
            "number of mismatches and ambiguous bases they contain. This " +
            "program will not adjust the read flags for paired reads if one " +
            "read of the pair is removed and the other retained.")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output BAM file containing reads that pass')
    parser.add_argument('--edit_max', action = 'store', metavar = 'Y',
                        type = int, default = 0,
                        help = ('If the edit distance between the read ' +
                                'sequence and the reference sequence falls ' +
                                'within [X, Y], keep the read. Otherwise ' +
                                'remove it. Default = 0'))
    parser.add_argument('--edit_min', action = 'store', metavar = 'X',
                        type = int, default = 0,
                        help = ('If the edit distance between the read ' +
                                'sequence and the reference sequence falls ' +
                                'within [X, Y], keep the read. Otherwise ' +
                                'remove it. Default = 0'))
    args = parser.parse_args()

    if args.edit_max < args.edit_min:
        exit("Invalid edit-distance parameters. Max: " +
             str(args.edit_max) + ", min: " + str(args.edit_min))

    return args

def filter_reads(args):
    with pysam.AlignmentFile(args.input, "rb") as input_file,\
         pysam.AlignmentFile(args.output, "wb", template = input_file) as output_file:

        for read in input_file.fetch(until_eof = True):
            if not read.is_unmapped and has_valid_edit_distance(read,
                    args.edit_min, args.edit_max):
                output_file.write(read)

def has_valid_edit_distance(read, lo, hi):
    score = int(read.get_tag("NM"))
    return score >= lo and score <= hi

if __name__ == "__main__":
    main()
