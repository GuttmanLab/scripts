from collections import defaultdict
import argparse
import pysam
import re

def main():
    args = parse_arguments()
    barcode_groups = get_barcode_groups(args)
    write_barcode_groups_to_file(barcode_groups, args)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a barcode-groups file from a BAM file.')
    parser.add_argument('-i', '--input',
                        metavar = "FILE",
                        action = "store",
                        help = "The input BAM file.")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        help = "The output barcode-groups file.")
    parser.add_argument('-n', '--num_barcodes',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        help = "The number of barcodes contained in the name " +
                               "of name of the BAM record.")
    return parser.parse_args()

def get_barcode_groups(args):
    barcode_groups = defaultdict(set)
    total = 0
    pattern = re.compile('::' + args.num_barcodes * '\[(\w+)\]')
    with pysam.AlignmentFile(args.input, "rb") as bamfile:
        for read in bamfile.fetch(until_eof = True):
            total += 1
            chromosome = read.reference_name
            position = read.reference_start
            name = read.query_name
            match = pattern.search(name)
            barcode_group = ".".join(match.groups())
            barcode_groups[barcode_group].add((chromosome, position))
    print_number_of_duplicates(barcode_groups, total)
    print "File has " + str(len(barcode_groups)) + " barcode-groups."
    return barcode_groups

def print_number_of_duplicates(barcode_groups, total):
    non_duplicates = 0
    for barcode_key, positions in barcode_groups.iteritems():
        non_duplicates += len(positions)
    print "Duplicates: " + str(total - non_duplicates)

def write_barcode_groups_to_file(barcode_groups, args):
    with open(args.output, 'w') as outfile:
        for barcode, positions in barcode_groups.iteritems():
            strings = [(lambda x: x[0] + ":" + str(x[1]))(pos) for pos in positions]
            outfile.write(barcode + "\t" + "\t".join(strings) + "\n")

if __name__ == "__main__":
    main()
