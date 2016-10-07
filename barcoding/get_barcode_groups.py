from collections import defaultdict
import argparse
import pysam
import re

parser = argparse.ArgumentParser(description = 'Generates a "raw contacts" file from a BAM file.')

parser.add_argument('-i', '--input', metavar = FILE, action = "store",
                    help = "The input BAM file. The pysam library requires this to be sorted and indexed.")

parser.add_argument('-o', '--output', metavar = FILE, action = "store",
                    help = "The output raw-contacts file.")

parser.add_argument('-n', '--num_barcodes', metavar = N, action = "store", type = int, default = 5,
                    help = "The number of barcodes contained in the name of each BAM record. (default = 5)")

args = parser.parse_args()

# barcode_groups :: Dict Protocol (Dict Barcode (Set Position))
barcode_groups = defaultdict(lambda: defaultdict(set))
total = 0
multimappers = 0
pattern = re.compile('::' + args.num_barcodes * '\[(\w+)\]')

# Read through BAM file
with pysam.AlignmentFile(args.input, "rb") as samfile:
    for read in samfile.fetch():

        total += 1

        chromosome = read.reference_name
        position   = read.reference_start
        mapq       = read.mapping_quality
        name       = read.query_name

        # Skip multimappers, defined by MAPQ score
        if mapq < 2:
            multimappers += 1
            continue

        match = pattern.search(name)

        # Duplicates are skipped automatically since the datastructure is a set.
        protocol = match.groups()[0][0]
        barcode = ".".join(match.groups()[1:])
        barcode_groups[protocol][barcode].add(chromosome + ":" + str(position))

# Count the number of non-duplicates by getting the size of each
# protocol-barcode set
non_duplicates = 0
for protocol_key, barcode_dict in barcode_groups.iteritems():
    for barcode_key, position_set in barcode_dict.iteritems():
        non_duplicates += len(position_set)

print "Multimappers: " + str(multimappers)
print "Duplicates: " + str(total - non_duplicates - multimappers)

# Print everything to file
for protocol in barcode_groups:
    with open(args.output + '_' + protocol + '.barcode_groups', 'w') as outfile:
        for barcode, positions in barcode_groups[protocol].iteritems():
            outfile.write(protocol + "\t" + barcode + "\t" + "\t".join(positions) + "\n")
    print "Protocol " + protocol + " has " + str(len(barcode_groups[protocol])) + " barcode groups."
