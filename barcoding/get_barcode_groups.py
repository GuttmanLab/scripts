from collections import defaultdict
import pysam
import re
import sys

# Change as needed. This variable defines a regex pattern to match against.
NUMBER_BARCODES = 5

input_file = sys.argv[1]

# barcode_groups :: Dict Protocol (Dict Barcode (Set Position))
barcode_groups = defaultdict(lambda : defaultdict(set))
total = 0
multimappers = 0
pattern = re.compile('::' + NUMBER_BARCODES * '\[(\w+)\]')

# Read through BAM file
samfile = pysam.AlignmentFile(input_file, "rb")
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

samfile.close()

# Count the number of non-duplicates by getting the size of each
# protocol-barcode set
non_duplicates = 0
for protocol_key, barcode_dict in barcode_groups.iteritems():
    for barcode_key, position_set in barcode_dict.iteritems():
        non_duplicates += len(position_set)

print "Multimappers: " + str(multimappers)
print "Duplicates: " + str(total - non_duplicates - multimappers)

# Print everything to files
for protocol in barcode_groups:
    with open(input_file + '_' + protocol + '.barcode_groups', 'w') as outfile:
        for barcode, positions in barcode_groups[protocol].iteritems():
            outfile.write(protocol + "\t" + barcode + "\t" + "\t".join(positions) + "\n")
    print "Protocol " + protocol + " has " + str(len(barcode_groups[protocol])) + " barcode groups."
