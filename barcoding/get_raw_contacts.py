import sys
from collections import Counter
from collections import defaultdict
from itertools import combinations

barcode_groups = sys.argv[1]
bin_size = int(sys.argv[2])
max_size = int(sys.argv[3])
output_dir = sys.argv[4]

LOG_STEP = 100000 # Output progress every LOG_STEP barcodes

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']

#if human:
#    chromosomes.extend(['chr20', 'chr21', 'chr22'])

for chromosome in chromosomes:
    contacts = defaultdict(lambda : Counter())
    total_groups = 0
    total_singles = Counter()

    with open(barcode_groups, 'r') as f:
        file_barcode_count = 0   # Number of barcodes (lines) in file
        chrom_barcode_count = 0  # Number of barcodes which contain a read
                                 # on this chromosome
        for line in f:
            file_barcode_count += 1
            chrom_barcode_count += 1
            reads = line.split()[2:]

            # Skip "mega-clusters"
            if len(reads) > max_size:
                continue

            bins = set() # Unique bins on this chromosome for this barcode

            for read in reads:
                chrom, position = read.split(':')
                if chrom == chromosome:
                    read_bin = int(position) // bin_size
                    bins.add((chrom, read_bin))

            # Tally contacts in this line
            for bin1, bin2 in combinations(bins, 2):
                contacts[bin1][bin2] += 1
                contacts[bin2][bin1] += 1
                total_singles[bin1] += 1
                total_singles[bin2] += 1
                total_groups += 1

            # Log progress
            if file_barcode_count % LOG_STEP == 0:
                print(chromosome + ":\t" + str(file_barcode_count))

    print("Chromosome: " + chromosome)
    print("Total number of barcodes in file: " + str(file_barcode_count))
    print("Total number of barcodes in chromosome: " + str(chrom_barcode_count))
    print("Total number of pairwise contacts within chromosome: " + str(total_groups))

    # Output to file
    with open(output_dir + '/' + barcode_groups + '.' + str(bin_size) + '.' + str(max_size) + '.' + chromosome + '.raw_contacts', 'w') as f:

        for (chrom1, position1), inner_dict in contacts.items():
            for (chrom2, position2), count in inner_dict.items():
                f.write(chrom1 + ' ' + str(position1) + '\t' +
                        chrom2 + ' ' + str(position2) + '\t' +
                        str(count) +'\t' +
                        str(total_singles[(chrom1, position1)]) +
                        '\n')
