from collections import Counter
from collections import defaultdict
from itertools import combinations
import argparse

def main():
    args = parse_arguments()
    contacts = get_contacts(args)
    write_contacts_to_file(contacts, args)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates an intrachromosomal raw-contacts file ' +
        'from a barcode-groups file.')

    parser.add_argument('-i', '--input',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The input barcode-groups file.')

    parser.add_argument('-o', '--output',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The output raw-contacts file.')

    parser.add_argument('-c', '--chromosome',
                        metavar = 'STRING',
                        action = 'store',
                        help = 'The chromosome to consider.')

    parser.add_argument('--min_cluster_size', 
                        metavar = 'INT',
                        type = int,
                        default = 2,
                        action = 'store',
                        help = 'The minimum read-cluster size to consider. ' + 
                               'Ignore clusters smaller than INT (default = 2).')

    parser.add_argument('--max_cluster_size',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        help = 'The maximum read-cluster size to consider. ' +
                               'Ignore clusters larger than INT.')

    parser.add_argument('-r', '--resolution', 
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        help = 'The resolution in bp.')

    parser.add_argument('--granularity',
                        metavar = 'INT',
                        type = int,
                        default = 1000,
                        action = 'store',
                        help = 'The resolution at which to consider two ' +
                               'reads to be equal (default = 1000).')

    return parser.parse_args()


def get_contacts(args):

    contacts = defaultdict(lambda : Counter())
    total_groups = 0

    with open(args.input, 'r') as f:
        for line in f:
            reads = line.split()[1:]

            # Check if cluster size is valid
            if not args.min_cluster_size <= len(reads) <= args.max_cluster_size:
                continue

            granules = set()

            # For reads on this chromosome, granulate the read positions.
            for read in reads:
                chrom, position = read.split(':')
                if chrom == args.chromosome:
                    granule = int(position) // args.granularity
                    granules.add(granule)

            # Scale-up the granulated read positions to the bin-resolution and
            for granule1, granule2 in combinations(granules, 2):
                bin1 = (granule1 * 1000) // args.resolution
                bin2 = (granule2 * 1000) // args.resolution
                if bin1 != bin2:
                    contacts[bin1][bin2] += 1
                    contacts[bin2][bin1] += 1

    return contacts
                
def write_contacts_to_file(contacts, args):
    with open(args.output, 'w') as f:
        for position1, second_positions in contacts.items():
            for position2, count in second_positions.items():
                f.write(args.chromosome + ' ' + str(position1) + '\t' +
                        args.chromosome + ' ' + str(position2) + '\t' +
                        str(count) + '\n')

if __name__ == "__main__":
    main()
