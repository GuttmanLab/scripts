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

    parser.add_argument('-i', '--input', metavar = 'FILE', action = 'store',
        help = 'The input barcode-groups file.')

    parser.add_argument('-o', '--output', metavar = 'FILE', action = 'store',
        help = 'The output raw-contacts file.')

    parser.add_argument('-c', '--chromosome', metavar = 'STRING',
        action = 'store', help = 'The chromosome to consider.')

    parser.add_argument('--max_cluster_size', metavar = 'INT', type = int,
        action = 'store', help = 'The maximum read-cluster size to ' +
        'consider. Ignore clusters larger than INT.')

    parser.add_argument('-r', '--resolution', metavar = 'INT', type = int,
        action = 'store', help = 'The resolution in bp.')

    return parser.parse_args()

def get_contacts(args):
    contacts = defaultdict(lambda : Counter())
    total_groups = 0
    with open(args.input, 'r') as f:
        for line in f:
            reads = line.split()[1:]
            if len(reads) > args.max_cluster_size:
                continue
            bins = set()
            for read in reads:
                chrom, position = read.split(':')
                if chrom == args.chromosome:
                    read_bin = int(position) // args.resolution
                    bins.add(read_bin)
            for bin1, bin2 in combinations(bins, 2):
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
