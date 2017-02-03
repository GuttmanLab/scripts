from collections import Counter
from collections import defaultdict
import argparse

def main():
    args = parse_arguments()
    (matrix, positions) = construct_matrix(args)
    write_matrix_to_file(matrix, positions, args)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a matrix text file from a raw-contacts file.')
    parser.add_argument('-i', '--input',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The input raw-contacts file.')
    parser.add_argument('-o', '--output',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The output matrix file.')
    parser.add_argument('-c', '--chromosome',
                        metavar = 'STRING',
                        action = 'store',
                        help = 'The name of the chromosome.')
    parser.add_argument('-r', '--resolution',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        help = 'The resolution of the raw-contacts file.')
    parser.add_argument('-a', '--assembly',
                        metavar = 'ASSEMBLY',
                        choices = ['mm9', 'mm10', 'hg19', 'hg38'],
                        default = 'mm9',
                        action = 'store',
                        help = 'The genome assembly')
    return parser.parse_args()

def construct_matrix(args):
    positions = set()
    contacts = defaultdict(lambda : Counter())
    with open(args.input, 'r') as f:
        for line in f:
            [_, pos1, _, pos2, value] = line.rstrip().split()
            pos1 = int(pos1) * args.resolution
            pos2 = int(pos2) * args.resolution
            positions.add(pos1)
            positions.add(pos2)
            contacts[pos1][pos2] = value
            contacts[pos2][pos1] = value
    positions = sorted(positions)
    return (contacts, positions)

def write_matrix_to_file(matrix, positions, args):
    with open(args.output, 'w') as f:
        for position in positions:
            f.write('\t' + args.assembly + '|' + args.chromosome + ':' +
                    str(position))
        f.write('\n')
        for pos1 in positions:
            f.write(args.assembly + '|' + args.chromosome + ':' + str(pos1))
            for pos2 in positions:
                f.write('\t' + str(matrix[pos1][pos2]))
            f.write('\n')

if __name__ == "__main__":
    main()
