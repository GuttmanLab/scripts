from collections import defaultdict
import argparse
import sys

def main():
    args = parse_arguments()
    chromosome_sizes = load_size_file(args.sizes)
    print_tiling(chromosome_sizes, args.resolution)

def parse_arguments():
    parser = ParserWithDefaultError(description = 
        "Generates an SAF file for a no-overlap tiling of a genome. " +
        "Prints to STDOUT")

    parser.add_argument("-r", "--resolution", metavar = "NUM",
        action = "store", type = int, default = 1000000,
        help = "Tile size. (default = 1 Mb)")

    parser.add_argument("-s", "--sizes", metavar = "FILE",
        action = "store",
        help = "Size file. See /storage/Genomes/mm9/sizes for an example")

    return parser.parse_args()

class ParserWithDefaultError(argparse.ArgumentParser):
    def error(self, message):
        print("error: " + message, file = sys.stderr)
        self.print_help()
        sys.exit(2)

def load_size_file(size_path):
    chromosome_sizes = defaultdict(int)
    with open(size_path, "r") as f:
        for line in f:
            chromosome, size = line.split('\t')
            chromosome_sizes[chromosome] = int(size)
    return chromosome_sizes

# GFF file intervals are 1-based, closed.
# Since featureCounts only uses GFF and its own SAF format, assume the interval
# conventions are the same.
def print_tiling(chromosome_sizes, resolution):
    print_saf_header()
    default_strand = "+"
    for chromosome, size in chromosome_sizes.items():
        num_iterations = size // resolution + 1
        for i in range(1, size, resolution):
            start = i
            end = i + resolution - 1
            geneid = chromosome + ":" + str(start) + "-" + str(end)
            print("\t".join((geneid, chromosome, str(start), str(end), default_strand)))

def print_saf_header():
    print("\t".join(("GeneID", "Chr", "Start", "End", "Strand")))

if __name__ == "__main__":
    main()
