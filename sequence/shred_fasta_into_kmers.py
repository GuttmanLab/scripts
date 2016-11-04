import argparse

def main():
    args = parse_arguments()
    shred_fastas_and_print(args)

class FastaRecord:
    def __init__(self, name, sequence):
        self.name = name.split()[0]  # Remove all characters after first whitespace
        self.sequence = sequence

    def shred(self, kmer_size, discard_too_short):
        kmers = []
        too_short = len(self.sequence) < kmer_size
        if too_short and not discard_too_short:
            kmers.append(FastaRecord(self.name + ":0-" + str(len(self.sequence)), self.sequence))
        elif not too_short:
            for start_coord in range(len(self.sequence) - kmer_size + 1):
                end_coord = start_coord + kmer_size
                subname = self.name + ":" + str(start_coord) + "-" + str(end_coord)
                subsequence = self.sequence[start_coord:end_coord]
                kmers.append(FastaRecord(subname, subsequence))
        return kmers

    def append(self, additional_sequence):
        self.sequence += additional_sequence

    def print(self):
        print(self.name)
        print(self.sequence)

def print_kmers(kmers):
    for kmer in kmers:
        kmer.print()

def shred_fastas_and_print(args):

    current_fasta = None

    with open(args.input, "r") as input_file:
        for line in input_file:
            if line.startswith(">"):
                if current_fasta:
                    kmers = current_fasta.shred(args.kmer_size, args.discard_too_short)
                    print_kmers(kmers)
                current_fasta = FastaRecord(line.rstrip(), "")
            else:
                current_fasta.append(line.rstrip())

        # shred the last sequence in the file
        if current_fasta:
            kmers = current_fasta.shred(args.kmer_size, args.discard_too_short)
            print_kmers(kmers)

def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "The FASTA file to shred",
                        action = "store")
    parser.add_argument("-k", "--size", help = "The size of the k-mers",
                        action = "store", type = int, default = 35, dest = "kmer_size")
    parser.add_argument("-d", "--discard", action = "store_true", dest = "discard_too_short",
                        help = "Discards FASTA records if length is less than k-mer size")
    return parser.parse_args()

if __name__ == "__main__":
    main()
