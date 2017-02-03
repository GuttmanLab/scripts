import argparse
import statistics
import subprocess

def main():
    args = parse_arguments()
    biases = calculate_bias_factors(args)
    ice(biases, args)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Runs iterative correction on a matrix text file.')

    parser.add_argument('-i', '--input',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The input matrix file.')

    parser.add_argument('-o', '--output',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The output ICEd-matrix file.')

    parser.add_argument('--ice',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'Location of ICE executable.')

    parser.add_argument('--bias',
                        metavar = 'FILE',
                        action = 'store',
                        help = 'The intermediate file of bias factors.')

    parser.add_argument('--iterations',
                        metavar = 'INT',
                        default = "100",
                        help = 'Number of ICE iterations.')

    return parser.parse_args()

def calculate_bias_factors(args):
    has_column_header = "1"
    has_row_header = "1"
    num_lines = count_lines_in_file(args.input, ignoreheader = True)
    subprocess.check_call([args.ice, args.input, str(num_lines),
            args.iterations, has_row_header, has_column_header, args.bias])
    return parse_bias_file(args.bias)

def count_lines_in_file(fname, ignoreheader = False):
    count = 0
    with open(fname, 'r') as f:
        if ignoreheader:
            next(f)
        for line in f:
            count += 1
    return count

def parse_bias_file(fname):
    biases = []
    with open(fname) as f:
        for line in f:
            biases.append(float(line.strip()))
    return biases

def ice(biases, args):
    med_diag_val = get_median_diagonal_value(args)
    with open(args.input, 'r') as f, open(args.output, 'w') as outf:
        outf.write(next(f))
        row = 0
        for line in f:
            outf.write(line.split()[0])
            values = line.split()[1:]
            for col in range(0, len(values)):
                val = float(values[col])
                if val > 0:
                    val /= (biases[row] * biases[col])
                    val = 1 if val >= med_diag_val else val / med_diag_val
                outf.write("\t" + str(val))
            outf.write("\n")
            row += 1    
            
def get_median_diagonal_value(args):
    diagonal_values = []
    with open(args.input, 'r') as f:
        next(f)
        row = 0
        for line in f:
            values = line.split()[1:]
            for col in range(0, len(values)):
                if abs(row - col) == 1:
                    diagonal_values.append(float(values[col]))
            row += 1
    return statistics.median(diagonal_values)

if __name__ == "__main__":
    main()
