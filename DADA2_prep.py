#!bin/bash
"""
Description:
This is a python based wrapper for adapter trimming with cutadapt.
After running this, run the DADA2 pipeline using the DADA2_notebook.Rmd in this directory.

Docs @ github.com/jimbopants/DADA2DADA2_processing
Get help by running `python DADA2_prep.py --h`

Argparse Options:
  -h, --help            show this help message and exit
  --primer_set {amoA,nxrB,16S_515F_926R}
                        Predefined primer sets
  --raw_dir RAW_DIR     Directory with raw fastq reads
  --out_dir OUT_DIR     Output directory
  --fwd FWD             New forward primer (Ignored if using the primer_set
                        option)
  --rev REV             New reverse primer (Ignored if using the primer_set
                        option)
  --verify_only         Verifies the first sample contains the specified
                        primers, prints cutadapt output and exits.
  --names               Prints the trimmed names for each sample. Does not
                        actually trim reads
                        
Quick note on cutadapt:
Paired end reads are written with the 2nd read not reverse complemented, so use original primer and check 5' end of both reads
cutadapt -g trims 5' adapter.
I don't use cutadapt's paired end option because I remember it causing an issue last time
Instead I just do each file individually and discard pairs with missing read IDs in DADA2.

author: JG
email: jimbogrif@gmail.com
Created: 1/8/18
Updated: 10/9/18
"""
## Imports:
import subprocess
import glob
import argparse
import sys
import os
import shutil

# Check dependencies and print header:
print('\nDADA2 Cutadapt Preparation Script')
print('v1.0 Written by Jim Griffin, 1/8/18')
print('See docs at github.com/jimbopants/DADA2_processing\n')
if shutil.which('cutadapt') is None:
    print('Cutadapt not in PATH, check environment')
    sys.exit(0)

def main():
    args = parse_arguments()
    primers = check_and_set_primers(args)
    os.makedirs(args.out_dir, exist_ok=True)

    # Get reads, checks that there are an even number of files (1 fwd/rev per sample) otherwise exits.
    files = []
    for root, dirs, dir_files in os.walk(args.raw_dir, topdown=False):
        for name in dir_files:
            if '.DS' in name:
                pass
            else:
                files.append(os.path.join(root, name))
    if len(files)%2 != 0:
        print('Uneven number of files. Check raw reads')
        sys.exit()

    # Cutadapt trim reads
    name_prefix = ['F_', 'R_']
    read_index = 0
    for i in files:
        cutadapt_process(i, args.out_dir, primers, name_prefix, read_index)
        read_index +=1
        read_index = read_index %2
        # Break early if just running the first 2 samples for verification:
        if args.verify_only:
            if read_index == 0:
                break
    if args.names:
        for i in files:
            print(split_name(i))

## Subfunctions:
def split_name(file_str):
    name = file_str.rsplit('/',1)[1]
    name = name.split('_', 1)[0]
    return name

def cutadapt_process(file, out_dir, primers, name_prefix, read_index):
    """
    Uses subprocess to trim 5' adapters from the reads in file,
    writes only the reads where the adapter was trimmed.
    """
    name = split_name(file)
    output_file = out_dir + name_prefix[read_index] + name + '.fastq'
    trim_command = 'cutadapt -g {0} -o {1} --discard-untrimmed {2}'.format(primers[read_index], output_file, file)
    subprocess.check_call(trim_command, shell=True)
    return

def parse_arguments():
    parser = argparse.ArgumentParser()
    # Default options:
    parser.add_argument("--primer_set", help="Predefined primer sets", choices=['amoA', 'nxrB', '16S_515F_926R'])
    parser.add_argument('--raw_dir', help='Directory with raw fastq reads', required=True)
    parser.add_argument('--out_dir', help='Output directory')
    # Allow manual primer input:
    parser.add_argument('--fwd', help='New forward primer (Ignored if using the primer_set option)')
    parser.add_argument('--rev', help='New reverse primer (Ignored if using the primer_set option)')
    ## Simple run-tests
    # Option to check just the first sample for the correct primer
    parser.add_argument('--verify_only',
    help='Verifies the first sample contains the specified primers, prints cutadapt output and exits.',
    action='store_true')
    # Option to print the result of the name prefix trimming
    parser.add_argument('--names',
    help='Prints the trimmed names for each sample. Does not actually trim reads',
    action='store_true')
    # Print help if no options given.
    if len(sys.argv)==1:
        parser.print_help()
        print("\n\nNeed command line input\n\n")
        sys.exit(1)
    # Parse Command Line Arguments:
    try:
        result = parser.parse_args()
        return result
    except Exception as e:
        parser.print_help()
        print(e)
        sys.exit(0)

def check_and_set_primers(args):
    """
    After parsing arguments, this function checks that a primer set was specified
    and returns the actual primer sets if a default primerset was specified at the command line.
    """
    if args.primer_set is None:
        if all([args.fwd, args.rev]):
            return args.fwd, args.rev
        else:
            print("\n\nNeed to enter a fwd/rev primer or a valid primer set name. Try --help at command line\n\n")
            sys.exit(1)
    primer_dict = {
    'amoA' : ['GGGGTTTCTACTGGTGGT', 'CCCCTCKGSAAAGCCTTCTTC'],
    'nxrB' : ['TACATGTGGTGGAACA', 'CGGTTCTGGTCRATCA'],
    '16S_515F_926R' : ['GTGYCAGCMGCCGCGGTAA', 'CCGYCAATTYMTTTRAGTTT'] # 2016+ version.
    }
    return primer_dict[args.primer_set]

if __name__ == "__main__":
    main()
