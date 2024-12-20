#! /usr/bin/env python3
#Adds missing exon number to the CDS in the final output of Transdecoder
import argparse
import re

parser = argparse.ArgumentParser(description='Add exon number to predicted CDS from TransDecoder')
parser.add_argument('--gff', type=str, required=True,
                    help='Path to final TransDecoder genome.gff3 output')
args = parser.parse_args()

def find_exons(GFF):
    with open(GFF,'r') as input:
        with open('exon_cds.gff3', 'w') as output:
            previous_line = input.readline()
            output.write(previous_line)
            for line in input:
                if "\tCDS\t" in line:
                    match = re.search(r'exon(\d+)', previous_line)
                    exon_number = match.group(1)

                    string = f";exon_number={exon_number}"
                    output.write(line.rstrip()+string+"\n")

                    previous_line = line
                else:
                    output.write(line)
                    previous_line = line

if __name__ == "__main__":
    find_exons(args.gff)