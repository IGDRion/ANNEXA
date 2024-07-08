#! /usr/bin/env python3
# Fix exon numbers on negative strand in Stringtie output
import re

def reverse_exon_numbers(input_file, output_file):
    current_transcript = None
    current_strand = None
    exons = []
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            
            if fields[2] == 'transcript':
                # Write previous transcript's exons if any
                if current_transcript and exons:
                    write_exons(outfile, exons, current_strand)
                    exons = []
                
                # Write the transcript line as is
                outfile.write(line)
                current_transcript = fields[8].split('transcript_id')[1].split(';')[0].strip(' "')
                current_strand = fields[6]
            
            elif fields[2] == 'exon':
                exons.append(line)
            
            else:
                # Write any other type of line as is
                outfile.write(line)
        
        # Write the last transcript's exons
        if exons:
            write_exons(outfile, exons, current_strand)

def write_exons(outfile, exons, strand):
    total_exons = len(exons)
    if strand == '-':
        for i, exon in enumerate(exons, 1):
            new_exon = re.sub(r'exon_number "\d+"', f'exon_number "{total_exons - i + 1}"', exon)
            outfile.write(new_exon)
    else:
        for exon in exons:
            outfile.write(exon)

# Usage
input_file = 'novel.full.gtf'
output_file = 'fixed_novel.full.gtf'
reverse_exon_numbers(input_file, output_file)