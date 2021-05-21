#! /usr/bin/env python3
import argparse
from GTF import GTF

def restore_from_ref(ref, gtf):
    genes = {}
    transcripts = {}

    with open(ref) as fd:
        for record in GTF.parse(fd, by_line=True):
            if record.feature == "gene":
                genes[record["gene_id"]] = record.attributes
            elif record.feature == "transcript":
                transcripts[record["transcript_id"]] = record.attributes
            
    with open(gtf) as fd:
        for record in GTF.parse(fd, by_line=True): 
            if record["transcript_id"].startswith("tx."):
                record.attributes.update(genes[record["gene_id"]])
                record["transcript_biotype"] = record["gene_biotype"]
                if "gene_source" in record:
                    record["transcript_source"] = "bambu"
            else:
                record.source = "ensembl"
                record.attributes.update(transcripts[record["transcript_id"]])

            print(record)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Restore ref attributes in novel gtf.')
    parser.add_argument('-gtf', type=str, help='gtf file to prettify', required=True)
    parser.add_argument('-ref', type=str, help='ref annotation where informations are', required=True)
    args = parser.parse_args()

    restore_from_ref(args.ref, args.gtf)
