#! /usr/bin/env python3
from GTF import GTF
import argparse

def main(ref, merged):
    ref_genes = {}
    
    for gene in GTF.parse(ref):
        ref_genes[gene["gene_id"]] = (gene.start, gene.end)


    for gene in GTF.parse(merged):
        ext_5 = ""
        ext_3 = ""
        if gene.start < ref_genes[gene["gene_id"]][0]:
            ext_5 = ",5'_ext"
        if gene.end > ref_genes[gene["gene_id"]][1]:
            ext_3 = ",3'_ext"

        if ext_3 != "" or ext_5 != "":
            print(f"{gene['gene_id']}{ext_5}{ext_3}")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
    parser.add_argument(
        "-ref", help="Path to reference annotation", type=argparse.FileType("r"), required=True
    )
    parser.add_argument(
        "-merged", help="Path to your novel merged file", type=argparse.FileType("r"), required=True
    )

    args = parser.parse_args()
    main(args.ref, args.merged)
