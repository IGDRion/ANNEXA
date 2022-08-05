#! /usr/bin/env python3
from GTF import GTF


if __name__ == "__main__":
    #######################################################
    import argparse

    parser = argparse.ArgumentParser(description="Check GTF format.")
    parser.add_argument(
        "gtf",
        help="Path to gtf to check for ANNEXA.",
        type=argparse.FileType("r"),
    )
    args = parser.parse_args()
    #######################################################

    for record in GTF.parse_by_line(args.gtf):
        if record.feature == "gene" or record.feature == "transcript":
            continue

        #######################################################
        # Check each record has gene_id and transcript_id attributes
        if not "gene_id" in record:
            exit(f"No gene_id feature in this line:\n{str(record)}")
        if not "transcript_id" in record:
            exit(f"No gene_id feature in this line:\n{str(record)}")

        #######################################################
        # Replace if gene_id starts with "gene." to avoid conflicts with bambu
        g_id = record["gene_id"]
        if g_id.startswith("gene."):
            record["gene_id"] = g_id.replace("gene", "known_gene")

        # Replace if transcript_id starts with "tx." to avoid conflicts with bambu
        t_id = record["transcript_id"]
        if t_id.startswith("tx."):
            record["transcript_id"] = t_id.replace("tx", "known_tx")

        #######################################################
        # Check if gene_biotype in each transcripts and exons
        if not "gene_biotype" in record:
            record["gene_biotype"] = "NA"
        g_biotype = record["gene_biotype"]

        # Check for RefSeq gene_biotype format
        if g_biotype == "mRNA":
            g_biotype = "protein_coding"
        if g_biotype == "lnc_RNA":
            g_biotype = "lncRNA"
        record["gene_biotype"] = g_biotype

        #######################################################
        if not "transcript_biotype" in record:
            record["transcript_biotype"] = g_biotype
        t_biotype = record["transcript_biotype"]

        # Check for RefSeq transcript_biotype format
        if t_biotype == "mRNA":
            t_biotype = "protein_coding"
        elif t_biotype == "lnc_RNA":
            t_biotype = "lncRNA"
        record["transcript_biotype"] = t_biotype

        print(str(record))
