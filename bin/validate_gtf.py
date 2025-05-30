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
        if "gene_biotype" in record:
            g_biotype = record["gene_biotype"]
        elif "gene_type" in record:
            g_biotype = record["gene_type"]
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
        # g_id = record["gene_id"]
        # if g_id.startswith("gene."):
        #     record["gene_id"] = g_id.replace("gene", "gene_known")
        #
        # Replace if transcript_id starts with "tx." to avoid conflicts with bambu
        # t_id = record["transcript_id"]
        # if t_id.startswith("tx."):
        #     record["transcript_id"] = t_id.replace("tx", "transcript_known")

        #######################################################
        # Check if gene_biotype in each transcripts and exons
        if not "gene_biotype" in record and "gene_type" not in record:
            record["gene_biotype"] = g_biotype

        # Check for RefSeq gene_biotype format
        if g_biotype == "mRNA":
            g_biotype = "protein_coding"
        if g_biotype == "lnc_RNA" or g_biotype == "ncRNA":
            g_biotype = "lncRNA"
        record["gene_biotype"] = g_biotype

        #######################################################
        if "transcript_biotype" in record:
            t_biotype = record["transcript_biotype"]
        elif "transcript_type" in record:
            record.attributes["transcript_biotype"] = record.attributes["transcript_type"]
            del record.attributes["transcript_type"]
            t_biotype = record["transcript_biotype"]
        else:
            record["transcript_biotype"] = g_biotype
            t_biotype = record["transcript_biotype"]

        # Check for RefSeq transcript_biotype format
        if t_biotype == "mRNA":
            t_biotype = "protein_coding"
        elif t_biotype == "lnc_RNA" or t_biotype == "ncRNA":
            t_biotype = "lncRNA"
        record["transcript_biotype"] = t_biotype

        print(str(record))
