#! /usr/bin/env python3
from GTF import GTF


if __name__ == "__main__":
    # Args
    import argparse

    parser = argparse.ArgumentParser(
        description="Utility tool to format your GTF file for ANNEXA."
    )
    parser.add_argument(
        "gtf",
        help="Path for gtf file you want to format for ANNEXA.",
        type=argparse.FileType("r"),
    )
    args = parser.parse_args()

    gtf = GTF.parse(args.gtf)

    for gene in gtf.values():
        if not "gene_id" in gene:
            exit(f"No gene_id feature in this line:\n{str(gene)}")

        # Add gene_id in each transcripts and exons
        g_id = gene["gene_id"]

        # Replace if gene_id starts with "gene." to avoid conflicts with bambu
        if g_id.startswith("gene."):
            gene["gene_id"] = g_id.replace("gene", "KNOWN_GENE")

        for transcript in gene.transcripts:
            if "gene_id" not in transcript:
                transcript["gene_id"] = g_id
        for exon in gene.exons:
            if "gene_id" not in exon:
                exon["gene_id"] = g_id

        # Check if gene_biotype in each transcripts and exons
        if not "gene_biotype" in gene:
            exit(f"No gene_biotype feature in this line:\n{str(gene)}")
        g_biotype = gene["gene_biotype"]

        # Check for RefSeq gene_biotype format
        if g_biotype == "mRNA":
            gene["gene_biotype"] = "protein_coding"
        elif g_biotype == "lnc_RNA":
            gene["gene_biotype"] = "lncRNA"

        for transcript in gene.transcripts:
            if "gene_biotype" not in transcript:
                transcript["gene_biotype"] = g_biotype
        for exon in gene.exons:
            if "gene_biotype" not in transcript:
                exon["gene_biotype"] = g_biotype

        # Remove transcript_id in gene
        if "transcript_id" in gene:
            del gene["transcript_id"]

        # Check each transcripts
        for transcript in gene.transcripts:
            # Check for transcript_id in tx and exons
            if not "transcript_id" in transcript:
                exit(f"No transcript_id feature in this line:\n{str(transcript)}")
            tx_id = transcript["transcript_id"]

            # Check if transcript_id doesn't starts with "tx."
            if tx_id.startswith("tx."):
                transcript["transcript_id"] = tx_id.replace("tx.", "KNOWN_TX.")

            for exon in transcript.exons:
                if "transcript_id" not in exon:
                    exon["transcript_id"] = tx_id

            # Check for transcript_biotype in tx and exons
            if not "transcript_biotype" in transcript:
                # If no transcript_biotype, use gene_biotype
                tx_biotype = transcript["gene_biotype"]
            else:
                tx_biotype = transcript["transcript_biotype"]

            # Check for RefSeq transcript_biotype format
            if tx_biotype == "mRNA":
                transcript["transcript_biotype"] = "protein_coding"
            elif tx_biotype == "lnc_RNA":
                transcript["transcript_biotype"] = "lncRNA"

            for exon in transcript.exons:
                if (
                    "transcript_biotype" not in exon
                    or exon["transcript_biotype"] != transcript["transcript_biotype"]
                ):
                    exon["transcript_biotype"] = transcript["transcript_biotype"]

        print(gene.format_to_gtf())
