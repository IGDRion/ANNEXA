#! /usr/bin/env python3 
from GTF import GTF
import warnings

def parse_gene_counts(file):
    """
    Parse gene_counts from bambu in a dict, ex for line:
    ENSG0000003    25   16   3   0
    { "ENSG00000003": {
        "counts": 44,
        "validates" : 3
      }
    }
    """
    genes = {}
    file.readline()  # Skip header

    for line in file:
        line = line.rstrip().split("\t")
        int_line = [float(count) for count in line[1:]]

        # Total expression over samples
        sum_counts = sum(int_line)

        # Number of samples which has counts != 0
        validates = len([sample for sample in int_line if sample != 0])

        # Store in dic both value
        genes[line[0]] = {"validates": validates, "counts": sum_counts}
    return genes


def get_ref_length(file):
    ref = {}
    for gene in GTF.parse(file):
        ref[gene["gene_id"]] = {"start": gene.start, "end": gene.end}
    return ref


def qc_gtf(gtf, gene_counts, ref, tx_discovery):
    ref_start_end = get_ref_length(ref)
    gene_counts = parse_gene_counts(gene_counts)
    missing_genes = []

    if tx_discovery == "bambu":
        gene_prefix,tx_prefix = "BambuGene","BambuTx"
    else:
        gene_prefix,tx_prefix = "MSTRG.","MSTRG."

    # CSV headers
    gene_str = f"gene_id,gene_biotype,nb_transcripts,length,ext_5,ext_3,discovery,validate_by,presents_in_sample\n"
    transcript_str = (
        f"transcript_id,gene_id,transcript_biotype,nb_exons,length,gene_discovery,tx_discovery\n"
    )
    exon_str = "exon_biotype,length,discovery\n"

    biotypes = set(("protein_coding", "lncRNA", "lnc_RNA", "ncRNA"))

    for gene in GTF.parse(gtf).values():
        g_id = gene["gene_id"]
        # Verify that gene is not missing in counts_gene.txt
        # Some genes (probably duplicates) can be omitted by Stringtie
        if g_id not in gene_counts:
            missing_genes.append(gene.get_attributes())
            warnings.warn("Missing gene: "+g_id, stacklevel=3)
            continue

        g_biotype = gene["gene_biotype"]
        if gene["gene_biotype"] not in biotypes:
            continue

        g_status = "novel" if g_id.startswith((gene_prefix,'unstranded.Gene')) else "known"
        g_count = gene_counts[g_id]["counts"]  # Counts in all samples
        g_samples = gene_counts[g_id]["validates"]  # Found in x samples
        g_nb_tx = len(gene.transcripts)  # Number of isoforms

        # Compute genomic extension with start/end in ref
        ext_5 = 0
        ext_3 = 0
        if not g_id.startswith((gene_prefix,'unstranded.Gene')):
            if gene.strand == "+":
                ext_5 = ref_start_end[gene["gene_id"]]["start"] - gene.start
                ext_3 = gene.end - ref_start_end[gene["gene_id"]]["end"]
            else:
                ext_3 = ref_start_end[gene["gene_id"]]["start"] - gene.start
                ext_5 = gene.end - ref_start_end[gene["gene_id"]]["end"]

        # Gene length = longest transcript e.g with longest sum of exon length
        length = max(
            [sum([len(exon) for exon in transcript.exons]) for transcript in gene.transcripts]
        )

        # Create a csv row
        gene_str += f"{g_id},{g_biotype},{g_nb_tx},{length},{ext_5},{ext_3},{g_status},{g_count},{g_samples}\n"

        for transcript in gene.transcripts:
            tx_id = transcript["transcript_id"]
            tx_biotype = transcript["transcript_biotype"]
            tx_status = "novel" if tx_id.startswith(tx_prefix) else "known"
            tx_nb_exons = len(transcript.exons)
            tx_length = sum([len(exon) for exon in transcript.exons])

            # Create csv row
            transcript_str += (
                f"{tx_id},{g_id},{tx_biotype},{tx_nb_exons},{tx_length},{g_status},{tx_status}\n"
            )

            for exon in transcript.exons:
                ex_length = len(exon)
                exon_str += f"{tx_biotype},{ex_length},{tx_status}\n"

    # Record missing genes
    if len(missing_genes) > 0:
        warnings.warn(str(len(missing_genes))+" gene(s) skipped. Please see missing_genes.txt", stacklevel=3)
        with open('missing_genes.txt', 'w') as f:
            f.write("# Some genes were missing in the final output. These are likely to be duplicated" + "\n" + \
                    "# genes in the reference annotations and were omitted by Stringtie."+ "\n" \
                    + "# These genes are recorded here." + "\n")
            for gene in missing_genes:
                f.write(str(gene)+'\n')

    # Each csv is stored in string, return 3 strings as tuple
    return gene_str, transcript_str, exon_str


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
    parser.add_argument(
        "-gtf",
        help="Path to your mRNA only GTF file.",
        type=argparse.FileType("r"),
        required=True,
    )
    parser.add_argument(
        "-c_gene",
        help="Path to bambu or stringtie2 counts (genes) file.",
        type=argparse.FileType("r"),
        required=True,
    )
    parser.add_argument(
        "-ref",
        help="Path to reference annotation file.",
        type=argparse.FileType("r"),
        required=True,
    )
    parser.add_argument(
        "-prefix",
        help="Output prefix.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-tx_discovery",
        help="Quantification method. Choices: bambu, stringtie2",
        type=str,
        choices=["bambu", "stringtie2"],
        required=True,
    )
    args = parser.parse_args()

    gene, transcript, exon = qc_gtf(args.gtf, args.c_gene, args.ref, args.tx_discovery)
    with open(f"{args.prefix}.gene.stats", "w") as fd:
        fd.write(gene)
    with open(f"{args.prefix}.transcript.stats", "w") as fd:
        fd.write(transcript)
    with open(f"{args.prefix}.exon.stats", "w") as fd:
        fd.write(exon)
