#! /usr/bin/env python3
from GTF import GTF


def parse_gene_counts(file):
    genes = {}
    samples = file.readline().rstrip().split("\t")
    for line in file:
        line = line.rstrip().split("\t")
        genes[line[0]] = line[1:]
    return genes, samples

def parse_transcripts_counts(file):
    transcripts = {}
    file.readline()
    for line in file:
        line = line.rstrip().split("\t")
        transcripts[line[0]] = line[2:]
    return transcripts

def get_ref_length(file):
    ref = {}
    for record in GTF.parse(file, by_line=True):
        if record.feature == "gene":
            ref[record["gene_id"]] = {"start": record.start, "end": record.end}
    return ref

def qc_gtf(gtf, gene_counts, transcript_counts, ref):
    ref_start_end = get_ref_length(ref)
    gene_counts, samples = parse_gene_counts(gene_counts)
    transcript_counts = parse_transcripts_counts(transcript_counts)
    samples = ",".join(samples)

    gene_str = (
        f"gene_id,gene_biotype,nb_transcripts,length,ext_5,ext_3,discovery,{samples}\n"
    )
    transcript_str = (
        f"gene_id,transcript_id,transcript_biotype,nb_exons,length,gene_discovery,discovery,{samples}\n"
    )
    exon_str = "exon_biotype,length,discovery\n"

    for gene in GTF.parse(gtf):
        if gene["gene_biotype"] not in ("protein_coding", "lncRNA"):
            continue

        gene_id = gene["gene_id"]
        ext = "NA"
        if not gene_id.startswith("gene."):
            if gene.strand == "+":
                ext_5 = ref_start_end[gene["gene_id"]]["start"] - gene.start
                ext_3 = gene.end - ref_start_end[gene["gene_id"]]["end"] 
            else:
                ext_3 = ref_start_end[gene["gene_id"]]["start"] - gene.start
                ext_5 = gene.end - ref_start_end[gene["gene_id"]]["end"] 

        length = max(
            [
                sum([len(exon) for exon in transcript.exons])
                for transcript in gene.transcripts
            ]
        )
        status = "novel" if gene_id.startswith("gene.") else "known"
        g_c = ",".join(gene_counts[gene_id])
        gene_str += f'{gene_id},{gene["gene_biotype"]},{len(gene.transcripts)},{length},{ext_5},{ext_3},{status},{g_c}\n'

        for transcript in gene.transcripts:
            tx_status = (
                "novel" if transcript["transcript_id"].startswith("tx.") else "known"
            )
            length = sum([len(exon) for exon in transcript.exons])
            t_c = ",".join(transcript_counts[transcript["transcript_id"]])
            transcript_str += f'{transcript["gene_id"]},{transcript["transcript_id"]},{transcript["transcript_biotype"]},{len(transcript.children)},{length},{status},{tx_status},{t_c}\n'

            for exon in transcript.children:
                exon_str += f'{exon["gene_biotype"]},{len(exon)},{status}\n'

    return gene_str, transcript_str, exon_str


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
    parser.add_argument(
        "-gtf",
        help="Path to your mRNA only GTF file.",
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "-c_gene",
        help="Path to bambu counts (genes) file.",
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "-c_tx",
        help="Path to bambu counts (transcripts) file.",
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "-ref",
        help="Path to reference annotation file.",
        type=argparse.FileType("r"),
    )
    args = parser.parse_args()

    gene, transcript, exon = qc_gtf(args.gtf, args.c_gene, args.c_tx, args.ref)
    with open("gene.stats", "w") as fd:
        fd.write(gene)
    with open("transcript.stats", "w") as fd:
        fd.write(transcript)
    with open("exon.stats", "w") as fd:
        fd.write(exon)
