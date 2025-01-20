process EXTRACT_QUANTS {
    input:
    val bam_length
    path ctab
    path stringtie_merged

    output:
    path "counts_gene.txt", emit: gene_counts
    path "counts_transcript.txt", emit: tx_counts
    path "empty.ndr", emit: ndr

    script:
    def mean_length = bam_length.sum() / bam_length.size()
    """
    extract_stringtie_counts.R ${mean_length} ${ctab} ${stringtie_merged}

    # Create empty NDR file for TFKMERS to have same workflow as bambu in filtering step (placeholder)
    touch empty.ndr
    """
}