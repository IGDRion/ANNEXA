process EXTRACT_QUANTS {
    conda (params.enable_conda ? "$baseDir/environment.yml" : null)
    container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
    publishDir "$params.outdir/stringtie2", mode: 'copy', pattern: '*.txt'
    if (params.filter == false){
        publishDir "$params.outdir/final", mode: 'copy', pattern: 'counts_transcript.txt', saveAs: {filename -> 'counts_transcript.full.txt'}
  }

    input:
    val bam_length
    path ctab
    path stringtie_merged

    output:
    path "counts_gene.txt", emit: gene_counts
    path "counts_transcript.txt", emit: tx_counts
    path "empty.ndr", emit: ndr
    path "rec_ndr.txt", emit: rec_ndr

    script:
    def mean_length = bam_length.collect { it as double }.sum() / bam_length.size()
    """
    extract_stringtie_counts.R \
    --av_length ${mean_length} \
    --ctabs ${ctab} \
    --ref ${stringtie_merged}

    # Create empty NDR file for TFKMERS to have same workflow as bambu in filtering step (placeholder)
    echo "1" > rec_ndr.txt
    """
}