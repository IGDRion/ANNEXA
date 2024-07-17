process GTFSORT {
    conda (params.enable_conda ? "bioconda::gtfsort=0.2.2" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/gtfsort:0.2.2--h4ac6f70_0':
        'biocontainers/gtfsort:0.2.2--h4ac6f70_0' }"
    publishDir "$params.outdir/transdecoder", mode: 'copy', pattern: 'novel_tx_CDS.gtf'
    //publishDir "$params.outdir/${params.tx_discovery}", mode: 'copy', saveAs: 'extended_annotations.gtf', overwrite: true
    cpus params.maxCpu
    memory params.maxMemory

    input:
    path fixed_novel
    path exon_cds

    output:
    path "novel_tx_CDS.gtf", emit: gtf

    script:
    """
    #Add CDS to novel gtf only (final output in /results/transdecoder)
    
    # Merge
    cat ${fixed_novel} ${exon_cds} > merged.gtf

    gtfsort \
        -i merged.gtf \
        -o novel_tx_CDS.gtf \
        -t ${params.maxCpu}
    """
}