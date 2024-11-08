process GTFSORT {
    conda (params.enable_conda ? "bioconda::gtfsort=0.2.2" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/gtfsort:0.2.2--h4ac6f70_0':
        'biocontainers/gtfsort:0.2.2--h4ac6f70_0' }"
    publishDir "$params.outdir/transdecoder", mode: 'copy', pattern: 'novel.full.gtf'
    publishDir "$params.outdir/final", mode: 'copy', pattern: 'novel.full.gtf', saveAs: {filename -> 'novel.full.gtf'}, overwrite: true

    input:
    path fixed_novel
    path exon_cds

    output:
    path "novel.full.gtf", emit: gtf

    script:
    """    
    cat ${fixed_novel} ${exon_cds} > merged.gtf

    gtfsort \
        -i merged.gtf \
        -o novel.full.gtf \
        -t ${params.maxCpu}
    """
}