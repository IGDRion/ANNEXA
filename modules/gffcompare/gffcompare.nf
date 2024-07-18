process GFFCOMPARE {
  conda (params.enable_conda ? "bioconda::gffcompare" : null)
  container "${ workflow.containerEngine == 'singularity' ?
      'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0' :
      'biocontainers/gffcompare:0.12.6--h9f5acd7_0' }"
  cpus params.maxCpu

    input:
    path reference_gtf
    path fasta
    path merged_gtf

    output:
    path("*.annotated.gtf"), emit: class_code_gtf
    path("gffcmp.tracking"), emit: tracking_file
  
    script:
    """
    gffcompare \
    -r ${reference_gtf} \
    -s ${fasta} \
    ${merged_gtf}
    """
}