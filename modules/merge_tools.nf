process MERGE_TOOLS {
  conda (params.enable_conda ? "bioconda::gffcompare" : null)
  container "${ workflow.containerEngine == 'singularity' ?
      'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0' :
      'quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0' }"
  cpus params.maxCpu

  input:
  path bambu_gtf, stageAs: 'bambu.extended_annotations.gtf'
  path stringtie_gtf, stageAs: 'stringtie.extended_annotations.gtf'

  output:
  path("gffcmp.combined.gtf"), emit: merged_tools
  path("gffcmp.tracking"), emit: tracking
  
  script:
  """
  gffcompare ${bambu_gtf} ${stringtie_gtf}
  """
}
