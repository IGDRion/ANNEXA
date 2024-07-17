process ADD_CLASS_CODE {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/final", mode: 'copy', overwrite: true

  input:
  file class_code_gtf
  file extended_annotation

  output:
  path extended_annotation, emit: extended_annotation_class_code

  script:
  """
  class_code.R ${class_code_gtf} ${extended_annotation}
  # Remove header created by gtfsort
  sed -i 1,3d ${extended_annotation}
  """
}