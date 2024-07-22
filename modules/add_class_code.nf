process ADD_CLASS_CODE {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/final", mode: 'copy', overwrite: true
  tag "$gtf"

  input:
  file class_code_gtf
  file gtf

  output:
  file gtf

  script:
  """
  class_code.R ${class_code_gtf} ${gtf}
  ## Remove header created by gtfsort
  sed -i 1,3d ${gtf}
  """
}