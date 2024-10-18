process RESTRAND_NOVEL_ISOFORMS {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/final", mode: 'copy', saveAs: {filename -> "${gtf}"}, overwrite: true
  tag "$gtf"

  input:
  file gtf

  output:
  path "restranded.${gtf}"

  script:
  """
  restrand_isoforms.R ${gtf} ${params.tx_discovery} "restranded.${gtf}"
  """
}