process RESTRAND_ISOFORMS {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"

  input:
  file gtf

  output:
  path "restranded.${gtf}"

  shell:
  '''
  restrand_isoforms.R !{gtf} !{params.tx_discovery} "restranded.!{gtf}"
  sed -i '/;\s*$/!s/\s*$/;/' "restranded.!{gtf}"
  '''
}