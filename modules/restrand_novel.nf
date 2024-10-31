process RESTRAND_NOVEL {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/final", mode: 'copy', saveAs: {filename -> 'novel.full.gtf'}, overwrite: true

  input:
  file gtf

  output:
  path "${gtf}"

  shell:
  '''
  restrand_isoforms.R !{gtf} !{params.tx_discovery} "restranded.!{gtf}"
  mv restranded.!{gtf} !{gtf}
  sed -i '/;\s*$/!s/\s*$/;/' !{gtf}
  '''
}