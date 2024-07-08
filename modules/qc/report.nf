process REPORT {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/qc/${prefix}", mode: 'copy'

  input:
  file gtf
  file ref
  file counts_gene
  val prefix

  output:
  file "${prefix}.gene.stats"
  file "${prefix}.transcript.stats"
  file "${prefix}.exon.stats"
  file "${prefix}.annexa.qc.pdf"
  path "missing_genes.txt", optional: true

  script:
  """
  qc_gtf.py -gtf ${gtf} \
    -c_gene ${counts_gene} \
    -ref ${ref} \
    -prefix ${prefix} \
    -tx_discovery ${params.tx_discovery}
  qc.R ${prefix} ${workflow.manifest.version}
  """
}
