process REPORT {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "mlorthiois/annexa:${workflow.revision? workflow.revision: "main"}"
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

  """
  qc_gtf.py -gtf ${gtf} \
    -c_gene ${counts_gene} \
    -ref ${ref} \
    -prefix ${prefix}
  qc.R ${prefix}
  """
}
