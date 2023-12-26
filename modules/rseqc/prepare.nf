process PREPARE_RSEQC {
  input:
  file novel
  file known

  output:
  path "*.rseqc.gtf"

  """
  grep "protein_coding" ${novel} > novel_mRNA.rseqc.gtf
  grep "lncRNA" ${novel} > novel_lncRNA.rseqc.gtf
  grep "protein_coding" ${known} > known_mRNA.rseqc.gtf
  grep "lncRNA" ${known} > known_lncRNA.rseqc.gtf
  """
}
