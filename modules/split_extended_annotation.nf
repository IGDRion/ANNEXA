process SPLIT_EXTENDED_ANNOTATION {
  input:
  file extended_annotation

  output:
  path 'novel_genes.gtf', emit: novel_genes
  path 'novel_isoforms.gtf', emit: novel_isoforms
  path 'novel.gtf', emit: novel

  script:
  """
  grep -e "BambuTx" -e 'transcript_id "MSTRG' ${extended_annotation} | awk '\$3=="exon"' > novel.gtf
  
  # Check if novel.gtf is empty
  if [ ! -s novel.gtf ]; then
    echo "novel.gtf is empty. No novel transcripts detected, halting pipeline."
    exit 1
  fi

  grep -e "BambuGene" -e 'gene_id "MSTRG' -e "unstranded.Gene" novel.gtf > novel_genes.gtf
  grep -v -e "BambuGene" -e 'gene_id "MSTRG' -e "unstranded.Gene" novel.gtf > novel_isoforms.gtf
  """
}