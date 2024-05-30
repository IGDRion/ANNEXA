process SPLIT_EXTENDED_ANNOTATION {
  input:
  file extended_annotation

  output:
  path 'novel_genes.gtf', emit: novel_genes
  path 'novel_isoforms.gtf', emit: novel_isoforms

  shell:
  if (params.tx_discovery == "bambu")
    '''
    grep "BambuTx" !{extended_annotation} | awk '$3=="exon"' > novel.gtf
    grep -e "BambuGene" -e "unstranded.Gene" novel.gtf > novel_genes.gtf
    grep -v -e "BambuGene" -e "unstranded.Gene" novel.gtf > novel_isoforms.gtf
    '''
  else if (params.tx_discovery == "stringtie2")
    '''
    grep 'transcript_id "MSTRG' !{extended_annotation} | awk '$3=="exon"' > novel.gtf
    grep -e 'gene_id "MSTRG' -e "unstranded.Gene" novel.gtf > novel_genes.gtf
    grep -v -e 'gene_id "MSTRG' -e "unstranded.Gene" novel.gtf > novel_isoforms.gtf
    '''
}