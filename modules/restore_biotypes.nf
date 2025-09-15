process RESTORE_BIOTYPE {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"

  input:
  file ref
  file novel_isoforms

  output:
  path "restored.gtf"

  script:
  """
  grep -e "BambuGene" -e 'gene_id "MSTRG' -e "unstranded.Gene" ${novel_isoforms} > novel_genes.gtf
  grep -v -e "BambuGene" -e 'gene_id "MSTRG' -e "unstranded.Gene" ${novel_isoforms} > novel_isoforms.gtf
  restore_ref_attributes.py -gtf novel_isoforms.gtf -ref ${ref} > restored.isoforms.gtf

  cat novel_genes.gtf restored.isoforms.gtf > restored.gtf
  """
}