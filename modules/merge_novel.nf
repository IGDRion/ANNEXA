process MERGE_NOVEL {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file novel_genes
  file novel_isoforms

  output:
  path "novel.full.gtf"

  """
  cat ${novel_genes} ${novel_isoforms} | GTF.py format > novel.full.gtf
  """
}
