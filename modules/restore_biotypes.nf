process RESTORE_BIOTYPE {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"

  input:
  file ref
  file novel_isoforms

  output:
  path "novel.isoforms.gtf"

  """
  restore_ref_attributes.py -gtf ${novel_isoforms} -ref $ref > novel.isoforms.gtf
  """
}
