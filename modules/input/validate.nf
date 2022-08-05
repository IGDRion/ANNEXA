process VALIDATE_INPUT_GTF {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"

  input:
  file gtf

  output:
  file 'input.formatted.gtf'

  """
  validate_gtf.py ${gtf} > input.formatted.gtf
  """
}
