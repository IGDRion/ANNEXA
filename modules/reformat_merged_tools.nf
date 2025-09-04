process REFORMAT_MERGED_TOOLS {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"

  input:
  file merged_tools
  file tracking

  output:
  path "reformated.gtf", emit: reformated_merged_gtf

  script:
  """
  reformat_merge.py
  """
}