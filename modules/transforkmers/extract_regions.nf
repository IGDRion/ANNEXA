process EXTRACT_TSS_REGIONS {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"

  input:
  file novel_gtf

  output:
  path "tss_regions.bed6", emit: tss_regions

  """
  cat ${novel_gtf} | extract_tss.py -l 512 > tss_regions.bed6
  """
}
