process MERGE_ANNOTATIONS {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file novel
  file ref
  val origin

  output:
  path "${params.prefilter_ndr ? 'prefilter.' : ''}extended_annotations.${origin}.gtf", emit: gtf

  script:
  def output_prefix = params.prefilter_ndr ? 'prefilter.' : ''
  """
  cat ${novel} ${ref} | GTF.py format > ${output_prefix}extended_annotations.${origin}.gtf
  """
}
