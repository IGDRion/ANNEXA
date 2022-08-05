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
  path "extended_annotations.${origin}.gtf"

  """
  cat ${novel} ${ref} | GTF.py format > extended_annotations.${origin}.gtf
  """
}
