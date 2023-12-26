process GENEPRED_TO_BED {
  conda (params.enable_conda ? "bioconda::ucsc-genepredtobed=377" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/ucsc-genepredtobed:377--ha8a8165_5' : 
                'quay.io/biocontainers/ucsc-genepredtobed:377--ha8a8165_5' }"

  input:
  file genePred

  output:
  path "${genePred.simpleName}.bed12"

  """
  genePredToBed ${genePred} ${genePred.simpleName}.bed12
  """
}
