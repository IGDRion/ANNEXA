process MERGE {
  conda (params.enable_conda ? "bioconda::stringtie=3.0.0" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A3.0.0--h29c0135_0' :
                'quay.io/biocontainers/stringtie:3.0.0--h29c0135_0'}"
  cpus params.maxCpu

  input:
  path stringtie_gtf
  file ref

  output:
  path 'stringtie_merged.gtf', emit: stringtie_merged_gtf

  script:
  """
  stringtie \
    -L \
    -p ${params.maxCpu} \
    --merge ${stringtie_gtf} \
    -G ${ref} \
    -o stringtie_merged.gtf
  """
}