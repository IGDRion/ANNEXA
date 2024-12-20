process MERGE {
  conda (params.enable_conda ? "bioconda::stringtie" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A2.2.3--h43eeafb_0' :
                'quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0'}"
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