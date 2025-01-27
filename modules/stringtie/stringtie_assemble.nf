process ASSEMBLE {
  conda (params.enable_conda ? "bioconda::stringtie=3.0.0" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A3.0.0--h29c0135_0' :
                'quay.io/biocontainers/stringtie:3.0.0--h29c0135_0'}"
  cpus params.maxCpu

  input:
  path bam
  path ref
  
  output:
  path "${bam.baseName}_stringtie_assemble.gtf", emit: stringtie_assemble_gtf

  script:
  """
  stringtie \
    -L \
    -p ${params.maxCpu} \
    -G ${ref} \
    -o ${bam.baseName}_stringtie_assemble.gtf \
    ${bam}
  """
}