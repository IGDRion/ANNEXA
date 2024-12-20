process ASSEMBLE {
  conda (params.enable_conda ? "bioconda::stringtie" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A2.2.3--h43eeafb_0' :
                'quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0'}"
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