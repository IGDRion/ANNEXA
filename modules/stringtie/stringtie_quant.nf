process STRINGTIE_QUANTIFY {
  conda (params.enable_conda ? "bioconda::stringtie" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A2.2.3--h43eeafb_0' :
                'quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0'}"
  cpus params.maxCpu
  tag 'bam'

  input:
  path bam
  path merged_gtf

  output:
  path 'stringtie_quant.gtf', emit: stringtie_quant_qtf

  script:
  """
  stringtie \
    -L \
    -eB \
    -p ${params.maxCpu} \
    -G ${merged_gtf} \
    -o stringtie_quant.gtf \
    ${bam}
  """
}