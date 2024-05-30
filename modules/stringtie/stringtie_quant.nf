process STRINGTIE_QUANTIFY {
  conda (params.enable_conda ? "bioconda::stringtie" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A2.2.3--h43eeafb_0' :
                'quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0'}"
  cpus params.maxCpu

  input:
  path '*'
  path merged_gtf

  output:
  path 'stringtie_quant.gtf', emit: stringtie_quant_qtf
  path 'e_data.ctab'
  path 'e2t.ctab'
  path 'i_data.ctab'
  path 'i2t.ctab'
  path 't_data.ctab'

  script:
  """
  stringtie \
    -L \
    -eB \
    -p ${params.maxCpu} \
    -G ${merged_gtf} \
    -o stringtie_quant.gtf \
    *.bam
  """
}