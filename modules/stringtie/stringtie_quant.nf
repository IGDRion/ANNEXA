process QUANTIFY {
  conda (params.enable_conda ? "bioconda::stringtie" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A2.2.3--h43eeafb_0' :
                'quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0'}"
  cpus params.maxCpu
  tag "$bam"

  input:
  path(bam)
  path merged_gtf

  output:
  path "${bam.simpleName}.ctab", emit: ctab

  script:
  """
  stringtie \
    -L \
    -eB \
    -G ${merged_gtf} \
    -o ${bam.simpleName}_quant.gtf \
    -p ${params.maxCpu} \
    ${bam}
  
  mv t_data.ctab ${bam.simpleName}.ctab
  """
}