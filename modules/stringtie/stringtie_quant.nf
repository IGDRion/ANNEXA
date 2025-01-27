process QUANTIFY {
  conda (params.enable_conda ? "bioconda::stringtie=3.0.0" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A3.0.0--h29c0135_0' :
                'quay.io/biocontainers/stringtie:3.0.0--h29c0135_0'}"
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