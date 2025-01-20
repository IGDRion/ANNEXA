process READ_LENGTH {
  conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/samtools%3A1.16.1--h6899075_0' : 
                'quay.io/biocontainers/samtools:1.16.1--h1170115_0' }"
  tag "$bam"

  input:
  file bam

  output:
  tuple path(bam), env(av_length), emit: bam_length
  env(av_length), emit: av_length
  
  script:
  """
  av_length=\$(samtools stats ${bam} | grep "average length" | awk '{print \$4}')
  """
}