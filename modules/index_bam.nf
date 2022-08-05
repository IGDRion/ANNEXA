process INDEX_BAM {
  conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' : 
                'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

  input:
  file bam

  output:
  file("*.bai")

  """
  samtools index $bam
  """
}
