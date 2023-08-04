process EXTRACT_TSS_SEQUENCES {
  conda (params.enable_conda ? "bioconda::bedtools=2.30" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' : 
                'quay.io/biocontainers/bedtools:2.30.0--h468198e_3' }"

  input:
  file tss_regions
  file fa

  output:
  path "tss.fa", emit: tss_sequences

  """
  bedtools getfasta -nameOnly -fi ${fa} -bed ${tss_regions} | tr a-z A-Z > tss.fa
  """
}
