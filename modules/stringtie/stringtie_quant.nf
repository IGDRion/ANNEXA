process QUANTIFY {
  conda (params.enable_conda ? "bioconda::stringtie" : null)
  container "${ workflow.containerEngine == 'Singularity' ? 
                'https://depot.galaxyproject.org/singularity/stringtie%3A2.2.3--h43eeafb_0' :
                'quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0'}"
  cpus params.maxCpu
  tag "$bam"

  input:
  tuple path(bam), val(av_length)
  path merged_gtf

  output:
  path "${bam.simpleName}.counts_transcripts.txt", emit: counts_transcript
  path "${bam.simpleName}.counts_gene.txt", emit: counts_gene

  script:
  """
  stringtie \
    -L \
    -eB \
    -G ${merged_gtf} \
    -o ${bam.simpleName}_quant.gtf \
    -p ${params.maxCpu} \
    ${bam}
  
  # Extract raw counts from stringtie with prepDE
  mkdir ${bam.simpleName}
  mv ${bam.simpleName}_quant.gtf ${bam.simpleName}

  prepDE.py \
  -i . \
  -g ${bam.simpleName}.counts_gene.txt \
  -t ${bam.simpleName}.counts_transcripts.txt \
  -l ${av_length}
  """
}