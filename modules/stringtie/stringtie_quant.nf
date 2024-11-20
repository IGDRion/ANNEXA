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
  path "${bam.simpleName}.counts_transcript.txt", emit: counts_transcript
  path "${bam.simpleName}.counts_gene.txt", emit: counts_gene

  script:
  """
  stringtie \
    -L \
    -eB \
    -G ${merged_gtf} \
    -o stringtie_quant.gtf \
    -A gene_abundance.tab \
    -p ${params.maxCpu} \
    ${bam}

    # Extract tx and gene counts
    awk -v bam="${bam.simpleName}" 'BEGIN {OFS="\t"; print "transcript_id", "gene_id", bam}  NR>1 {print \$6, \$9, \$11}' t_data.ctab > ${bam.simpleName}.counts_transcript.txt
    awk -v bam="${bam.simpleName}" 'BEGIN {OFS="\t"; print "gene_id", bam}  NR>1 {print \$1, \$7}' gene_abundance.tab > ${bam.simpleName}.counts_gene.txt
  """
}