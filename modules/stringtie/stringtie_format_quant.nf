process FORMAT_QUANT {
  tag "$gene_count"

  input:
  path gene_count
  path transcript_count
  path merged_gtf

  output:
  path "reformated_output/${gene_count}", emit: counts_gene
  path "reformated_output/${transcript_count}", emit: counts_transcript

  script:
  """
  mkdir reformated_output

  reformat_stringtie_transcripts.sh ${transcript_count} ${merged_gtf} reformated_output/${transcript_count}
  sed -i 's/\\r\$//' reformated_output/${transcript_count}

  sed 's/\\([^,|]*\\)|[^,]*/\\1/' ${gene_count} > reformated_output/${gene_count}
  sed -i 's/,/\\t/' reformated_output/${gene_count}
  sed -i 's/[[:space:]]\\+/\\t/g' reformated_output/${gene_count}
  """
}