process BAMBU {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/bambu", mode: 'copy'
  cpus params.maxCpu
  memory params.maxMemory

  input:
  // Collect parallel bam channel into 1
  file '*' 
  file ref
  file fa

  output:
  path 'extended_annotations.gtf', emit: bambu_gtf
  path 'counts_transcript.txt', emit: tx_counts
  path 'counts_gene.txt', emit: gene_counts
  path 'bambu_ndr.csv', emit: ndr
  path 'rec_ndr.txt', emit: rec_ndr
  path 'extended_annotations.NDR_filter.gtf', emit: ndr_filter_bambu_gtf, optional: true

  script:
  """
  run_bambu.R \
    --tag=. \
    --ncore=${params.maxCpu} \
    --annotation=${ref} \
    --fasta=${fa} \
    --bambu_strand=${params.bambu_strand} \
    --bambu_singleexon=${params.bambu_singleexon} \
    --bambu_rec_ndr=${params.bambu_rec_ndr} \
    *.bam

  sed -i 's/*/./g' extended_annotations.gtf

  if [ "${params.bambu_rec_ndr}" == "false" ]; then
    echo "1" > rec_ndr.txt
  fi
  """
}