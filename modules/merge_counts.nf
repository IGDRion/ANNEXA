process MERGE_COUNTS {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  if (params.filter == false){
    publishDir "$params.outdir/final", mode: 'copy', pattern: 'combined_tx.txt', saveAs: {filename -> 'counts_transcript.full.txt'}
  }
  
  input:
  path bambu_gene,      stageAs: 'bambu_gene.txt'
  path stringtie_gene,  stageAs: 'stringtie_gene.txt'
  path bambu_tx,        stageAs: 'bambu_tx.txt'
  path stringtie_tx,    stageAs: 'stringtie_tx.txt'
  path bambu_cc,        stageAs: 'bambu_cc.gtf'
  path stringtie_cc,    stageAs: 'stringtie_cc.gtf'

  output:
  path "combined_genes.txt", emit: gene_counts
  path "combined_tx.txt", emit: tx_counts
  path "combined_cc.gtf", emit: class_code_gtf

  script:
  """
  grep "MSTRG" stringtie_gene.txt > stringtie_genes.txt
  cat bambu_gene.txt stringtie_genes.txt > combined_genes.txt

  grep "MSTRG" stringtie_tx.txt > stringtie_txs.txt
  cat bambu_tx.txt stringtie_txs.txt > combined_tx.txt

  grep "MSTRG" stringtie_cc.gtf > stringtie_ccs.gtf
  cat bambu_cc.gtf stringtie_ccs.gtf > combined_cc.gtf
  """
}