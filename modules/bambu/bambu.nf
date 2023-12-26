process BAMBU {
  conda (params.enable_conda ? "bioconda::bioconductor-bambu=3.0.8" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.0.8--r42hc247a5b_0' :
                'quay.io/biocontainers/bioconductor-bambu:3.0.8--r42hc247a5b_0' }"
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

  """
  run_bambu.R \
    --tag=. \
    --ncore=${params.maxCpu} \
    --annotation=${ref} \
    --fasta=${fa} \
    *.bam

  sed -i 's/*/./g' extended_annotations.gtf
  """
}
