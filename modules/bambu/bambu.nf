process BAMBU {
  conda (params.enable_conda ? "bioconda::bioconductor-bambu=3.0.8" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.4.0--r43hf17093f_1' :
                'quay.io/biocontainers/bioconductor-bambu:3.4.0--r43hf17093f_1 ' }"
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

  shell:
  '''
  run_bambu.R \
    --tag=. \
    --ncore=!{params.maxCpu} \
    --annotation=!{ref} \
    --fasta=!{fa} \
    --bambu_strand=!{params.bambu_strand} \
    *.bam

  sed -i 's/*/./g' extended_annotations.gtf

  # Optional: Prefilter bambu output to keep tx with NDR < threshold
  if [ !{params.prefilter_ndr} == 'true' ]; then
    awk -F',' '$3 < !{params.bambu_threshold} {print $1}' bambu_ndr.csv > valid_transcripts.txt

    awk '
    BEGIN {
        while (getline < "valid_transcripts.txt") {
            valid[$1] = 1
        }
    }
    {
        if ($0 ~ /transcript_id/) {
            split($0, a, "transcript_id \\"")
            split(a[2], b, "\\"")
            id = b[1]
            if (id in valid || id !~ /BambuTx/) {
                print
            }
        } else {
            print
        }
    }' extended_annotations.gtf > final_output.gtf

    mv final_output.gtf extended_annotations.gtf
    rm valid_transcripts.txt
  fi
  '''
}