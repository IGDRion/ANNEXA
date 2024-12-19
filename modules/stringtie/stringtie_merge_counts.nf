process MERGE_COUNTS {
  publishDir "$params.outdir/stringtie2", mode: 'copy', pattern: '*.txt'
  if (params.filter == false){
    publishDir "$params.outdir/final", mode: 'copy', pattern: 'counts_transcript.txt', saveAs: {filename -> 'counts_transcript.full.txt'}
  }

  input:
  path gene_counts
  path tx_counts

  output:
  path "counts_gene.txt", emit: gene_counts
  path "counts_transcript.txt", emit: tx_counts
  path "empty.ndr", emit: ndr

  shell:
  '''
  # Sort list of input files alphanumerically
  gene_counts=$(echo !{gene_counts} | tr ' ' '\\n' | sort | xargs)
  tx_counts=$(echo !{tx_counts} | tr ' ' '\\n' | sort | xargs)
  
  # Merge the individual outputs of prepDE of each .bam into a single file

  paste \${gene_counts} \
  | awk '{printf("%s ",$1); for (i=2;i<=NF;i+=2){printf("%s ",$i)}print "\\n"}' \
  | grep -v -e "^$"  \
  | awk -v OFS='\\t' '{$1=$1}1' > counts_gene.txt

  paste \${tx_counts} \
  | awk '{printf("\\n%s %s ",$1,$2); for (i=3;i<=NF;i+=3){printf("%s ",$i)}}' \
  | grep -v -e "^$"  \
  | awk -v OFS='\\t' '{$1=$1}1' > counts_transcript.txt

  # Rename first line
  sed -i '1s/^Geneid/transcript_id/' counts_transcript.txt

  # Sort genes and tx rows alphanumerically
  (head -n 1 counts_transcript.txt && tail -n +2 counts_transcript.txt | sort -k1,1) > counts_transcript.txt.new \
  && mv counts_transcript.txt.new counts_transcript.txt

  (head -n 1 counts_gene.txt && tail -n +2 counts_gene.txt | sort -k1,1) > counts_gene.txt.new \
  && mv counts_gene.txt.new counts_gene.txt

  # Create empty NDR file for TFKMERS to have same workflow as bambu in filtering step 
  # (placeholder, not used)
  touch empty.ndr
  '''
}