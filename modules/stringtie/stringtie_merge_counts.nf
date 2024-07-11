process MERGE_COUNTS {
  publishDir "$params.outdir/stringtie2", mode: 'copy', pattern: '*.txt'

  input:
  path gene_counts
  path tx_counts

  output:
  path "counts_gene.txt", emit: gene_counts
  path "counts_transcript.txt", emit: tx_counts
  path "empty.ndr", emit: ndr

  shell:
  '''
  # Merge the individual outputs of featurecount of each .bam into a single file

  paste !{gene_counts} \
  | awk '{printf("%s ",$1); for (i=2;i<=NF;i+=2){printf("%s ",$i)}print "\\n"}' \
  | grep -v -e "^$"  \
  | awk -v OFS='\\t' '{$1=$1}1' > counts_gene.txt

  paste !{tx_counts} \
  | awk '{printf("\\n%s %s ",$1,$2); for (i=3;i<=NF;i+=3){printf("%s ",$i)}}' \
  | grep -v -e "^$"  \
  | awk -v OFS='\\t' '{$1=$1}1' > counts_transcript.txt

  # Create empty NDR file for TFKMERS to have same workflow as bambu in filtering step 
  # (placeholder, not used)
  touch empty.ndr
  '''
}
