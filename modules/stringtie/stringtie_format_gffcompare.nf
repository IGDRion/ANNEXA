process FORMAT_GFFCOMPARE {

  input:
  path extended_annotations
  path merged_gtf
  path tracking_file

  output:
  

  shell:
  '''
  # Add information for stringtie2 process
  # Reformat the output of gffcompare to correctly match novel isoforms to known genes
  # Takes the transcript_id identified by Stringtie and assigns it to reference gene_id

  awk 'BEGIN{
      while(getline<"gffcmp.tracking">0){
          if ($4 !="u" && $4 !="r"){
              split($3,gn,"|");
              split($5,tx,"|"); 
              final["\\""tx[2]"\\";"]="\\""gn[1]"\\";"
              }
          }
  } {
      if ($12 in final){
          $10=final[$12]; print $0} else {print $0}
  }' !{merged_gtf} | gtf2gtf_cleanall.sh  > extended_annotations_preaclean.gtf

  # Match correct ref_gene_id to gene_id to some overlapping genes in the reference annotation

  awk '{if ($3 == "transcript" && $13=="ref_gene_id" && $10!=$14) {
      $10 = $14;
      print $0
  } else if ($3 == "exon" && $15=="ref_gene_id" && $10!=$16) {
      $10 = $16;
      print $0
  } else {print $0}
  }' extended_annotations_preaclean.gtf | gtf2gtf_cleanall.sh > extended_annotations.gtf

  # Remove header lines (command and version)
  sed -i 1,2d extended_annotations.gtf
  '''
}
