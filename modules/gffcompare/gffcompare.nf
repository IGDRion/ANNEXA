process GFFCOMPARE {
  conda (params.enable_conda ? "bioconda::gffcompare" : null)
  container "${ workflow.containerEngine == 'singularity' ?
      'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0' :
      'biocontainers/gffcompare:0.12.6--h9f5acd7_0' }"
  publishDir "$params.outdir/stringtie2", mode: 'copy', pattern: 'extended_annotations.gtf'
  cpus params.maxCpu

  input:
    input:
    path reference_gtf
    path fasta
    path merged_gtf

    output:
    path("*.annotated.gtf"), emit: annotated_gtf
    path("*.loci")         , emit: loci
    path("*.stats")        , emit: stats
    path("*.tracking")     , emit: tracking
    path("extended_annotations.gtf"), emit: stringtie_gtf

    shell:
    '''
    gffcompare \
    -r !{reference_gtf} \
    -s !{fasta} \
    !{merged_gtf}

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
    '''
}