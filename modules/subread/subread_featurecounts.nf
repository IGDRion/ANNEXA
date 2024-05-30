process SUBREAD_FEATURECOUNTS {
    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'Singularity' ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'quay.io/biocontainers/subread:2.0.1--hed695b0_0' }"
    cpus params.maxCpu

    input:
    path bam
    path ref

    output:
    path "${bam.baseName}_counts_gene.txt"              , emit: gene_counts
    path "${bam.baseName}_counts_transcript.txt"        , emit: tx_counts
    path "${bam.baseName}_o_counts_gene.txt.summary"
    path "${bam.baseName}_o_counts_transcript.txt.summary"
    
    shell:
    '''
    featureCounts \
        -L \
        -O \
        -g gene_id \
        -t exon \
        -T !{params.maxCpu} \
        -a !{ref} \
        -o !{bam.baseName}_o_counts_gene.txt \
        !{bam}

    featureCounts \
        -L \
        -O \
        --primary \
        --fraction \
        -F GTF \
        -g transcript_id \
        -t transcript \
        --extraAttributes gene_id \
        -T !{params.maxCpu} \
        -a !{ref} \
        -o !{bam.baseName}_o_counts_transcript.txt \
        !{bam}
    
    # Remove columns of output to have Bambu-like format
    #cut -f2,3,4,5,6 --complement o_counts_gene.txt > counts_gene.txt
    #cut -f2,3,4,5,6 --complement o_counts_transcript.txt > counts_transcript.txt

    # Remove first line (command info) to have correct header
    #sed -i '1d' counts_gene.txt
    #sed -i '1d' counts_transcript.txt
    
    # Same but awk version
    # Select columns, remove first line, restore tab separators
    awk 'NR>1 {printf("\\n%s ",$1); for (i=7;i<=NF;i++){printf("%s ",$i)}}' !{bam.baseName}_o_counts_gene.txt | grep -v -e '^$' | awk -v OFS="\\t" '$1=$1'> !{bam.baseName}_counts_gene.txt 
    awk 'NR>1 {printf("\\n%s ",$1); for (i=7;i<=NF;i++){printf("%s ",$i)}}' !{bam.baseName}_o_counts_transcript.txt | grep -v -e '^$' | awk -v OFS="\\t" '$1=$1' > !{bam.baseName}_counts_transcript.txt

    # Remove .bam from sample name in column header
    sed -i "s/\\.bam//" !{bam.baseName}_counts_gene.txt
    sed -i "s/\\.bam//" !{bam.baseName}_counts_transcript.txt
    '''
}
