process TRANSDECODER_PREDICT {
    conda (params.enable_conda ? "bioconda::transdecoder=5.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4' :
        'quay.io/comp-bio-aging/transdecoder' }"
    publishDir "$params.outdir/transdecoder", mode: 'copy'
    
    input:
    path novel_full_gtf
    path fa

    output:
    path 'novel.fasta.transdecoder.genome.gff3', emit: gff3

    script:
    """
    # Filter annotation to protein coding only
    grep "protein_coding" ${novel_full_gtf} > pt_coding.gtf

    # Construct transcript fasta file from genome
    gtf_genome_to_cdna_fasta.pl \
        pt_coding.gtf \
        ${fa} \
        > novel.fasta
    
    # Convert transcript structure GTF to alignment-GFF3 formatted file
    gtf_to_alignment_gff3.pl \
        pt_coding.gtf > novel.gff3
    
    # Extract long ORF and predict likely coding regions
    TransDecoder.LongOrfs \
        -S \
        -t novel.fasta
    
    TransDecoder.Predict \
        --single_best_only \
        -t novel.fasta
    
    # Generated genome-based coding region annotation file
    cdna_alignment_orf_to_genome_orf.pl \
        novel.fasta.transdecoder.gff3 \
        novel.gff3 \
        novel.fasta > novel.fasta.transdecoder.genome.gff3
    """
}