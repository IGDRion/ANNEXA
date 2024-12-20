process AGAT_CONVERTSPGFF2GTF {
    conda (params.enable_conda ? "bioconda::agat=1.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.0--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.0--pl5321hdfd78af_0' }"

    input:
    path gff3

    output:
    path 'converted_transdecoder.gtf', emit: td_gtf

    script:
    """
    agat_convert_sp_gff2gtf.pl \
        --gff ${gff3} \
        --gtf_version 3 \
        --output exon_cds.gtf

    # Remove redundant information from gtf
    sed -i 's/ID "\\(.*\\)\\.exon\\([0-9]*\\)"/exon_number "\\2"/g' exon_cds.gtf

    # Remove exons (already present in Stringtie output)
    awk '\$3=="CDS" || \$3=="five_prime_utr" || \$3=="three_prime_utr"' exon_cds.gtf > converted_transdecoder.gtf
    """
}