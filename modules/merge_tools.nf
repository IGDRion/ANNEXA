process MERGE_TOOLS {
    conda (params.enable_conda ? "bioconda::agat=1.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.0--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.0--pl5321hdfd78af_0' }"

    input:
    path bambu_gtf, stageAs: 'bambu.extended_annotations.gtf'
    path stringtie_gtf, stageAs: 'stringtie.extended_annotations.gtf'

    output:
    path("cleaned_merged.gtf"), emit: merged_tools
    
    shell:
    '''
    agat_sp_merge_annotations.pl \\
    -f bambu.extended_annotations.gtf \\
    -f stringtie.extended_annotations.gtf \\
    -o merged.gff

    agat_convert_sp_gff2gtf.pl \\
    -gff merged.gff \\
    -o merged.gtf

    awk 'BEGIN{FS=OFS="\\t"}
    {
    if(\$0 ~ /^#/){ print; next }
    if(\$3 == "gene") next

    split(\$9, a, ";");
    new_attrs = "";
    for(i in a) {
        gsub(/^[ \t]+|[ \t]+$/, "", a[i]);
        if(a[i] ~ /^gene_id/ || a[i] ~ /^transcript_id/ || a[i] ~ /^exon_number/ || a[i] ~ /^ref_gene_id/) {
            new_attrs = new_attrs a[i] "; "
        }
    }
    \$9 = new_attrs;
    print;
    }' merged.gtf > cleaned_merged.gtf
    '''
}
