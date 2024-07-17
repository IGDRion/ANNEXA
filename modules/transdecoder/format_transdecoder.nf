process FORMAT_TRANSDECODER {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"

  input:
  path td_gff
  path novel_full_gtf

  output:
  path 'exon_cds.gff3', emit: exon_gff3
  path 'fixed_novel.full.gtf', emit: fixed_novel

  script:
  """
    ### Transdecoder output
    # Restore exon_number in CDS and create exon_cds.gff3
    transdecoder_add_exon_to_CDS.py --gff ${td_gff}

    # Remove redundant information from final gff3 in gene_name and transcript_id
    sed -i 's/\\^[^ ]*/;/' exon_cds.gff3
    sed -i 's/\\.p[0-9]//'g exon_cds.gff3

    ### Stringtie output
    # Correct exon number in Stringtie output and merge new CDS
    stringtie_fix_exon_number.py ${novel_full_gtf} fixed_novel.full.gtf
  """
}