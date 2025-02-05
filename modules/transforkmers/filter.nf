process FILTER {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/python:3.10.4' : 
                'quay.io/biocontainers/python:3.10.4' }"
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file novel
  file counts_tx
  file tfkmers
  file bambu_ndr
  file rec_ndr

  output:
  path "novel.filter.gtf", emit: gtf
  path "counts_transcript.filter.txt"
  path "counts_transcript.full.txt"

  script:
  """
  cp ${counts_tx} counts_transcript.full.txt

  if [ "${params.bambu_rec_ndr}" == "true" ]; then
    bambu_threshold=\$(cat ${rec_ndr})
  else
    bambu_threshold=${params.bambu_threshold}
  fi

  filter_gtf_ndr.py \
    --gtf ${novel} \
    --counts_tx ${counts_tx} \
    --tfkmers ${tfkmers} \
    --bambu ${bambu_ndr} \
    --tfkmers-threshold ${params.tfkmers_threshold} \
    --bambu-threshold \$bambu_threshold \
    --operation ${params.operation} \
    --tx_discovery ${params.tx_discovery}

  GTF.py format -i unformat.novel.filter.gtf > novel.filter.gtf
  """
}
