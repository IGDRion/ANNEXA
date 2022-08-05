process PREDICT {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "mlorthiois/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/transforkmers", mode: 'copy'

  cpus params.maxCpu
  memory '16GB'

  input:
  file tss_sequences
  path tokenizer
  path model

  output:
  path "output.csv", emit: tss_prob

  """
  transforkmers predict \
    --model_path_or_name ${model} \
    --tokenizer ${tokenizer} \
    --input ${tss_sequences} \
    --quantize-model \
    --output . \
    --per_device_eval_batch_size 64
  """
}
