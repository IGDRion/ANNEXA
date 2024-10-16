process ADD_CLASS_CODE {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"
  publishDir "$params.outdir/final", mode: 'copy', saveAs: {filename -> "${gtf}"}, overwrite: true
  tag "$gtf"

  input:
  file class_code_gtf
  file gtf

  output:
  path "class_code.${gtf}"

  script:
  """
  class_code.R ${class_code_gtf} ${gtf} "class_code.${gtf}"

  # Remove header created by gtfsort
  sed -i 1,3d "class_code.${gtf}"

  # Add semicolon at end of tx lines
  sed -i '/\\ttranscript\\t/s/\$/;/' "class_code.${gtf}"
  """
}