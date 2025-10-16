process REFORMAT_MERGED_TOOLS {
  conda (params.enable_conda ? "$baseDir/environment.yml" : null)
  container "ghcr.io/igdrion/annexa:${workflow.revision? workflow.revision: "main"}"

  input:
  file merged_tools
  file tracking

  output:
  path "reformated.gtf", emit: reformated_merged_gtf

  shell:
  '''
  reformat_merge.R
  sed -i '/;\s*$/!s/\s*$/;/' reformated.gtf
  '''
}