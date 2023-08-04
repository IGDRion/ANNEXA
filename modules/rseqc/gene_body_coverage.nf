process GENEBODY_COVERAGE {
  conda (params.enable_conda ? "bioconda::rseqc=4.0" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/rseqc:4.0.0--py310h1425a21_2' : 
                'quay.io/biocontainers/rseqc:4.0.0--py310h1425a21_2' }"
  publishDir "${params.outdir}/qc/${prefix}/rseqc", mode: 'copy'

  input:
  file bed 
  file "*"
  file "*"
  val prefix

  output:
  file "${prefix}.${bed.simpleName}.geneBodyCoverage.curves.pdf"
  file "${prefix}.${bed.simpleName}.geneBodyCoverage.txt"

  """
  geneBody_coverage.py \
    -i ./ \
    -r ${bed} \
    -o ${prefix}.${bed.simpleName}
  """
}
