include { MERGE_ANNOTATIONS } from './merge_known_novel.nf'
include { REPORT            } from './report.nf'
include { RSEQC             } from '../rseqc/workflow.nf'

workflow QC {
  take:
    bam
    bai
    novel_gtf
    ref_gtf
    counts_gene
    origin

  main:
    MERGE_ANNOTATIONS(novel_gtf, ref_gtf, origin)

    REPORT(
      MERGE_ANNOTATIONS.out,
      ref_gtf,
      counts_gene,
      origin
    )

    if (params.withGeneCoverage) {
      RSEQC(
        bam,
        bai,
        novel_gtf,
        ref_gtf,
        origin
      )
    }
}
