include { EXTRACT_TSS_REGIONS   } from './extract_regions.nf'
include { EXTRACT_TSS_SEQUENCES } from './extract_sequences.nf'
include { PREDICT               } from './predict.nf'
include { FILTER                } from './filter.nf'

workflow TFKMERS {
  take:
    novel_gtf
    ref_fa
    bambu_ndr
    tokenizer
    model
    counts_tx

  main:
    EXTRACT_TSS_REGIONS(
      novel_gtf,
    )

    EXTRACT_TSS_SEQUENCES(
      EXTRACT_TSS_REGIONS.out.tss_regions,
      ref_fa
    )

    PREDICT(
      EXTRACT_TSS_SEQUENCES.out.tss_sequences,
      tokenizer,
      model
    )

    FILTER(
      novel_gtf,
      counts_tx,
      PREDICT.out,
      bambu_ndr
    )

  emit:
    gtf = FILTER.out.gtf
}
