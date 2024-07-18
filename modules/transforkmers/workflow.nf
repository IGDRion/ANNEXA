include { EXTRACT_TSS_REGIONS   } from './extract_regions.nf'
include { EXTRACT_TSS_SEQUENCES } from './extract_sequences.nf'
include { PREDICT               } from './predict.nf'
include { FILTER                } from './filter.nf'
include { ADD_CLASS_CODE        } from '../add_class_code.nf'

workflow TFKMERS {
  take:
    novel_gtf
    ref_fa
    bambu_ndr
    tokenizer
    model
    counts_tx
    class_code

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

    ADD_CLASS_CODE(class_code, FILTER.out.gtf)

  emit:
    gtf = FILTER.out.gtf
}
