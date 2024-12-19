include { ASSEMBLE                       } from './stringtie_assemble.nf'
include { MERGE                          } from './stringtie_merge.nf'
include { QUANTIFY                       } from './stringtie_quant.nf'
include { GFFCOMPARE                     } from '../gffcompare/gffcompare.nf'
include { FORMAT_GFFCOMPARE              } from './stringtie_format_gffcompare.nf'
include { READ_LENGTH                    } from './stringtie_read_length.nf'
include { FORMAT_QUANT                   } from './stringtie_format_quant.nf'
include { MERGE_COUNTS                   } from './stringtie_merge_counts.nf'

workflow STRINGTIE {
    take:
    samples
    input_gtf
    ref_fa

    main:
    ASSEMBLE(
        samples, 
        input_gtf)

    MERGE(
        ASSEMBLE.out.stringtie_assemble_gtf.collect(), 
        input_gtf)

    GFFCOMPARE(
        input_gtf, 
        ref_fa, 
        MERGE.out.stringtie_merged_gtf)

    FORMAT_GFFCOMPARE(
        MERGE.out.stringtie_merged_gtf,
        GFFCOMPARE.out.tracking_file)

    READ_LENGTH(samples)

    QUANTIFY(
        READ_LENGTH.out,
        MERGE.out.stringtie_merged_gtf)

    FORMAT_QUANT(
        QUANTIFY.out.counts_gene,
        QUANTIFY.out.counts_transcript,
        MERGE.out.stringtie_merged_gtf
    )

    MERGE_COUNTS(
        FORMAT_QUANT.out.counts_gene.collect(), 
        FORMAT_QUANT.out.counts_transcript.collect())

    emit:
    stringtie_gtf = FORMAT_GFFCOMPARE.out.stringtie_gtf
    class_code_gtf = GFFCOMPARE.out.class_code_gtf
    gene_counts = MERGE_COUNTS.out.gene_counts
    tx_counts = MERGE_COUNTS.out.tx_counts
    ndr = MERGE_COUNTS.out.ndr
}