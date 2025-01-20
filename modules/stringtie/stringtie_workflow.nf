include { ASSEMBLE                       } from './stringtie_assemble.nf'
include { MERGE                          } from './stringtie_merge.nf'
include { QUANTIFY                       } from './stringtie_quant.nf'
include { GFFCOMPARE                     } from '../gffcompare/gffcompare.nf'
include { FORMAT_GFFCOMPARE              } from './stringtie_format_gffcompare.nf'
include { READ_LENGTH                    } from './stringtie_read_length.nf'
include { EXTRACT_QUANTS                 } from './stringtie_extract_counts.nf'

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
        READ_LENGTH.out.bam_length,
        MERGE.out.stringtie_merged_gtf)
    
    EXTRACT_QUANTS(
        READ_LENGTH.out.av_length.collect(),
        QUANTIFY.out.collect(),
        MERGE.out.stringtie_merged_gtf
    )

    emit:
    stringtie_gtf = FORMAT_GFFCOMPARE.out.stringtie_gtf
    class_code_gtf = GFFCOMPARE.out.class_code_gtf
    gene_counts = EXTRACT_QUANTS.out.gene_counts
    tx_counts = EXTRACT_QUANTS.out.tx_counts
    ndr = EXTRACT_QUANTS.out.ndr
}