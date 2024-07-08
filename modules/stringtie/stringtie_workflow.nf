include { STRINGTIE_ASSEMBLE             } from './stringtie_assemble.nf'
include { STRINGTIE_MERGE                } from './stringtie_merge.nf'
include { STRINGTIE_QUANTIFY             } from './stringtie_quant.nf'
include { GFFCOMPARE                     } from '../gffcompare/gffcompare.nf'
include { SUBREAD_FEATURECOUNTS          } from '../subread/subread_featurecounts.nf'
include { MERGE_COUNTS                   } from '../merge_stringtie_quant.nf'

workflow STRINGTIE {
    take:
    samples
    input_gtf
    ref_fa

    main:
    STRINGTIE_ASSEMBLE(
        samples, 
        input_gtf)

    STRINGTIE_MERGE(
        STRINGTIE_ASSEMBLE.out.stringtie_assemble_gtf.collect(), 
        input_gtf)

    STRINGTIE_QUANTIFY(
        samples, 
        STRINGTIE_MERGE.out.stringtie_merged_gtf)

    GFFCOMPARE(
        input_gtf, 
        ref_fa, 
        STRINGTIE_MERGE.out.stringtie_merged_gtf)

    SUBREAD_FEATURECOUNTS(
        samples, 
        GFFCOMPARE.out.stringtie_gtf)

    MERGE_COUNTS(
        SUBREAD_FEATURECOUNTS.out.gene_counts.collect(), 
        SUBREAD_FEATURECOUNTS.out.tx_counts.collect())

    emit:
    extended_annotation = GFFCOMPARE.out.stringtie_gtf
    gene_counts = MERGE_COUNTS.out.gene_counts
    tx_counts = MERGE_COUNTS.out.tx_counts
    ndr = MERGE_COUNTS.out.ndr
}