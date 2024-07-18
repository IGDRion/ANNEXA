include { TRANSDECODER_PREDICT  } from './transdecoder_predict.nf'
include { FORMAT_TRANSDECODER   } from './format_transdecoder.nf'
include { AGAT_CONVERTSPGFF2GTF } from './convert_to_gtf.nf'
include { GTFSORT               } from './gtfsort.nf'
include { ADD_CLASS_CODE        } from '../add_class_code.nf'

workflow TRANSDECODER {
    take:
    novel_full_gtf
    ref_fa
    class_code

    main:
    TRANSDECODER_PREDICT(
        novel_full_gtf, 
        ref_fa)

    FORMAT_TRANSDECODER(
        TRANSDECODER_PREDICT.out.gff3, 
        novel_full_gtf)
        
    AGAT_CONVERTSPGFF2GTF(
        FORMAT_TRANSDECODER.out.exon_gff3)
        
    GTFSORT(
        AGAT_CONVERTSPGFF2GTF.out.td_gtf, 
        FORMAT_TRANSDECODER.out.fixed_novel)

    ADD_CLASS_CODE(
        class_code,
        GTFSORT.out,
    )

    emit:
    cds_gtf = GTFSORT.out.gtf
}