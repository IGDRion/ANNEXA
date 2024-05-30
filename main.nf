///////////////////////////////////////////////////////////////////////////
// PARSE ARGS
///////////////////////////////////////////////////////////////////////////
if (params.input) { input = file(params.input, checkIfExists: true) }
else { exit 1, "Input file not specified!" }

if (params.gtf) { ref_gtf = file(params.gtf, checkIfExists: true) }
else { exit 1, "Reference annotation file not specified!" }

if (params.fa) { ref_fa = file(params.fa, checkIfExists: true) }
else { exit 1, "Reference genome file not specified!" }

if (params.tx_discovery != 'bambu' && params.tx_discovery != 'stringtie2') {
  exit 1, "Please specify a valid quantification method ('bambu' (default) or 'stringtie2')."
}

if (params.filter) {
  if (params.tfkmers_model) {
    model = Channel.fromPath(params.tfkmers_model, checkIfExists: true)
  } else { exit 1, "Please specify a valid transforkmers model path."}

  if (params.tfkmers_tokenizer) {
    tokenizer = Channel.fromPath(params.tfkmers_tokenizer, checkIfExists: true)
  } else { exit 1, "Please specify a valid transforkmers tokenizer path."}
}

include { logHeader           } from './modules/header.nf'
log.info logHeader(params)

///////////////////////////////////////////////////////////////////////////
// WORKFLOW
///////////////////////////////////////////////////////////////////////////
nextflow.enable.dsl=2
include { VALIDATE_INPUT_GTF             } from './modules/input/validate.nf'
include { INDEX_BAM                      } from './modules/index_bam.nf'
include { BAMBU                          } from './modules/bambu/bambu.nf'
include { STRINGTIE_ASSEMBLE             } from './modules/stringtie/stringtie_assemble.nf'
include { STRINGTIE_MERGE                } from './modules/stringtie/stringtie_merge.nf'
include { STRINGTIE_QUANTIFY             } from './modules/stringtie/stringtie_quant.nf'
include { GFFCOMPARE                     } from './modules/gffcompare/gffcompare.nf'
include { SPLIT_EXTENDED_ANNOTATION      } from './modules/split_extended_annotation.nf'
include { SUBREAD_FEATURECOUNTS          } from './modules/subread/subread_featurecounts.nf'
include { MERGE_COUNTS                   } from './modules/merge_stringtie_quant.nf'
include { FEELNC_CODPOT                  } from './modules/feelnc/codpot.nf'
include { FEELNC_FORMAT                  } from './modules/feelnc/format.nf'
include { RESTORE_BIOTYPE                } from './modules/restore_biotypes.nf'
include { MERGE_NOVEL                    } from './modules/merge_novel.nf'
include { TFKMERS                        } from './modules/transforkmers/workflow.nf'
include { QC as QC_FULL; QC as QC_FILTER } from './modules/qc/workflow.nf'

workflow {
  ///////////////////////////////////////////////////////////////////////////
  // PROCESS INPUT FILES
  ///////////////////////////////////////////////////////////////////////////
  samples = Channel
    .fromPath(input)
    .splitCsv()
    .map { it ->
            workflow.profile.contains('test') ?
              file("${baseDir}/${it[0]}", checkIfExists: true) :
              file(it[0], checkIfExists: true) }

  VALIDATE_INPUT_GTF(ref_gtf)
  INDEX_BAM(samples)

  ///////////////////////////////////////////////////////////////////////////
  // NEW TRANSCRIPTS DISCOVERY
  ///////////////////////////////////////////////////////////////////////////
  if(params.tx_discovery == "bambu") {
    BAMBU(samples.collect(), VALIDATE_INPUT_GTF.out, ref_fa)
    SPLIT_EXTENDED_ANNOTATION(BAMBU.out.bambu_gtf)
  }
  else if (params.tx_discovery == "stringtie2") {
    STRINGTIE_ASSEMBLE(samples, VALIDATE_INPUT_GTF.out)
    STRINGTIE_MERGE(STRINGTIE_ASSEMBLE.out.stringtie_assemble_gtf.collect(), VALIDATE_INPUT_GTF.out)
    STRINGTIE_QUANTIFY(samples, STRINGTIE_MERGE.out.stringtie_merged_gtf)
    GFFCOMPARE(VALIDATE_INPUT_GTF.out, ref_fa, STRINGTIE_MERGE.out.stringtie_merged_gtf)
    SUBREAD_FEATURECOUNTS(samples, GFFCOMPARE.out.stringtie_gtf)
    MERGE_COUNTS(SUBREAD_FEATURECOUNTS.out.gene_counts.collect(), SUBREAD_FEATURECOUNTS.out.tx_counts.collect())
    SPLIT_EXTENDED_ANNOTATION(GFFCOMPARE.out.stringtie_gtf)
  }

  ///////////////////////////////////////////////////////////////////////////
  // EXTRACT AND CLASSIFY NEW TRANSCRIPTS, AND PERFORM QC
  ///////////////////////////////////////////////////////////////////////////
  FEELNC_CODPOT(VALIDATE_INPUT_GTF.out, ref_fa, SPLIT_EXTENDED_ANNOTATION.out.novel_genes)
  FEELNC_FORMAT(FEELNC_CODPOT.out.mRNA, FEELNC_CODPOT.out.lncRNA)
  RESTORE_BIOTYPE(VALIDATE_INPUT_GTF.out, SPLIT_EXTENDED_ANNOTATION.out.novel_isoforms)
  MERGE_NOVEL(FEELNC_FORMAT.out, RESTORE_BIOTYPE.out)

  if(params.tx_discovery == "bambu") {
    ch_gene_counts = BAMBU.out.gene_counts
    ch_tx_counts = BAMBU.out.tx_counts
    ch_ndr = BAMBU.out.ndr
  }
  else if (params.tx_discovery == "stringtie2") {
    ch_gene_counts = MERGE_COUNTS.out.gene_counts
    ch_tx_counts = MERGE_COUNTS.out.tx_counts
    ch_ndr = MERGE_COUNTS.out.ndr
  }

  QC_FULL(samples, 
          INDEX_BAM.out, 
          MERGE_NOVEL.out, 
          VALIDATE_INPUT_GTF.out, 
          ch_gene_counts,
          "full")

  ///////////////////////////////////////////////////////////////////////////
  // FILTER NEW TRANSCRIPTS, AND QC ON FILTERED ANNOTATION
  ///////////////////////////////////////////////////////////////////////////
  if(params.filter) {
    TFKMERS(MERGE_NOVEL.out, 
            ref_fa, 
            ch_ndr, 
            tokenizer, 
            model, 
            ch_tx_counts)
    QC_FILTER(samples,
              INDEX_BAM.out, 
              TFKMERS.out.gtf, 
              VALIDATE_INPUT_GTF.out, 
              ch_gene_counts, 
              "filter")
  }
}
