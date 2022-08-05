///////////////////////////////////////////////////////////////////////////
// PARSE ARGS
///////////////////////////////////////////////////////////////////////////
if (params.input) { input = file(params.input, checkIfExists: true) }
else { exit 1, "Input file not specified!" }

if (params.gtf) { ref_gtf = file(params.gtf, checkIfExists: true) }
else { exit 1, "Reference annotation file not specified!" }

if (params.fa) { ref_fa = file(params.fa, checkIfExists: true) }
else { exit 1, "Reference genome file not specified!" }


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
include { BAMBU_SPLIT_RESULTS            } from './modules/bambu/split.nf'
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
  BAMBU(samples.collect(), VALIDATE_INPUT_GTF.out, ref_fa)
  BAMBU_SPLIT_RESULTS(BAMBU.out.bambu_gtf)

  ///////////////////////////////////////////////////////////////////////////
  // EXTRACT AND CLASSIFY NEW TRANSCRIPTS, AND PERFORM QC
  ///////////////////////////////////////////////////////////////////////////
  FEELNC_CODPOT(VALIDATE_INPUT_GTF.out, ref_fa, BAMBU_SPLIT_RESULTS.out.novel_genes)
  FEELNC_FORMAT(FEELNC_CODPOT.out.mRNA, FEELNC_CODPOT.out.lncRNA)
  RESTORE_BIOTYPE(VALIDATE_INPUT_GTF.out, BAMBU_SPLIT_RESULTS.out.novel_isoforms)
  MERGE_NOVEL(FEELNC_FORMAT.out, RESTORE_BIOTYPE.out)

  QC_FULL(samples, 
          INDEX_BAM.out, 
          MERGE_NOVEL.out, 
          VALIDATE_INPUT_GTF.out, 
          BAMBU.out.gene_counts, 
          "full")

  ///////////////////////////////////////////////////////////////////////////
  // FILTER NEW TRANSCRIPTS, AND QC ON FILTERED ANNOTATION
  ///////////////////////////////////////////////////////////////////////////
  if(params.filter) {
    TFKMERS(MERGE_NOVEL.out, ref_fa, BAMBU.out.ndr, 
            tokenizer, model, BAMBU.out.tx_counts)
    QC_FILTER(samples,
              INDEX_BAM.out, 
              TFKMERS.out.gtf, 
              VALIDATE_INPUT_GTF.out, 
              BAMBU.out.gene_counts, 
              "filter")
  }
}
