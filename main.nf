///////////////////////////////////////////////////////////////////////////
// PARSE ARGS
///////////////////////////////////////////////////////////////////////////

// Help Message
include { paramsHelp } from 'plugin/nf-schema'
include { helpHeader           } from './modules/header.nf'
if (params.help) {
  log.info helpHeader()
  log.info paramsHelp("nextflow run IGDRion/ANNEXA -profile {test,docker,singularity,conda,slurm} --input samples.txt --gtf /path/to/ref.gtf --fa /path/to/ref.fa")
  exit 0
}

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

if (params.tx_discovery == "bambu"){
  log.warn """You are using Bambu as the transcript discovery tool. Please note that Bambu takes by default strandedness into account.
  Bambu using strandedness: ${params.bambu_strand}. Please verifiy that it corresponds to your sequencing protocol.
  You can change strandedness mode using the --bambu_strand parameter (options: true or false).
  """
}

///////////////////////////////////////////////////////////////////////////
// WORKFLOW
///////////////////////////////////////////////////////////////////////////
nextflow.enable.dsl=2
include { VALIDATE_INPUT_GTF             } from './modules/input/validate.nf'
include { INDEX_BAM                      } from './modules/index_bam.nf'
include { BAMBU                          } from './modules/bambu/bambu.nf'
include { STRINGTIE                      } from './modules/stringtie/stringtie_workflow.nf'
include { GFFCOMPARE                     } from './modules/gffcompare/gffcompare.nf'
include { RESTRAND_ISOFORMS              } from './modules/restrand_isoforms.nf'
include { SPLIT_EXTENDED_ANNOTATION      } from './modules/split_extended_annotation.nf'
include { FEELNC_CODPOT                  } from './modules/feelnc/codpot.nf'
include { FEELNC_FORMAT                  } from './modules/feelnc/format.nf'
include { RESTORE_BIOTYPE                } from './modules/restore_biotypes.nf'
include { MERGE_NOVEL                    } from './modules/merge_novel.nf'
include { TRANSDECODER                   } from './modules/transdecoder/transdecoder_workflow.nf'
include { RESTRAND_NOVEL                 } from './modules/restrand_novel.nf'
include { TFKMERS                        } from './modules/transforkmers/workflow.nf'
include { QC as QC_FULL; QC as QC_FILTER } from './modules/qc/workflow.nf'
include { ADD_CLASS_CODE                 } from './modules/add_class_code.nf'

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
    GFFCOMPARE(VALIDATE_INPUT_GTF.out, ref_fa, BAMBU.out.bambu_gtf)
    RESTRAND_ISOFORMS(BAMBU.out.bambu_gtf)
    SPLIT_EXTENDED_ANNOTATION(RESTRAND_ISOFORMS.out)
  }
  else if (params.tx_discovery == "stringtie2") {
    STRINGTIE(samples, VALIDATE_INPUT_GTF.out, ref_fa)
    RESTRAND_ISOFORMS(STRINGTIE.out.stringtie_gtf)
    SPLIT_EXTENDED_ANNOTATION(RESTRAND_ISOFORMS.out)
  }

  ///////////////////////////////////////////////////////////////////////////
  // EXTRACT AND CLASSIFY NEW TRANSCRIPTS
  ///////////////////////////////////////////////////////////////////////////
  FEELNC_CODPOT(VALIDATE_INPUT_GTF.out, ref_fa, SPLIT_EXTENDED_ANNOTATION.out.novel_genes)
  FEELNC_FORMAT(FEELNC_CODPOT.out.mRNA, FEELNC_CODPOT.out.lncRNA)
  RESTORE_BIOTYPE(VALIDATE_INPUT_GTF.out, SPLIT_EXTENDED_ANNOTATION.out.novel_isoforms)
  MERGE_NOVEL(FEELNC_FORMAT.out, RESTORE_BIOTYPE.out)

  if(params.tx_discovery == "bambu") {
    ch_gene_counts = BAMBU.out.gene_counts
    ch_tx_counts = BAMBU.out.tx_counts
    ch_ndr = BAMBU.out.ndr
    ch_rec_ndr = BAMBU.out.rec_ndr
    class_code = GFFCOMPARE.out.class_code_gtf
  }
  else if (params.tx_discovery == "stringtie2") {
    ch_gene_counts = STRINGTIE.out.gene_counts
    ch_tx_counts = STRINGTIE.out.tx_counts
    ch_ndr = STRINGTIE.out.ndr
    ch_rec_ndr = STRINGTIE.out.rec_ndr
    class_code = STRINGTIE.out.class_code_gtf
  }

  ///////////////////////////////////////////////////////////////////////////
  // PREDICT CDS ON NOVEL TRANSCRIPTS
  ///////////////////////////////////////////////////////////////////////////
  TRANSDECODER(MERGE_NOVEL.out.novel_full_gtf, ref_fa)
  RESTRAND_NOVEL(TRANSDECODER.out)
  ///////////////////////////////////////////////////////////////////////////
  // PERFORM QC ON FULL ANNOTATION
  ///////////////////////////////////////////////////////////////////////////
  QC_FULL(samples, 
          INDEX_BAM.out, 
          RESTRAND_NOVEL.out, 
          VALIDATE_INPUT_GTF.out, 
          ch_gene_counts,
          "full")

  ///////////////////////////////////////////////////////////////////////////
  // FILTER NEW TRANSCRIPTS, AND QC ON FILTERED ANNOTATION
  ///////////////////////////////////////////////////////////////////////////
  if(params.filter) {
    TFKMERS(RESTRAND_NOVEL.out, 
            ref_fa, 
            ch_ndr, 
            ch_rec_ndr, 
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

  ///////////////////////////////////////////////////////////////////////////
  // ADD GFFCOMPARE CLASS CODES TO FINAL GTFS
  ///////////////////////////////////////////////////////////////////////////
  final_gtf = RESTRAND_NOVEL.out.mix(QC_FULL.out.gtf)
  if (params.filter){
  final_gtf = RESTRAND_NOVEL.out.mix(QC_FULL.out.gtf,TFKMERS.out.gtf,QC_FILTER.out.gtf)
  }
  ADD_CLASS_CODE(class_code, final_gtf)
}