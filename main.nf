if (params.input) { ch_input = file(params.input, checkIfExists: true) }
else { exit 1, "Samplesheet file not specified!" }

if (params.gtf) { ch_ref = file(params.gtf, checkIfExists: true) }
else { exit 1, "Reference Annotation file not specified!" }

if (params.fa) { ch_fa = file(params.fa, checkIfExists: true) }
else { exit 1, "Reference Genome file not specified!" }

params.maxCpu = 8

// Prevent typo
if (params.profile) exit 1, "Maybe you want to use -profile instead of --profile?"

def logHeader() {
    // Log colors ANSI codes
    c_dim = "\033[2m";
    c_green = "\033[0;32m";
    c_purple = "\033[0;35m";
    c_reset = "\033[0m";

    return """-${c_dim}----------------------------------------${c_reset}-
${c_green}    ___    _   ___   _________  __ ___ 
   /   |  / | / / | / / ____/ |/ //   |
  / /| | /  |/ /  |/ / __/  |   // /| |
 / ___ |/ /|  / /|  / /___ /   |/ ___ |
/_/  |_/_/ |_/_/ |_/_____//_/|_/_/  |_|
                                       ${c_reset}
-${c_dim}----------------------------------------${c_reset}-
${c_purple}github.com/mlorthiois/ANNEXA${c_reset}
Reference Annotation: ${params.gtf}
Reference Genome    : ${params.fa}
Input Samplesheet   : ${params.input}
Gene coverage       : ${params.withGeneCoverage}
-${c_dim}----------------------------------------${c_reset}-
""".stripIndent()
}

log.info logHeader()

///////////////////////////////////////////////////////////////////////////
// Create channels from samples.txt, each bam independently
Channel
  .fromPath( ch_input )
  .splitCsv()
  .map { it -> workflow.profile.contains('test') ? file("${baseDir}/${it[0]}", checkIfExists: true) : file(it[0], checkIfExists: true) }
  .into { ch_bam; ch_bam_rseq; ch_bam_bai }

process FORMAT_INPUT_GTF {
  input:
  file gtf from ch_ref

  output:
  file 'input.formatted.gtf' into ch_format_ref

  """
  format_gtf.py ${gtf} > input.formatted.gtf
  """
}

process BAMBU {
  publishDir "$params.outdir/bambu", mode: 'copy'
  cpus params.maxCpu
  memory '40GB'

  input:
  // Collect parallel bam channel into 1
  file '*' from ch_bam.collect()
  file ref from ch_format_ref
  file fa from ch_fa

  output:
  file 'extended_annotations.gtf' into ch_bambu_gtf
  file 'counts_transcript.txt' into ch_bambu_tx
  file 'counts_gene.txt' into ch_bambu_gene

  """
  bambu_counts.R \
    --tag=. \
    --ncore=${params.maxCpu} \
    --annotation=${ref} \
    --fasta=${fa} \
    *.bam
  """
}


process SPLIT_BAMBU {
  input:
  file input from ch_bambu_gtf

  output:
  file("known.gtf") into ch_known
  file("novel.gtf") into ch_novel

  """
  sed 's/*/./g' ${input} > tmp.gtf
  split_merged_gtf.sh tmp.gtf
  """
}

process FEELNC {
    publishDir "$params.outdir", mode: 'copy', saveAs: { filename -> "feelnc/${filename.split('/')[-1]}" }
    memory '16 GB'

    input:
    file ref from ch_format_ref
    file fa from ch_fa
    file gtf from ch_novel

    output:
    file("feelnc_codpot_out/new.lncRNA.gtf") into ch_feelnc_codpot_lncRNA
    file("feelnc_codpot_out/new.mRNA.gtf") into ch_feelnc_codpot_mRNA

    shell:
    """
    FEELnc_codpot.pl \
        -i $gtf \
        -a $ref \
        -g $fa \
        -b transcript_biotype=protein_coding \
        -b transcript_status=KNOWN \
        --numtx=2000,2000 \
        --mode=shuffle \
        -o new
    """
}

process RESTORE_ATTRIBUTES_FROM_REF {
  input:
  file ref from ch_format_ref
  file known from ch_known

  output:
  file("known.restored.gtf") into ch_known_restored

  """
  restore_from_ref.py -gtf $known -ref $ref > known.restored.gtf
  """
}

process JOIN_FEELNC_WITH_BIOTYPE {
  publishDir "$params.outdir/feelnc", mode: 'copy'

  input:
  file mRNA from ch_feelnc_codpot_mRNA
  file lncRNA from ch_feelnc_codpot_lncRNA

  output:
  file("feelnc.combined.gtf") into ch_feelnc_combined

  """
  merge_feelnc.py -lncRNA ${lncRNA} -mRNA ${mRNA} > feelnc.combined.gtf
  """
}

process FORMAT_NOVEL_KNOWN {
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file known from ch_known_restored
  file novel from ch_feelnc_combined

  output:
  file("extended_annotations.gtf") into (ch_final, ch_final_rseq)

  """
  cat $known $novel | GTF.py format > extended_annotations.gtf
  """
}

process QC_EXTENDED_GTF {
  publishDir "$params.outdir/qc", mode: 'copy'

  input:
  file gtf from ch_final
  file ref from ch_format_ref
  file counts_gene from ch_bambu_gene

  output:
  file "gene.stats" into ch_gene_stats
  file "transcript.stats" into ch_tx_stats
  file "exon.stats" into ch_exon_stats
  file "qc_gtf.pdf" into ch_qc_pdf

  """
  qc_gtf.py -gtf ${gtf} \
    -c_gene ${counts_gene} \
    -ref ${ref}
  qc.R
  """
}

if (params.withGeneCoverage) {
  process PREPARE_RSEQC {
    input:
    file gtf from ch_final_rseq

    output:
    file "*.filter.gtf" into ch_gtf_filter_rseq

    """
    grep "FEELnc" $gtf | grep "protein_coding" > novel_mRNA.filter.gtf
    grep "FEELnc" $gtf | grep "lncRNA" > novel_lncRNA.filter.gtf
    grep -v "FEELnc" $gtf | grep "protein_coding" > known_mRNA.filter.gtf
    grep -v "FEELnc" $gtf | grep "lncRNA" > known_lncRNA.filter.gtf
    """
  }

  process CONVERT_TO_BED12 {
    input:
    // Split 4 files in 1 channel into 4 parallel channels
    file gtf from ch_gtf_filter_rseq.flatten()

    output:
    file "${gtf.simpleName}.bed12" into ch_rseq_bed12

    """
    gtfToGenePred -genePredExt ${gtf} convert.genePred
    genePredToBed convert.genePred ${gtf.simpleName}.bed12
    """
  }

  process CREATE_BAI {
    input:
    file bam from ch_bam_bai

    output:
    file("*.bai") into ch_bam_bai_rseq 

    """
    samtools index $bam
    """
  }

  process GENE_BODY_COVERAGE {
    publishDir "$params.outdir/RSeQC", mode: 'copy'

    input:
    file bed from ch_rseq_bed12
    // Collect parallel bam and bai channels inside 1
    file "*" from ch_bam_bai_rseq.collect()
    file "*" from ch_bam_rseq.collect()

    output:
    file "${bed.simpleName}.geneBodyCoverage.curves.pdf"
    file "${bed.simpleName}.geneBodyCoverage.txt"

    """
    geneBody_coverage.py \
      -i ./ \
      -r ${bed} \
      -o ${bed.simpleName}
    """
  }
}
