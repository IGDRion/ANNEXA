![ANNEXA logo](./img/ANNEXA_logo.png)


## Introduction

**ANNEXA** is an all-in-one reproductible pipeline, written using the [Nextflow](https://nextflow.io) workflow system, which allows users to analyze LR-RNAseq data (Long-Read RNASeq), and to reconstruct and quantify known and novel genes and transcript isoforms.

## Pipeline summary

![Metro map](./assets/metro_map.png)

ANNEXA works by using only three parameter files (a reference genome, a reference annotation and mapping files) and provides users with an extended annotation distinguishing between novel protein-coding (mRNA) versus long non-coding RNAs (lncRNA) genes. All known and novel gene/transcript models are further characterized through multiple features (length, number of spliced transcripts, normalized expression levels,...) available as graphical outputs.

1. Check if the input annotation contains all the information needed.
2. Transcriptome reconstruction and quantification with [bambu](https://github.com/GoekeLab/bambu) or [StringTie](https://github.com/gpertea/stringtie).
3. Novel classification with [FEELnc](https://github.com/tderrien/FEELnc).
4. Retrieve information from input annotation and format final gtf with 3 level structure: gene -> transcript -> exon.
5. Predict the CDS of novel protein-coding transcripts with [TransDecoder](https://github.com/TransDecoder/TransDecoder).
6. Classify novel transcripts with class codes from [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#transfrag-class-codes).
7. Filter novel transcripts based on [bambu NDR (Novel Discovery Rates)](https://github.com/GoekeLab/bambu) and/or [TransforKmers TSS validation](https://github.com/IGDRion/transforkmers) to assess fulllength transcripts.
8. Perform a quality control of both the full and filtered extended annotations (see [example](https://github.com/igdrion/ANNEXA/blob/master/examples/results/qc_gtf.pdf)).
9. Optional: Check gene body coverage with [RSeQC](http://rseqc.sourceforge.net/#genebody-coverage-py).

This pipeline has been tested with reference annotation from Ensembl and NCBI-RefSeq.

## System requirements
ANNEXA requires 0.5GB of disk space for installation, and requires at least 16 GB of RAM and 2 CPU threads to run the test example.
By default, it requires. 40 GB of RAM and 8 CPU threads. This can be changed in the options (see Performance options).

## Usage

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
(Optional) If you are using Linux, we recommend installing [Apptainer](https://apptainer.org/docs/user/main/quick_start.html) (formerly Singularity) to handle tools used by ANNEXA.

2. Test the pipeline on a small dataset

```sh
nextflow run IGDRion/ANNEXA \
    -profile test
```

Or, if you have installed Apptainer: 

```sh
nextflow run IGDRion/ANNEXA \
    -profile test,singularity
```

3. Run ANNEXA on your own data (change input, gtf, fa with path of your files).

```sh
nextflow run IGDRion/ANNEXA \
    -profile {test,docker,singularity,conda,slurm} \
    --input samples.txt \
    --gtf /path/to/ref.gtf \
    --fa /path/to/ref.fa
```

The input parameter takes a file listing the `bam` path files to analyze (see example below)

```
/path/to/1.bam
/path/to/2.bam
/path/to/3.bam
```

### Options

```sh
Required options
  --input               [string]  Path to file listing paths to bam files.
  --fa                  [string]  Path to reference genome.
  --gtf                 [string]  Path to reference annotation.

Profile options
  --profile test        [string]  Run annexa on toy dataset.
  --profile slurm       [string]  Run annexa on slurm executor.
  --profile singularity [string]  Run annexa in singularity container.
  --profile conda       [string]  Run annexa in conda environment.
  --profile docker      [string]  Run annexa in docker container.

Main options
  --tx_discovery        [string]  Specify which transcriptome reconstruction tool to use. (accepted: bambu, stringtie2) [default: bambu]
  --filter              [boolean] Perform or not the filtering step. [default: true]
  --withGeneCoverage    [boolean] Run RSeQC (can be long depending on annotation and bam sizes). [default: false]

Bambu options
  --bambu_strand        [boolean] Run bambu with stranded data [default: true]
  --bambu_singleexon    [boolean] Include single exon transcripts in Bambu output or not. These are known to have a high frequency of false positives.
                                  [default: true]
  --bambu_threshold     [integer] bambu NDR threshold below which new transcripts are retained. [default: 0.2]
  --bambu_rec_ndr       [boolean] Use NDR threshold recommended by Bambu instead of preset threshold. [default: false]

Filtering options
  --tfkmers_tokenizer   [string]  Path to TransforKmers tokenizer. Required if filter option is activated.
  --tfkmers_model       [string]  Path to TransforKmers model. Required if filter activated.
  --tfkmers_threshold   [integer] TransforKmers prediction threshold below which new transcripts are retained. [default: 0.2]
  --operation           [string]  Operation to retained novel transcripts. 'union' retain tx validated by either bambu or transforkmers, 'intersection' retain
                                  tx validated by both. (accepted: union, intersection) [default: intersection]

Performance options
  --maxCpu              [integer] Max cpu threads used by ANNEXA. [default: 8]
  --maxMemory           [integer] Max memory (in GB) used by ANNEXA. [default: 40]

Nextflow options
  --resume              [null]    Resume task from cached work (useful for recovering from errors when using singularity).
  --with-report         [null]    Create an HTML execution report with metrics such as resource usage for each workflow process.
```

If the filter argument is set to `true`, TransforKmers model and tokenizer paths have to be given. They can be either downloaded from the [TransforKmers official repository](https://github.com/IGDRion/TransforKmers) or trained in advance by yourself on your own data.

### Filtering step

By activating the filtering step (`--filter`), ANNEXA proposes to filter the generated extended annotation according to 2 methods:

1. By using the NDR proposed by [bambu](https://github.com/GoekeLab/bambu). This threshold includes several information such as sequence profile, structure (mono-exonic, etc) and quantification (number of samples, expression). Each transcript with an NDR below the classification threshold will be retained by ANNEXA (default: 0.2).

2. By analysing the Transcription Start Sites (TSS) of each new transcripts using the [TransforKmers](https://github.com/IGDRion/TransforKmers) deep-learning based tool. Each TSS validated below a certain threshold will be retained (default: 0.2). We already provide 2 trained models for filtering TSS with TransforKmers.

- A [human specific model](https://genostack-api-swift.genouest.org/v1/AUTH_07c8a078861e436ba41c4416a821e5d0/transforkmers/hsa_5prime_bert_6-512.zip?temp_url_sig=59e4bd439f42fc2bb8953e78eae82306466917d2&temp_url_expires=2661501621)
- A [dog specific model](https://genostack-api-swift.genouest.org/v1/AUTH_07c8a078861e436ba41c4416a821e5d0/transforkmers/dog_5prime_bert_6-512.zip?temp_url_sig=a5378b6f2cc9ffc10b8f5d4fa6e535070d22f845&temp_url_expires=2661844043)

To use them, extract the zip, and point `--tfkmers_model` and `--tfkmers_tokenizer` to the subdirectories.

The filtered annotation can be the `union` of these 2 tools, _i.e._ all the transcripts validated by one or both of these tools; or the `intersection`, _i.e._ the transcripts validated by both tools (the latter being the default). Please, feee free to see the [dedicated wiki page](https://github.com/IGDRion/ANNEXA/wiki/ANNEXA-wiki#fitlering-operations).

At the end, the QC steps are performed both on the full and filtered extended annotations.