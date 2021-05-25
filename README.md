# ANNEXA: Analysis of Nanopore with Nextflow for EXtended Annotation

[Nextflow](https://www.nextflow.io) pipeline to extend reference annotation with nanopore reads, classify novel genes (mRNAs vs lncRNAs), and count gene/transcript expression levels.

## Usage

Just use the following command (change gtf, fa with your files).

```sh
nextflow run mlorthiois/ANNEXA \
    --input input.csv \
    --gtf /path/to/ref.gtf \
    --fa /path/to/ref.fa \
    --sampleNumber 1 \
    --readCount 5
```

Here, input parameter takes file with BAM files and condition (separated by `,`) (See example below)

```
sample,condition
/path/to/1.bam,MM
/path/to/2.bam,MM
/path/to/3.bam,HS
/path/to/4.bam,HS
```

### Options

```
--profile test    : Run annexa on toy dataset
--profile slurm   : Run annexa on slurm executor
--readCount n     : Keep novel transcripts if min n read counts. 5 by default
--sampleNumber k  : Keep novel tanscripts with min k samples having n counts. 1 by default
--input           : Path to samplesheet
--fa              : Path to reference genome
--gtf             : Path to reference annotation
--withGeneCoverage: Run RSeQC (can be long depending on annotation and bam sizes). False by default
```

sampleNumber and readCount are params passed to [bambu](https://github.com/GoekeLab/bambu#advanced-options).

### Test the pipeline

A toy dataset is provided to test the pipeline. Data come from real data. To use it, just set the profile as test like in the command below :

```sh
nextflow run mlorthiois/ANNEXA \
    --profile test
```

## Steps

1. Bambu. Ouputs: extended annotation, gene and transcript counts
2. FEELnc. Outputs: Novel genes identified by bambu classified in protein coding or lncRNA
3. Extended annotation with novel transcriptsg/gene classified as protein coding/lncRNA
4. DESeq2 normalization of gene counts
5. QC on the final gtf with normalized gene counts
6. Gene body coverage with RSeQC
