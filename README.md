# ANNEXA: Analysis of Nanopore with Nextflow for EXtended Annotation

[Nextflow](https://www.nextflow.io) pipeline to extend reference annotation with nanopore reads, classify novel genes (mRNAs vs lncRNAs), and count gene/transcript expression levels.

## Usage

Just use the following command (change gtf, fa with your files).
sampleNumber (1 by default) and readCount (5 by default) are params passed to [bambu](https://github.com/GoekeLab/bambu#advanced-options).

```sh
nextflow run mlorthiois/ANNEXA \
    --input input.csv \
    --gtf /path/to/ref.gtf \
    --fa /path/to/ref.fa \
    --sampleNumber 1 \
    --readCount 5
```

This pipeline provide a profile to be runned on slurm executor. To use it, add `--profile slurm` in your command.

Here, input parameter takes file with BAM files and condition (separated by `,`) (See example below)

```
sample,condition
/path/to/1.bam,MM
/path/to/2.bam,MM
/path/to/3.bam,HS
/path/to/4.bam,HS
```

## Test

A toy dataset is provided to test the pipeline. Data come from real data, but only reads mapped to the chr10 are selected.
To use it, just set the profile as test like in the command below :

```sh
nextflow run mlorthiois/ANNEXA \
    --profile test
```

## Steps

1. Bambu. Ouputs: extended annotation, gene and transcript counts
2. FEELnc. Outputs: Novel genes identified by bambu classified in protein coding or lncRNA
3. Final gtf
4. DESeq2 normalization of gene counts
5. QC on the final gtf with normalized gene counts

