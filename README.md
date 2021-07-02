# ANNEXA: Analysis of Nanopore with Nextflow for EXtended Annotation

## Introduction

**ANNEXA** is an all-in-one reproductible pipeline, written in the [Nextflow](https://nextflow.io), which allows users to analyze LR-RNAseq sequences from Oxford Nanopore Technologies (ONT), and to reconstruct and quantify known and novel genes and isoforms. 

More specifically, ANNEXA works by using only three parameter files (a reference genome, a reference annotation and mapping files) and provides users with an extended annotation distinguishing between novel protein-coding (mRNA) versus long non-coding RNAs (lncRNA) genes. All known and novel gene/transcript models are further characterized through multiple features (length, number of spliced transcripts, normalized expression levels,...) available as graphical outputs.

## Usage

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Test the pipeline on a small dataset

```sh
nextflow run mlorthiois/ANNEXA \
    -profile test
```

3. Run ANNEXA on your own data (change input, gtf, fa with path of your files).

```sh
nextflow run mlorthiois/ANNEXA \
    --input samples.txt \
    --gtf /path/to/ref.gtf \
    --fa /path/to/ref.fa
```

The input parameter takes a file listing the bams to analyze (see example below)

```
/path/to/1.bam
/path/to/2.bam
/path/to/3.bam
```

### Options

```
--profile test    : Run annexa on toy dataset
--profile slurm   : Run annexa on slurm executor
--input           : Path to file listing paths to bam files
--fa              : Path to reference genome
--gtf             : Path to reference annotation
--withGeneCoverage: Run RSeQC (can be long depending on annotation and bam sizes). False by default
```

## Pipeline summary

This pipeline has been tested with reference annotation from Ensembl and NCBI-RefSeq.

1. Transcriptome reconstruction and quantification with [bambu](https://github.com/GoekeLab/bambu)
2. Novel classification with [FEELnc](https://github.com/tderrien/FEELnc)
3. Retrieve information from input annotation and format final gtf with 3 level structure: gene -> transcript -> exon.
4. Perform a quality control of the extended annotation (see [example](https://github.com/mlorthiois/ANNEXA/blob/master/examples/results/qc_gtf.pdf)).
5. Check gene body coverage with [RSeQC](http://rseqc.sourceforge.net/#genebody-coverage-py)

## Filtering extended annotation

At the end of the pipeline, depending on the quality report, you may want to filter some new genes. For this, ANNEXA provides a script that allows you to filter the annotation according to 2 criteria: 

- The **structure of the genes** (mono-isoform and/or mono-exonic): The script allows you to filter the annotation according to whether a new gene is mono-isoform or mono-exonic. For this, the `--filterStruct` option takes as argument `Y` or `N` for each structure (isoform then exon) if you want to filter or not the structure, and the `|` (or) or `&` (and) operator.

  Let's say you want to remove new genes that are mono-isoform AND mono-exon, use `--filterStruct "Y&Y"`. If you want to filter the genes that are mono-isoform OR mono-exonic, use `--filterStruct "Y|Y"`. If now you want to filter only mono-isoforms (mono-exonics), use `--filterStruct "Y&N"` (`--filterStruct "N&Y"` respectively). 

- The **quantification aspect of the genes** : The script also allows you to filter the annotation according to a quantitative criterion. In the same way as for structures, the `--filterQuant` option takes as argument a number or `NA` for the minimum number of reads to keep the gene, and the number of samples expressing that gene.

  Let's say you want to remove new genes validated by less than 50 reads AND present in less than 5 samples, use `--filterQuant "50&7"`. If you want to filter out genes validated by less than 50 reads OR present in less than 5 samples, use `--filterQuant "50|5"`. If you now only want to filter genes validated by less than 50 reads (or only in less than 5 samples), use `--filterQuant "50&NA"` (or `--filterQuant "NA&5"` respectively)


### Use case

1. This script needs `pandas` to work. Pandas is included in the conda environment of ANNEXA, I first need to activate this environment or install python:

```sh
# Method 1: Activate conda environment of ANNEXA
conda activate work/conda/*/

# Method 2: Install pandas
pip install pandas
```

2. Now, I want to remove mono-exonic genes and genes that are not validates by at least 30 reads and not present in 6 samples. I use the script `filter_gtf.py` stored in the bin subfolder of ANNEXA :

```sh
./bin/filter_gtf.py \
  --gtf results/final/extented_annotations.gtf \
  --gene_stats results/qc/gene.stats \
  --tx_stats results/qc/transcript.stats \
  --filterQuant="30&6" \
  --filterStruct "N&Y"
```

Help page:

```
usage: filter_gtf.py [-h] --gtf GTF --gene_stats GENE_STATS --tx_stats TX_STATS
                     [--filterStruct {Y&Y,Y|Y,Y&N,Y|N,N&Y,N|Y,N|N}] [--filterQuant FILTERQUANT] [-o OUTPUT]

Filter your extended GTF generated by ANNEXA. See https://github.com/mlorthiois/annexa for more informations

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Path to extended_annotations.gtf generated by ANNEXA
  --gene_stats GENE_STATS
                        Path to gene.stats generated by ANNEXA
  --tx_stats TX_STATS   Path to transcript.stats generated by ANNEXA
  --filterStruct {Y&Y,Y|Y,Y&N,Y|N,N&Y,N|Y,N|N}
                        Filter or not mono-isoform and/or mono-exonic genes.
  --filterQuant FILTERQUANT
                        Filter or not genes with less than X read counts and/or present in less than X examples.
  -o OUTPUT, --output OUTPUT
                        Directs the filtered gtf to file with the name of your choice
```

## Troubleshooting

GTF parsing can be hard, and ANNEXA need some informations from the input annotation. The first step is to check if the input annotation contains all the information needed.

For example, your GTF should contains all the 3 levels gene -> transcript -> exon, with the attributes gene_id and transcript correctly annotated. Your gene line also have to contains a gene_biotype attributes.

