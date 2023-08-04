#!/bin/bash

nextflow run main.nf \
    --input examples/input.csv \
    --gtf /path/to/ref.gtf \
    --fa /path/to/ref.fa
