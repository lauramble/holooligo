#! /bin/bash

VCF="FUT1.vcf"
POPULATION="population.tsv"

NAME=${VCF%.*}

Rscript preprocessing.R -i $VCF -p $POPULATION
Rscript PCA.R -i $VCF -p $POPULATION -t 0.8
Rscript structure.R -i $PWD/$NAME.recode.noheader.strct_in -p $POPULATION
Rscript haplotypes.R -i $VCF -g $PWD/$NAME.gds -p $POPULATION
