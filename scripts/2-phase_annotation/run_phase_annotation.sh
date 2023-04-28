#! /bin/bash

#FILE=$(basename $1)
FILE="FUT1.vcf"

PATH_TO_VEP="~/ensembl-vep-109/vep"
PATH_TO_PHASE="~/phase.2.1.1.linux/PHASE"

NAME=${FILE%.*}

# PHASING
perl vcf2phase.pl $FILE $NAME.input
$PATH_TO_PHASE -c $NAME.input $NAME.out
perl phase2vcf.pl $NAME.out $FILE $NAME.phased.vcf

# ANNOTATION
$PATH_TO_VEP --custom ../data/futs.sorted.gtf.gz \
			 --fasta ../../data/corrected_Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
			 --force_overwrite --input_file $NAME.phased.vcf \
			 --output_file $NAME.vep.vcf \
			 --vcf
