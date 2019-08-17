#!/bin/bash

# File of SNP with methylayion % in REF and ALT and TOTAL (4 columns)
touch ${SAMPLE}_geno.txt 

# Create an empty file to receive SAM-format for each allele
touch ${SAMPLE}_ref.sam ${SAMPLE}_alt.sam 

# Need SNP pos and the number of reads