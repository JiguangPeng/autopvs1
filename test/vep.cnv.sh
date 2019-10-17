#!/bin/bash

#vep --offline -i clinvar_lof_test.vcf -o output --refseq --use_given_ref --hgvs --hgvsg --canonical --symbol --check_existing
#vep --offline -i clinvar_lof_test.vcf -o output.vcf --vcf --refseq --use_given_ref --hgvs --hgvsg --canonical --symbol --check_existing
#vep --offline -i clinvar_lof_test.vcf -o STDOUT --tab --refseq --use_given_ref --hgvs --hgvsg --canonical --symbol --check_existing
vep --offline --refseq --use_given_ref \
    --species "homo_sapiens" \
    --assembly "GRCh37" \
    --fork 1 \
    --hgvs --hgvsg --canonical --symbol \
    --distance 0 \
    --exclude_predicted \
    --flag_pick \
    --hgvs \
    --hgvsg \
    --lookup_ref \
    --input_file $1 \
    --output_file STDOUT --no_stats \
    --numbers \
    --tab --fields "Uploaded_variation,SYMBOL,Feature,CANONICAL,PICK,EXON,INTRON,Consequence,CDS_position"
#    --everything --tab
