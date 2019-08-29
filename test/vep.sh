#!/bin/bash

vep --offline --refseq --use_given_ref \
    --species "homo_sapiens" \
    --assembly "GRCh37" \
    --fork 4 \
    --canonical \
    --flag_pick \
    --hgvs --hgvsg --symbol \
    --distance 500 \
    --exclude_predicted \
    --numbers \
    --lookup_ref \
    --input_file $1 \
    --output_file STDOUT --no_stats \
    --tab --fields "Uploaded_variation,SYMBOL,Feature,CANONICAL,PICK,Consequence,HGVSc,HGVSp,HGVSg,EXON,INTRON"
