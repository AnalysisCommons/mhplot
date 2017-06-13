#!/bin/bash

p_col="p"
gene_col="gene"
chr_col="chr"
pos_col="pos"
snp_col="Name"
chr_labels="NULL"
offset="0.01"
suggestive_line="-log10(1e-5)"
genomewide_line="-log10(5e-8)"
highlight="NULL"
logp="TRUE"
main_title="NULL"
seqmetaresults=( /home/mbrown/DNAnexus/mb_test_data_nochrpos.Rda )
mh_color=( gray50 black )
