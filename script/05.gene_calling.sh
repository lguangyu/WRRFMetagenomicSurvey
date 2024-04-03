#!/bin/bash
#SBATCH -pshort -N1 -c1
#SBATCH -o .log/gene_calling.log

assembly_dir="assembly"
output_dir="gene_calling"

mkdir -p $output_dir
$HOME/opt/prodigal/default/prodigal \
	-a $output_dir/prodigal.cds.faa \
	-d $output_dir/prodigal.cds.fna \
	-f gff \
	-p meta \
	-q \
	-i $assembly_dir/final.contigs.fasta \
	-o $output_dir/prodigal.gff

