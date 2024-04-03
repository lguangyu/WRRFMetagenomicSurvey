#!/bin/bash
#SBATCH -pshort -N1 -c4
#SBATCH -o .log/qc.trimmomatic.log

input_dir="fastq"
output_dir="trimmomatic"

mkdir -p $output_dir
$HOME/opt/java/jre1.8.0_281/bin/java \
	-jar $HOME/Jar/Trimmomatic-0.39/trimmomatic-0.39.jar \
	PE $input_dir/R1.fastq $input_dir/R2.fastq \
	$output_dir/R1.fastq \
	$output_dir/R1.unpaired.fastq \
	$output_dir/R2.fastq \
	$output_dir/R2.unpaired.fastq \
	-threads $SLURM_CPUS_PER_TASK \
	ILLUMINACLIP:$HOME/Jar/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:36

ln -sfT $output_dir qc

