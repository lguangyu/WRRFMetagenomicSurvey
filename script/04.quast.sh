#!/bin/bash
#SBATCH -pshort -N1 -c4
#SBATCH -o .log/quast.log

assembly_dir="assembly"
mapping_dir="mapping"
output_dir="quast"

. $HOME/.local/env/python-3.10.4-2023-01-27-venv-anvio-7/bin/activate

quast.py \
	-o $output_dir \
	-t $SLURM_CPUS_PER_TASK \
	--ref-bam $mapping_dir/mapped.sorted.rmdup.bam \
	$assembly_dir/final.contigs.fasta

deactivate

