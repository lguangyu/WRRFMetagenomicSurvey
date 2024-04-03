#!/bin/bash
#SBATCH -pshort -N1 -c12 --mem 280GB --time 24:00:00
#SBATCH -o .log/assembly.spades.log

input_dir="qc"
output_dir="spades"
tmp_dir="$output_dir/tmp"

. $HOME/.local/env/python-3.10.4-2023-01-27-venv-anvio-7/bin/activate

if [[ -d $output_dir ]]; then
	$HOME/opt/spades/3.15.5/bin/spades.py \
		--continue \
		-o $output_dir
else
	$HOME/opt/spades/3.15.5/bin/spades.py \
		--meta \
		-1 $input_dir/R1.fastq \
		-2 $input_dir/R2.fastq \
		-t $SLURM_CPUS_PER_TASK \
		-m 280 \
		-k 33,55,77,99 \
		-o $output_dir
fi

# reformat fasta to make anvio shutup
anvi-script-reformat-fasta \
	-l 500 \
	--simplify-names --prefix "$(basename $(pwd))" \
	-o $output_dir/final.contigs.fasta \
	$output_dir/contigs.fasta 

deactivate

ln -sfT $output_dir assembly
rm -rf $tmp_dir

