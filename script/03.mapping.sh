#!/bin/bash
#SBATCH -pshort -N1 -c12 --time 12:00:00
#SBATCH -o .log/mapping.log

qc_dir="qc"
assembly_dir="assembly"
output_dir="mapping"

mkdir -p $output_dir

. $HOME/.local/env/python-3.10.4-2023-01-27-venv-anvio-7/bin/activate

# build bowtie index
$HOME/opt/bowtie/2.4.5/bin/bowtie2-build \
	--threads $SLURM_CPUS_PER_TASK \
	$assembly_dir/final.contigs.fasta \
	$output_dir/contigs.fasta.bt2_idx

# align and filter
$HOME/opt/bowtie/2.4.5/bin/bowtie2 \
	--sensitive \
	--local \
	-p $SLURM_CPUS_PER_TASK \
	-x $output_dir/contigs.fasta.bt2_idx \
	-1 $qc_dir/R1.fastq \
	-2 $qc_dir/R2.fastq \
| $HOME/opt/htslib/default/bin/samtools view \
	-F 0x4 -b \
	> $output_dir/mapped.bam

deactivate

# sort
$HOME/opt/htslib/default/bin/samtools sort \
	-@ $SLURM_CPUS_PER_TASK \
	$output_dir/mapped.bam \
	> $output_dir/mapped.sorted.bam

# remove duplicates
$HOME/opt/java/jre1.8.0_281/bin/java \
	-jar $HOME/Jar/picard.jar \
	MarkDuplicates \
	REMOVE_DUPLICATES=ture \
	I=$output_dir/mapped.sorted.bam \
	O=$output_dir/mapped.sorted.rmdup.bam \
	M=$output_dir/marked_dup_metrics.txt

# index bam file
$HOME/opt/htslib/default/bin/samtools index \
	-@ $SLURM_CPUS_PER_TASK \
	$output_dir/mapped.sorted.rmdup.bam

# clean up intermediates
rm -f $output_dir/*.bt2_idx*
rm -f $output_dir/mapped.bam
rm -f $output_dir/mapped.sorted.bam

