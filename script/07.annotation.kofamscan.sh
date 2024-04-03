#!/bin/bash
#SBATCH -pshort -N1 -c16 --time 24:00:00
#SBATCH -o .log/annotation.kofamscan.log

data_dir="data"
data_kofam_db_dir="$data_dir/kofam_db"
gene_calling_dir="gene_calling"
tmpdir="$gene_calling_dir/kofamscan.tmp"

module load ruby/2.5.3

$HOME/opt/kofamscan/default/exec_annotation \
	-p $data_kofam_db_dir/profiles \
	-k $data_kofam_db_dir/ko_list \
	--cpu $SLURM_CPUS_PER_TASK \
	--tmp-dir $tmpdir \
	-f mapper \
	-o $gene_calling_dir/prodigal.cds.faa.kofamscan \
	$gene_calling_dir/prodigal.cds.faa

rm -rf $tmpdir
