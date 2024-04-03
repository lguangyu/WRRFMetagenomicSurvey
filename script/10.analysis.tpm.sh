#!/bin/bash
#SBATCH -pshort -N1 -c1 --time 24:00:00
#SBATCH -o .log/analysis.rpkm.log

gene_calling_dir="gene_calling"
mapping_dir="mapping"
analysis_dir="analysis"
tmpdir="$analysis_dir/rpkm.tmp"

mkdir -p $analysis_dir
mkdir -p $tmpdir

# bam index stat
$HOME/opt/htslib/default/bin/samtools view \
	-b -F 0x900 \
	$mapping_dir/mapped.sorted.rmdup.bam \
> $tmpdir/tmp.bam
$HOME/opt/htslib/default/bin/samtools index \
	$tmpdir/tmp.bam
$HOME/opt/htslib/default/bin/samtools idxstats \
	$tmpdir/tmp.bam \
> $tmpdir/tmp.bam.idxstats

# cds coverage
cut -f 1,4,5 $gene_calling_dir/prodigal.gff |
	grep -v '^#' \
> $analysis_dir/cds.bed
$HOME/opt/htslib/default/bin/samtools bedcov \
	-c \
	$analysis_dir/cds.bed \
	$tmpdir/tmp.bam |
$HOME/.local/env/python-3.11-venv-generic/bin/python3 \
	script/rename_bed_contig.py \
	-o $analysis_dir/coverage.tsv

# bam stats
$HOME/opt/htslib/default/bin/samtools stats \
	$tmpdir/tmp.bam |
grep -E '^SN\s+sequences' |
grep -oE '[0-9]+$' \
> $analysis_dir/total_mapped_reads.txt

# per-cds
$HOME/.local/env/python-3.11-venv-generic/bin/python3 \
	script/cds_tpm.py \
	-c $(cat $analysis_dir/total_mapped_reads.txt) \
	-o $analysis_dir/cds_tpm.tsv \
	$analysis_dir/coverage.tsv

# per-module
$HOME/.local/env/python-3.11-venv-generic/bin/python3 \
	./script/ko_rpkm.py \
	-c $analysis_dir/cds_tpm.tsv \
	-k $gene_calling_dir/prodigal.cds.faa.kofamscan \
	-o $analysis_dir/ko_tpm.tsv

$HOME/.local/env/python-3.11-venv-generic/bin/python3 \
	./script/module_rpkm.py \
	-k $analysis_dir/ko_tpm.tsv \
	-m ../../supp/kegg.module_ko.json \
> $analysis_dir/module_tpm.tsv

# per-contig rpkm
$HOME/.local/env/python-3.11-venv-generic/bin/python3 \
	script/contig_rpkm.py \
	-c $(cat $analysis_dir/total_mapped_reads.txt) \
	-o $analysis_dir/contig_rpkm.tsv \
	$tmpdir/tmp.bam.idxstats


rm -rf $tmpdir
