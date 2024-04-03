#!/bin/bash
#SBATCH -o .log/analysis.kaiju.log
#SBATCH -pshort -N1 -c12
#SBATCH --mem 224G --time 24:00:00

data_dir="../../data"
mapping_dir="mapping"
analysis_dir="analysis"
tmpdir="$analysis_dir/kaiju.tmp"

mkdir -p $analysis_dir

rm -rf $tmpdir/
mkdir $tmpdir

nodes_dmp="$data_dir/kaiju_db/nodes.dmp"
names_dmp="$data_dir/kaiju_db/names.dmp"
kaiju_fmi="$data_dir/kaiju_db/kaiju_db_nr.fmi"

# extract mapped reads
$HOME/opt/htslib/default/bin/samtools fastq \
	-o $tmpdir/kaiju.fastq \
	$mapping_dir/mapped.sorted.rmdup.bam

# run kaiju
$HOME/opt/kaiju/default/bin/kaiju \
	-z $SLURM_CPUS_PER_TASK \
	-t $nodes_dmp \
	-f $kaiju_fmi \
	-i $tmpdir/kaiju.fastq \
	-o $analysis_dir/kaiju.out

# add names
$HOME/opt/kaiju/default/bin/kaiju-addTaxonNames \
	-t $nodes_dmp \
	-n $names_dmp \
	-r superkingdom,phylum,class,order,family,genus,species \
	-i $analysis_dir/kaiju.out \
	-o $analysis_dir/kaiju.names.out

# table summary
for rank in {order,family,genus,species}; do
	$HOME/opt/kaiju/default/bin/kaiju2table \
		-t $nodes_dmp \
		-n $names_dmp \
		-r $rank \
		-o $analysis_dir/kaiju.$rank.tsv \
		$analysis_dir/kaiju.out
done

rm -rf $tmpdir
