#!/bin/bash
#SBATCH -pshort -N1 -c16 --mem 128GB --time 24:00:00
#SBATCH -o .log/annotation.blastp.log

data_dir="data"
data_ncbi_tax_dump_dir="$data_dir/ncbi_tax_dump"
data_blast_db_dir="$data_dir/blast_db"
assembly_dir="assembly"
gene_calling_dir="gene_calling"
tmpdir="$gene_calling_dir/blstp.tmp"
mkdir -p $tmpdir

# done in gene calling dir
$HOME/opt/diamond/diamond-default blastp \
	-q $gene_calling_dir/prodigal.cds.faa \
	-d $data_blast_db_dir/nr \
	-o $gene_calling_dir/prodigal.cds.faa.blastp \
	-f 6 qseqid sseqid pident length mismatch evalue bitscore \
	--max-target-seqs 1 \
	-p $SLURM_CPUS_PER_TASK \
	--tmpdir $tmpdir

# add taxonomy info
accs_list=$tmpdir/accs_list
tax_table=$tmpdir/tax_table
contig_gene_taxid_map=$tmpdir/contig_gene_taxid_map
cut -f 2 $gene_calling_dir/prodigal.cds.faa.blastp > $accs_list

$HOME/opt/ncbi/blast+-2.11.0/bin/blastdbcmd \
	-db $data_blast_db_dir/nr \
	-entry_batch $accs_list \
	-out $tax_table \
	-outfmt "%a	%t	%S	%T"

# add taxonomy to blast output
. $HOME/.local/env/python-3.11-venv-generic/bin/activate

script/blastp_table_add_tax.py \
	-b $gene_calling_dir/prodigal.cds.faa.blastp \
	-t $tax_table \
	-o $gene_calling_dir/prodigal.cds.faa.blastp.add_tax
cut -f 1,10 $gene_calling_dir/prodigal.cds.faa.blastp.add_tax > $contig_gene_taxid_map
./script/contig_tax_lca.py \
	-c $contig_gene_taxid_map \
	-n $data_ncbi_tax_dump_dir/nodes.dmp \
	-a $data_ncbi_tax_dump_dir/names.dmp \
	-b 0.5 \
	-o $gene_calling_dir/prodigal.cds.faa.blastp.add_tax.lca \
	-t $gene_calling_dir/prodigal.cds.faa.blastp.add_tax.lca.newick
./script/count_contig_group_total_length.py \
	-f $assembly_dir/final.contigs.fasta \
	-g $gene_calling_dir/prodigal.cds.faa.blastp.add_tax.lca \
	-o $gene_calling_dir/prodigal.cds.faa.blastp.add_tax.lca.tax_total_size

deactivate
rm -rf $tmpdir

