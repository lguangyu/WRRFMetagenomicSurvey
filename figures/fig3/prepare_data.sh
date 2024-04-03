#!/bin/bash

mkdir -p prepare_data.tmp

for sample in {DU,HD,MD,SC,UB,WR}; do
	# ppk
	grep 'polyphosphate kinase' \
		../../sample/$sample/gene_calling/prodigal.cds.faa.blastp.add_tax |
	grep -v 'polyphosphate kinase 2' |
	cut -f1 \
	> prepare_data.tmp/$sample.ppk.ids
	# ppx
	grep 'exopolyphosphatase' \
		../../sample/$sample/gene_calling/prodigal.cds.faa.blastp.add_tax |
	cut -f1 \
	> prepare_data.tmp/$sample.ppx.ids

	#
	for g in {ppk,ppx}; do
		for t in {faa,fna}; do
			/opt/seqtk/1.3-2018-06-17/seqtk subseq \
				../../sample/$sample/gene_calling/prodigal.cds.$t \
				prepare_data.tmp/$sample.$g.ids \
				> data/$sample.$g.$t
		done
		cut -f1,10 \
			../../sample/$sample/gene_calling/prodigal.cds.faa.blastp.add_tax |
		# add tax
		table_extract_lines_by_value \
			-f prepare_data.tmp/$sample.$g.ids -k 1 |
		python3 add_tax_name.py \
		> data/$sample.$g.tax
	done
done

for g in {ppk,ppx}; do
	for t in {faa,fna,tax}; do
		cat data/*.$g.$t > data/$g.$t
	done
done

grep 'Accumulibacter' reannotate/ppk.faa.blastp.add_tax |
cut -f1 \
> prepare_data.tmp/acc.ids

/opt/seqtk/1.3-2018-06-17/seqtk subseq \
	data/ppk.fna \
	prepare_data.tmp/acc.ids \
> prepare_data.tmp/acc_ppk.fasta

# ppk tree
cat prepare_data.tmp/acc_ppk.fasta \
	../../supp/ppk1_cap.ref_tree/pplacer.refpkg/ppk1.muscle.clean.fasta \
> data/acc.ppk.with_ref.fasta

/opt/mafft/7.490-2021-10-30/bin/mafft \
	--thread 8 \
	--auto \
	data/acc.ppk.with_ref.fasta \
> data/acc.ppk.with_ref.fasta.mafft.fasta

export LC_ALL=C
/opt/pplacer/1.1.alpha19/pplacer \
	-c ../../supp/ppk1_cap.ref_tree/pplacer.refpkg \
	-p \
	-o data/acc.ppk.fna.mafft.fasta.jplace \
	data/acc.ppk.with_ref.fasta.mafft.fasta

rm -rf prepare_data.tmp/
