#!/bin/bash
#SBATCH -pshort -N1 -c1
#SBATCH -o /dev/null

tar -jcvf gene_calling.tar.bz2 \
	gene_calling/
