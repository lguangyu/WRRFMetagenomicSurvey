#!/bin/bash
#SBATCH -pshort -N1 -c1
#SBATCH -o /dev/null

tar -jcvf quast.tar.bz2 \
	quast/
