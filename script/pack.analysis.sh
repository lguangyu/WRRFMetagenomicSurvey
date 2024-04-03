#!/bin/bash
#SBATCH -pshort -N1 -c1
#SBATCH -o /dev/null

tar -jcvf analysis.tar.bz2 \
	analysis/
