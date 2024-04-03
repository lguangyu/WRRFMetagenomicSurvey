#!/bin/bash
#SBATCH -pshort -N1 -c1
#SBATCH -o /dev/null

tar -jcvf trimmomatic.tar.bz2 \
	trimmomatic/
