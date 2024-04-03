#!/bin/bash
#SBATCH -pshort -N1 -c1
#SBATCH -o /dev/null

tar -jcvf mapping.tar.bz2 \
	mapping/
