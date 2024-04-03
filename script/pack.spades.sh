#!/bin/bash
#SBATCH -pshort -N1 -c1
#SBATCH -o /dev/null

tar -jcvf spades.tar.bz2 \
	spades/
