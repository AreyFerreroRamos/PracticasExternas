#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: March 23, 2023. Version: 1.
# Description:
#	Paramethers:
#		-
#	Output:
#		-

if [ $# -ne 1 ]
then
	echo -e "Error: The number of parameters is incorrect. The only parameter has to be a directory." >&2
	exit 1
fi

if [ ! -d "$1" ]
then
	echo -e "Error: The input must be a directory corresponding to a vertebrate species that contains the files of every individual of the specie used in this study." >&2
	exit 2
fi

IFS=$'\n'
for file_name in "$1"/*
do
	alpha_diversity=0
	for line in $(cat $file_name | tail -n +2)
	do
		abundance=$(echo $line | cut -f2 -d$'\t')
		let alpha_diversity=$alpha_diversity+$abundance*$(bc -l <<< "scale=0; l($abundance)")
	done
	let alpha_diversity=-alpha_diversity
	echo $alpha_diversity
done

exit 0
