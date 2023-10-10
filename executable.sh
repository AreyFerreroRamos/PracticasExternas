#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: March 31, 2023. Version: 3.
# Description: This script runs the command that runs the project.
#	Input:
#		-The type of data structure that we need
#		-The type of plot that we want to display.
#	Output:
#		-A plot that show the results

if [ ! -d input_files ]
then
	echo -e "Error: The directory input_files that contains the input data for the project doesn't exist." >&2
	exit 1
fi

if [ ! -f input_files/count_Genus_all.tsv ]
then
	echo -e "Error: The file count_Genus_all.tsv, that contains the diversity of the individuals of all the vertebrate species in the study, doesn't exist." >&2
	exit 2
fi

if [ ! -f input_files/metadata.csv ]
then
	echo -e "Error: The file metadata.csv, that contains the metadata of all the samples in the study, doesn't exist." >&2
	exit 3
fi

if [ ! -f input_files/sp_code.txt ]
then
	echo -e "Error: The file sp_code.txt, that contains a list of the scientific name that corresponds to the code that identifies a vertebrate species, doesn't exist" >&2
	exit 4
fi

python3 main.py input_files/count_Genus_all.tsv input_files/metadata.csv input_files/sp_code.txt "$1" "$2" "$3"

exit 0
