#!/bin/bash

# Author: Arey Ferrero Ramos.
# Data: March 28, 2023. Version: 2.
# Description: Given a directory, this script copys all the files from that directory which have the information in wich we are interested (specified by parameter) in a local directory inside our project.
# 	Input:
#		-A directory with a list of files.
#		-The type of information in which we are interested.
#	Output:
#		-A directory with all the files from the input directory that contains the information in which we are interested.

if [ $# -eq 0 ]
then
	echo -e "copy_count_tables.sh: $0 directory search_information" >&2
	exit 1
fi

if [ $# -gt 2 ]
then
	echo -e "Error: The number of parameters is incorrect. The input must be a directory with a list of files and the information that indicated in which files we are interested in." >&2
	exit 2
fi

if [ ! -d "$1" ]
then
	echo -e "Error: The first input is not a directory. It must be a directory with a list of files." >&2
	exit 3
fi

if [ -z "$(ls -A -- "$1")" ]
then
	echo -e "Error: The first input is a directory, but it is empty. The directory must contains a list of files." >&2
	exit 4
fi

if [ -z "$2" ]
then
	echo -e "Error: The second input is empty. It must be a string with the information that indicated in which files of the directory we are interested in." >&2
	exit 5
fi

mkdir -p ../Tables/$2

for file_name in "$1"/*
do
	if [ ! -z $(echo "$file_name" | grep -E ".*$2.*") ]
	then
		cp -p "$file_name" ../Tables/$2
	fi
done

exit 0
