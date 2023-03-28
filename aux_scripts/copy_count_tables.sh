#!/bin/bash

# Author: Arey Ferrero Ramos.
# Data: March 28, 2023. Version: 1.
# Description: Given a directory, this script copys all the 'count_genus' files of that directory in a local directory inside the project."
# 	Input:
#		-A directory with a list of files.
#	Output:
#		-A directory with all the 'count_genus' files from the input directory.

if [ $# -eq 0 ]
then
	echo -e "copy_count_tables.sh: $0 directory" >&2
	exit 1
fi

if [ $# -gt 1 ]
then
	echo -e "Error: The number of parameters is incorrect. The input must be only one directory with a list of files." >&2
	exit 2
fi

if [ ! -d "$1" ]
then
	echo -e "Error: The input is not a directory. It must be a directory with a list of files." >&2
	exit 3
fi

if [ -z "$(ls -A -- "$1")" ]
then
	echo -e "Error: The input directory is empty. The directory must contains a list of files." >&2
	exit 4
fi

mkdir -p ../Tables

for file_name in "$1"/*
do
	if [ ! -z $(echo "$file_name" | grep -E ".*count_genus.*") ]
	then
		cp -p "$file_name" ../Tables 
	fi
done

exit 0
