#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: March 11, 2023. Version: 1.
# Description: Given a compressed file with the extension '.spring', this file is descompressed in a directory with the same name as the file (without the extension '.spring') and then, we put on an output file
#		called 'data.txt' all the mimetypes of the files that were contained in the .spring file.
#	Paramethers:
#		-A compressed file with the extension '.spring' (spring compressor) with the data that is need to be processed.
#	Output:
#		-A file that contains the mimetypes of all the files that are in the compressed file that is given by parameter.

if [ $# -ne 1 ]
then
	echo -e "Error: The number of parameters is incorrect. The only parameter has to be a compressed file." >&2
	exit 1
fi

if [ ! -f "$1" ] || [ $(echo "$1" | cut -f2 -d'.') != 'spring' ]
then
	echo -e "Error: The input must be a compressed file with the extension '.spring' (spring compressor)." >&2
	exit 2
fi

DIR="$(dirname "$(realpath "$0")")"

dir_name="$(echo "$1" | cut -f1 -d'.')"
mkdir "$dir_name"
tar -xf "$1" -C "$dir_name"
IFS=$'\n'
for file_name in "$DIR"/"$dir_name"/*
do
	file --mime-type "$file_name" >> data.txt
done
exit 0
