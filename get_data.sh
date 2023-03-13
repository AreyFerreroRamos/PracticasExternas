#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: March 11, 2023. Version: 1.
# Description:.
#	Paramethers:
#		-A compressed file with the extension '.spring' (spring compressor) with the data that is need to be processed.
#	Output:
#		-A file that contains the mimetypes of all the files that are in the compressed file that is given by parameter.

if [ $# -ne 1 ]
then
	echo -e "Error: The number of parameters is incorrect. The only parameter has to be a compressed file." >&2
	exit 1
fi

if [ $(echo $1 | cut -f2 -d'.') != 'spring' ]
then
	echo -e "Error: The input must be a compressed file with the extension '.spring' (spring compressor)." >&2
	exit 2
fi

dir_name=$(echo $1 | cut -f1 -d'.')
mkdir $dir_name
tar -xf $1 -C $dir_name
IFS=$'\n'
for file_name in $(ls -l $dir_name | tail -n +2 | tr -s ' ' | cut -f9 -d' ')
do
	file --mime-type $dir_name/$file_name >> data.txt
done
exit 0
