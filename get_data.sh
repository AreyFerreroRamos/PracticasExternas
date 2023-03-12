#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: March 11, 2023. Version: 1.
# Description:.
#	Paramethers:
#		- Compressed file with the data that is need to be processed.
#	Output:
#		-

if [ $# -eq 1 ]
then
	if [ $(echo $1 | cut -f2 -d'.') = 'spring' ]
	then
		dir_name=$(echo $1 | cut -f1 -d'.')
		if [ ! -d $dir_name ]
		then
			mkdir $dir_name
			tar -xf $1 -C $dir_name
		fi
		IFS=$'\n'
		for file_name in $(ls -l $dir_name | tail -n +2 | tr -s ' ' | cut -f9 -d' ')
		do
			file --mime-type $dir_name/$file_name
		done
		exit 0
	else
		echo -e "Error: The input must be a compressed file with the extension '.spring' (spring compressor)." >&2
		exit 1
	fi
else
	echo -e "Error: The number of parameters is incorrect. A compressed file is need to be given as parameter." >&2
	exit 2
fi
