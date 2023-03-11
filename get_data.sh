#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: March 11, 2023. Version: 1.
# Description:.
#	Paramethers:
#		- File with the data that is need to process.
#	Output:
#		-

if [ $# -eq 1 ]
then
	if [ $(echo $1 | cut -f2 -d'.') = 'spring' ]
	then
		name_dir=$(echo $1 | cut -f1 -d'.')
		if [ ! -d $name_dir ]
		then
			mkdir $name_dir
			tar -xvf $1 -C $name_dir
		fi
		IFS=$'\n'
		echo $name_dir
		for file in $(ls -l $name_dir | tail -n +2 | tr -s ' ' | cut -f9 -d' ')
		do
			echo $file
			#file --mime-type $file		# La comanda no funciona encara.
		done
	else
		echo -e "Error: The file that contains the data must be created with the spring compressor." >&2
		exit 1
	fi
	exit 0
else
	echo -e "Error: The number of parameters is incorrect. A compressed file is need to be given as parameter." >&2
	exit 2
fi
