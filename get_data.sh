#!/bin/bash

# Author: Arey Ferrero Ramos.
# Data: March 13, 2023. Version: 1.
# Description:.
#	Parameters:
#		-
#	Output:
#		-

if [ $# -ne 1 ]
then
	echo -e "Error: The number of parameters is incorrect. The only parameter has to be " >&2
	exit 1
fi

if [ $(echo $1 | cut -f2 -d'.') != 'json' ]
then
	echo -e "The input must be a " >&2
	exit 2
fi


exit 0
