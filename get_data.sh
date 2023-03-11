#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: March 11, 2023. Version: 1.
# Description:.
#	Paramethers:
#		- Directory with the data is needed to process.
#	Output:
#		- 

IFS=$'\n'
for file in $(ls -l $1 | tail -n +2 | tr -s ' ' | cut -f9 -d' ')
do
	print $file
done
