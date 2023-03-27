# Author: Arey Ferrero Ramos.
# Date: March 24, 2023. Version: 1
# Description: This script calculates the alpha diversity of a the wild and captivity individuals from a species of vertebrate
#        from those used in the study.
#   Parameters:
#       -A first file that contains a table with the number of species per genus in every individual of that species of
#        vertebrate used in that study.
#       -A second file with the metadata of all the samples in the study. We need it for telling wild from captive individuals.
#   Output:
#       -The alpha diversity of the wild individuals in the species of vertebrate.
#       -The alpha diversity of the captive individuals in the species of vertebrate.

import pandas as pd
import math
import sys
import os

if len(sys.argv) != 3:
    print("Error: The number of parameters is incorrect. Two files are needed.")
    exit()

if not os.path.isfile(sys.argv[1]):
    print("Error: The first parameter must be a file corresponding to a species of vertebrate.")
    exit()

if not os.path.isfile(sys.argv[2]):
    print("Error: The second parameter must be a file corresponding to a the metadata of all the samples in the study.")
    exit() 

df_vertebrate = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)

acum_wild = 0
acum_captivity = 0
num_wild = 0
num_captivity = 0

for individual in df_vertebrate:
    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrate[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus
    row = 0
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
        else:
            row += 1
    if sample_type == "Wild":
        acum_wild += num_bacterial_species_per_individual * math.log(num_bacterial_species_per_individual)
        num_wild += 1
    else:
        acum_captivity += num_bacterial_species_per_individual * math.log(num_bacterial_species_per_individual)
        num_captivity += 1

print("Alpha_diversity wild: "+str(round(acum_wild / num_wild))+" species.")
print("Alpha_diversity captivity: "+str(round(acum_captivity / num_captivity))+" species.")