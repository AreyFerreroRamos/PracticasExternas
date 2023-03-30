# Author: Arey Ferrero Ramos.
# Date: March 24, 2023. Version: 1.
# Description: This script calculates the alpha diversity of the wild and captivity individuals from vertebrate species from
#        those used in the study.
#   Parameters:
#       -A first file that contains a table with the number of species per genus in every individual of that vertebrate species
#        used in that study.
#       -A second file with the metadata of all the samples in the study. We need it for telling wild from captive individuals.
#   Output:
#       -The alpha diversity of the wild individuals in vertebrate species.
#       -The alpha diversity of the captive individuals in vertebrate species.
#       -A boxplot that shows the distribution of alpha diversities in both wild and captive individuals in vertebrate species.

from matplotlib import pyplot as plt
from matplotlib import gridspec
import pandas as pd
import math
import sys
import os

if len(sys.argv) != 3:
    print("Error: The number of parameters is incorrect. Two files are needed.")
    exit()

if not os.path.isfile(sys.argv[1]):
    print("Error: The first parameter must be a file corresponding to a vertebrate species.")
    exit()

if not os.path.isfile(sys.argv[2]):
    print("Error: The second parameter must be a file corresponding to a the metadata of all the samples in the study.")
    exit() 

df_vertebrates = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)

alpha_diversities = {'Wild': {}, 'Captivity': {}}

for individual in df_vertebrates: 
    row = 1
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
            specie = df_metadata.loc[row, df_metadata.columns[2]]
        else:
            row += 1

    if specie not in alpha_diversities[sample_type]:
        alpha_diversities[sample_type][specie] = []

    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus
    alpha_diversity = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        if num_bacterial_species_per_genus != 0:
            alpha_diversity += (num_bacterial_species_per_genus / num_bacterial_species_per_individual) * math.log(num_bacterial_species_per_genus / num_bacterial_species_per_individual)

    alpha_diversities[sample_type][specie].append(round(- alpha_diversity, 2))

for sample_type in alpha_diversities:
    print(sample_type)
    for specie in alpha_diversities[sample_type]:
        print(specie, alpha_diversities[sample_type][specie])

plt.boxplot([alpha_diversities['Wild']['MYTR'], alpha_diversities['Captivity']['MYTR']], labels=['Wild', 'Captivity'], patch_artist=True)

plt.title("Bacterial diversity in animal species")
plt.xlabel("Sample type")
plt.ylabel("Alpha diversity")

plt.show()