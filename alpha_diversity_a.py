# Author: Arey Ferrero Ramos.
# Date: March 24, 2023. Version: 3.
# Description: This script calculates the alpha diversity of the wild and captivity individuals from vertebrate species from
#        those used in the study.
#   Parameters:
#       -A first file that contains a table with the number of species per genus in every individual of that vertebrate species
#        used in that study.
#       -A second file with the metadata of all the samples in the study. We need it for telling wild from captive individuals.
#       -A third file with the scientific name that corresponds to the code that identifies a vertebrate species.
#   Output:
#       -The alpha diversity of the wild individuals in all vertebrate species.
#       -The alpha diversity of the captive individuals in all vertebrate species.
#       -A boxplot that shows the distribution of alpha diversities in both wild and captive individuals in vertebrate species.

from matplotlib import pyplot as plt
from matplotlib import gridspec
import pandas as pd
import math
import sys
import os

if len(sys.argv) != 5:
    print("Error: The number of parameters is incorrect. Three files are needed.")
    exit()

if not os.path.isfile(sys.argv[1]):
    print("Error: The first parameter must be a file corresponding to all vertebrate species.")
    exit()

if not os.path.isfile(sys.argv[2]):
    print("Error: The second parameter must be a file corresponding to the metadata of all the samples in the study.")
    exit()

if not os.path.isfile(sys.argv[3]):
    print("Error: The third parameter must be a file with a list of the scientific name that corresponds to the code that identifies a vertebrate species.")
    exit()

df_vertebrates = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)
f_codes_vertebrates = open(sys.argv[3], 'r')

alpha_diversities = {}

for individual in df_vertebrates: 
    row = 1
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            specie = df_metadata.loc[row, df_metadata.columns[2]]
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
        else:
            row += 1

    if specie not in alpha_diversities:
        alpha_diversities[specie] = {'Wild': [], 'Captivity': []}

    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus
    alpha_diversity = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        if num_bacterial_species_per_genus != 0:
            alpha_diversity += (num_bacterial_species_per_genus / num_bacterial_species_per_individual) * math.log(num_bacterial_species_per_genus / num_bacterial_species_per_individual)

    alpha_diversities[specie][sample_type].append(round(0 - alpha_diversity, 2))

print('Wild')
for specie in alpha_diversities:
    print(specie, alpha_diversities[specie]['Wild'])
print('Captivity')
for specie in alpha_diversities:
    print(specie, alpha_diversities[specie]['Captivity'])

for vertebrate_specie in f_codes_vertebrates:
    if sys.argv[4] == vertebrate_specie.split()[0]:
        name_specie = vertebrate_specie.split()[1]

figure = plt.figure()
spec = gridspec.GridSpec(nrows=5, ncols=5, figure=figure)

ax_list = []
row = 0
for specie in alpha_diversities:
    column = 0
    while column < 5:
        ax_box = figure.add_subplot(spec[row, column])
        ax_box.boxplot([alpha_diversities[sys.argv[4]]['Wild'], alpha_diversities[sys.argv[4]]['Captivity']], labels=['Wild', 'Captivity'])
        #ax_box.title(str(name_specie.replace('_', ' ', 1)))
        #ax_box.xlabel("Sample type")
        #ax_box.ylabel("Alpha diversity")
        column += 1
    row += 1

#plt.boxplot([alpha_diversities[sys.argv[4]]['Wild'], alpha_diversities[sys.argv[4]]['Captivity']], labels=['Wild', 'Captivity'])

#plt.title(name_specie.replace('_', ' ', 1))
#plt.xlabel("Sample type")
#plt.ylabel("Alpha diversity")

plt.suptitle("Bacterial diversity in vertebrate species")
plt.show()