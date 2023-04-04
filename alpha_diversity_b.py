# Author: Arey Ferrero Ramos.
# Date: March 24, 2023. Version: 4.
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
#       -A boxplot that shows the distribution of alpha diversities in both wild and captive individuals in all vertebrate species.

import pandas as pd
import math
import sys
import os
import support_functions
import show_functions

if len(sys.argv) != 4:
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

alpha_diversities_individual = {}
alpha_diversities_specie = {}

previous_specie = ""

for individual in df_vertebrates: 
    row = 1
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            specie = df_metadata.loc[row, df_metadata.columns[2]]
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
        else:
            row += 1

    if specie not in alpha_diversities_individual:
        alpha_diversities_individual[specie] = {'Wild': [], 'Captivity': []}
        alpha_diversities_specie = {'Wild': [], 'Captivity': []}
        num_bacterial_species_per_specie = {'Wild': 0, 'Captivity': 0}
        if previous_specie != specie:
            if previous_specie:
                # Recorregut per alpha_diversities_specie.
            previous_specie = specie
        
        #show_functions.print_specie(specie, f_codes_vertebrates)

    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus
    num_bacterial_species_per_specie[sample_type] += num_bacterial_species_per_individual
    
    alpha_diversity_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        if num_bacterial_species_per_genus != 0:
            alpha_diversity_individual += (num_bacterial_species_per_genus / num_bacterial_species_per_individual) * math.log(num_bacterial_species_per_genus / num_bacterial_species_per_individual)

    alpha_diversities_individual[specie][sample_type].append(round(0 - alpha_diversity_individual, 2))

support_functions.fill_table(alpha_diversities_specie, alpha_diversities_individual, 'Wild')
support_functions.fill_table(alpha_diversities_specie, alpha_diversities_individual, 'Captivity')

show_functions.print_table(alpha_diversities_individual, 'Wild')
show_functions.print_table(alpha_diversities_individual, 'Captivity')

#show_functions.print_table(alpha_diversities_specie, 'Wild')
#show_functions.print_table(alpha_diversities_specie, 'Captivity')

#show_functions.show_plot('Boxplot', alpha_diversities_individual, '')
#show_functions.show_plot('Histogram', alpha_diversities_specie, 'Wild')
#show_functions.show_plot('Histogram', alpha_diversities_specie, 'Captivity')
