# Author: Arey Ferrero Ramos.
# Date: March 24, 2023. Version: 5.
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

import support_functions as spf
import show_functions as shf
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

alpha_diversities_individual = {}
alpha_diversities_specie = {}
relative_abundances = {}

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
        alpha_diversities_specie[specie] = {}
        
        #shf.print_specie(specie, f_codes_vertebrates)
    
    relative_abundances[individual] = []
    num_bacterial_species_per_individual = alpha_diversity = num_zeros = num_genus = 0
    
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus
    
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        if num_bacterial_species_per_genus != 0:
            relative_abundances[individual].append(num_bacterial_species_per_genus / num_bacterial_species_per_individual)
            alpha_diversity += (num_bacterial_species_per_genus / num_bacterial_species_per_individual) * math.log(num_bacterial_species_per_genus / num_bacterial_species_per_individual)
        else:
            num_zeros += 1
        num_genus += 1
    
    alpha_diversities_individual[specie][sample_type].append(round(0 - alpha_diversity, 4))
    
    if (individual == sys.argv[4]):
        print(individual+" ("+specie+", "+sample_type+")\t"+str(round(num_zeros / num_genus * 100, 2))+"% zeros.")
        shf.show_histogram(relative_abundances, individual)

#spf.calculate_alpha_diversity_specie(alpha_diversities_specie, alpha_diversities_individual)

#shf.print_alpha_diversities(alpha_diversities_individual)
#shf.print_alpha_diversities(alpha_diversities_specie)

#shf.show_boxplot(alpha_diversities_individual)
