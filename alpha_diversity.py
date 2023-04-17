import calculation_functions as calculation
import support_functions as support
import show_functions as show
import pandas as pd
import math
import sys
import os

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

relative_abundances = {}
for bacterial_genus in df_vertebrates.index:
    relative_abundances[bacterial_genus] = 0

num_zeros = num_genus = num_individuals = 0

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
    
    num_bacterial_species_per_individual = alpha_diversity = pos = 0
    
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus
    
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        if num_bacterial_species_per_genus != 0:
            alpha_diversity += (num_bacterial_species_per_genus / num_bacterial_species_per_individual) * math.log(num_bacterial_species_per_genus / num_bacterial_species_per_individual)
            relative_abundances[df_vertebrates.index[pos]] += num_bacterial_species_per_genus / num_bacterial_species_per_individual
        else:
            num_zeros += 1
        num_genus += 1
        pos += 1

    num_individuals += 1
    alpha_diversities_individual[specie][sample_type].append(round(0 - alpha_diversity, 4))

calculation.normalize_relative_abundances(relative_abundances, num_individuals)
calculation.t_test(alpha_diversities_individual)

show.alpha_diversities(alpha_diversities_individual)
#show.t_test(calculation.t_test(alpha_diversities_individual))

show.pyplot_boxplot(alpha_diversities_individual, f_codes_vertebrates)

print("Total zeros: "+str(round(num_zeros / num_genus * 100, 2))+"%.")
show.histogram(support.to_array(relative_abundances))

f_codes_vertebrates.close()
