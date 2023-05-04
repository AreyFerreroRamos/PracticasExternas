import calculation_functions as calculation
import support_functions as support
import show_functions as show
import pandas as pd
import numpy as np
import math
import sys
import os


if len(sys.argv) != 7:
    print("Error: The number of parameters is incorrect. Three files are needed.")
    exit()

if not os.path.isfile(sys.argv[1]):
    print("Error: The first parameter must be a file corresponding to all vertebrate species.")
    exit()

if not os.path.isfile(sys.argv[2]):
    print("Error: The second parameter must be a file corresponding to the metadata of all the samples in the study.")
    exit()

if not os.path.isfile(sys.argv[3]):
    print("Error: The third parameter must be a file with a list of the scientific name that corresponds to the code "
          "that identifies a vertebrate species.")
    exit()

df_vertebrates = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)

alpha_diversities_individual = {}

matrix_individuals_genus = np.empty((0, 0))
matrix_vertebrates_genus = np.empty((0, 0))

relative_abundances = {}
for bacterial_genus in df_vertebrates.index:
    relative_abundances[bacterial_genus] = 0

num_zeros = num_abundances = num_individuals = num_genus = num_species = num_wild = num_captivity = 0

for individual in df_vertebrates:
    row = 1
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            specie = df_metadata.loc[row, df_metadata.columns[2]]
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
        else:
            row += 1

    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus
        if df_vertebrates.columns.get_loc(individual) == 0:
            num_genus += 1

    if matrix_individuals_genus.size == 0:
        matrix_individuals_genus.resize((1, num_genus))
    else:
        matrix_individuals_genus.resize((matrix_individuals_genus.shape[0] + 1, matrix_individuals_genus.shape[1]))

    if specie not in alpha_diversities_individual:
        alpha_diversities_individual[specie] = {'Wild': [], 'Captivity': []}
        if matrix_vertebrates_genus.size == 0:
            matrix_vertebrates_genus.resize((2, num_genus))
        else:
            calculation.normalize_matrix_vertebrates(matrix_vertebrates_genus, num_species, num_wild, num_captivity)
            num_species += 2
            num_wild = num_captivity = 0
            matrix_vertebrates_genus.resize((matrix_vertebrates_genus.shape[0] + 2, matrix_vertebrates_genus.shape[1]))

    alpha_diversity = pos = column_genus = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
        if relative_abundance != 0:
            alpha_diversity += relative_abundance * math.log(relative_abundance)
            relative_abundances[df_vertebrates.index[pos]] += relative_abundance
        else:
            num_zeros += 1
        num_abundances += 1
        pos += 1

        matrix_individuals_genus[num_individuals][column_genus] = relative_abundance
        matrix_vertebrates_genus[num_species + support.offset(sample_type)][column_genus] += relative_abundance
        column_genus += 1

    if sample_type == 'Wild':
        num_wild += 1
    else:
        num_captivity += 1
    num_individuals += 1
    alpha_diversities_individual[specie][sample_type].append(round(0 - alpha_diversity, 4))

if sys.argv[4] == "alpha-diversities":
    show.alpha_diversities(alpha_diversities_individual)
else:
    ploter = show.select_ploter(sys.argv[6])
    if sys.argv[4] == "histogram":
        calculation.normalize_relative_abundances(relative_abundances, num_individuals)
        print("Total zeros: "+str(round(num_zeros / num_abundances * 100, 2))+"%.")
        ploter.histogram(support.to_array(relative_abundances))
    elif sys.argv[4] == "boxplot":
        calculation.t_test(alpha_diversities_individual)
        ploter.boxplot(alpha_diversities_individual, sys.argv[3])
    else:
        calculation.normalize_matrix_vertebrates(matrix_vertebrates_genus, num_species, num_wild, num_captivity)

        if sys.argv[4] == "dendrogram":
            ploter.dendrogram(matrix_individuals_genus, 'Individuals')
            ploter.dendrogram(matrix_vertebrates_genus, 'Vertebrate species')
        elif sys.argv[4] == "heatmap":
            ploter.heatmap(matrix_individuals_genus)
            ploter.heatmap(matrix_vertebrates_genus)
        elif sys.argv[4] == "clustermap":
            if sys.argv[5] == "individuals":
                calculation.log_matrix(matrix_individuals_genus)
                ploter.cluster_map(matrix_individuals_genus, '')
            elif sys.argv[5] == "vertebrates":
                calculation.log_matrix(matrix_vertebrates_genus)
                ploter.cluster_map(matrix_vertebrates_genus, '')
            elif sys.argv[5] == "discrete-individuals":
                calculation.discretize_matrix(matrix_individuals_genus, 0.0000001)
                ploter.cluster_map(matrix_individuals_genus, 'viridis')
            elif sys.argv[5] == "discrete-vertebrates":
                calculation.discretize_matrix(matrix_vertebrates_genus, 0.0000001)
                ploter.cluster_map(matrix_vertebrates_genus, 'viridis')
            elif sys.argv[5] == "log-fold-vertebrates":
                matrix_vertebrates_genus_log_fold = calculation.generate_log_fold_matrix(matrix_vertebrates_genus)
                ploter.cluster_map(matrix_vertebrates_genus_log_fold, 'RdBu')

