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


def relative_abundances_bacterial_genus_dictionary(df_vertebrates):
    relative_abundances = {}
    for bacterial_genus in df_vertebrates.index:
        relative_abundances[bacterial_genus] = 0

    num_zeros = num_abundances = num_individuals = 0
    for individual in df_vertebrates:
        num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

        pos = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            if relative_abundance != 0:
                relative_abundances[df_vertebrates.index[pos]] += relative_abundance
            else:
                num_zeros += 1
            num_abundances += 1
            pos += 1

        num_individuals += 1

    calculation.normalize_relative_abundances_dictionary(relative_abundances, num_individuals)
    return num_zeros / num_abundances, relative_abundances


def alpha_diversities_dictionary(df_vertebrates, df_metadata):
    alpha_diversities = {}

    for individual in df_vertebrates:
        specie, sample_type = support.get_specie_sample_type(individual, df_metadata)

        if specie not in alpha_diversities:
            alpha_diversities[specie] = {'Wild': [], 'Captivity': []}

        num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

        alpha_diversity = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            if relative_abundance != 0:
                alpha_diversity += relative_abundance * math.log(relative_abundance)

        alpha_diversities[specie][sample_type].append(0 - alpha_diversity)

    return alpha_diversities


def abundances_individuals_matrix(df_vertebrates):
    individuals_matrix = np.empty((df_vertebrates.shape[1], df_vertebrates.shape[0]))

    num_individuals = 0
    for individual in df_vertebrates:
        num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

        column_genus = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            individuals_matrix[num_individuals][column_genus] = relative_abundance
            column_genus += 1

        num_individuals += 1

    return individuals_matrix


def abundances_vertebrates_matrix(df_vertebrates):
    vertebrates_matrix = np.zeros((2, df_vertebrates.shape[0]))
    specie_ant, _ = support.get_specie_sample_type(df_vertebrates.columns[0], df_metadata)

    num_species = num_wild = num_captivity = 0
    for individual in df_vertebrates:
        specie, sample_type = support.get_specie_sample_type(individual, df_metadata)

        if specie_ant != specie:
            calculation.normalize_matrix_vertebrates(vertebrates_matrix, num_species, num_wild, num_captivity)
            num_species += 2
            num_wild = num_captivity = 0
            vertebrates_matrix.resize((vertebrates_matrix.shape[0] + 2, vertebrates_matrix.shape[1]))
            specie_ant = specie

        num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

        column_genus = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            vertebrates_matrix[num_species + support.offset(sample_type)][column_genus] += relative_abundance
            column_genus += 1

        if sample_type == "Wild":
            num_wild += 1
        else:
            num_captivity += 1

    calculation.normalize_matrix_vertebrates(vertebrates_matrix, num_species, num_wild, num_captivity)
    return vertebrates_matrix


def vertebrates_abundances_dictionary(df_vertebrates):
    vertebrates_relative_abundances = {}

    for individual in df_vertebrates:
        specie, sample_type = support.get_specie_sample_type(individual, df_metadata)

        if specie not in vertebrates_relative_abundances:
            vertebrates_relative_abundances[specie] = {'Wild': {}, 'Captivity': {}}

        vertebrates_relative_abundances[specie][sample_type][individual] = []

        num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            vertebrates_relative_abundances[specie][sample_type][individual].append(relative_abundance)

    return vertebrates_relative_abundances


def vertebrates_distances_list(vertebrates_relative_abundances):
    vertebrate_distances = {}
    average_distances = {}

    for specie in vertebrates_relative_abundances:
        vertebrate_distances[specie] = {}
        average_distances[specie] = {}

        vertebrate_distances[specie]['Wild'], average_distances[specie]['Wild'] = calculation.intra_distances(
            vertebrates_relative_abundances[specie]['Wild'])
        vertebrate_distances[specie]['Captivity'], average_distances[specie]['Captivity'] = calculation.intra_distances(
            vertebrates_relative_abundances[specie]['Captivity'])
        vertebrate_distances[specie]['Wild-Captivity'], average_distances[specie]['Wild-Captivity'] = calculation.\
            inter_distances(vertebrates_relative_abundances[specie]['Wild'],
                            vertebrates_relative_abundances[specie]['Captivity'])

    return vertebrate_distances, average_distances


df_vertebrates = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)

if sys.argv[4] == "alpha-diversities" or sys.argv[4] == "boxplot-grid" or sys.argv[4] == "scatterplot":
    alpha_diversities_list = alpha_diversities_dictionary(df_vertebrates, df_metadata)

if sys.argv[4] == "distances" or sys.argv[4] == "scatterplot":
    vertebrates_relatives_abundances = vertebrates_abundances_dictionary(df_vertebrates)
    vertebrates_distances, average_distances = vertebrates_distances_list(vertebrates_relatives_abundances)

if sys.argv[4] == "individuals":
    abundances_matrix = abundances_individuals_matrix(df_vertebrates)
elif sys.argv[4] == "vertebrates":
    abundances_matrix = abundances_vertebrates_matrix(df_vertebrates)


ploter = show.select_ploter('seaborn')

if sys.argv[4] == "alpha-diversities":
    ploter.alpha_diversities(alpha_diversities_list)

elif sys.argv[4] == "histogram":
    zeros_ratio, global_relative_abundances = relative_abundances_bacterial_genus_dictionary(df_vertebrates)
    print("Total zeros: " + str(round(zeros_ratio * 100, 2)) + "%.")
    ploter.histogram(support.to_array(global_relative_abundances))

elif sys.argv[4] == "boxplot-grid":
    calculation.t_test(alpha_diversities_list)
    ploter.boxplot_grid(alpha_diversities_list, sys.argv[3])

elif sys.argv[4] == "distances":
    base_line = calculation.random_pairs(vertebrates_relatives_abundances, 10000)
    ploter.boxplot(vertebrates_distances, base_line, sys.argv[3])
    ploter.stemplot(average_distances, sys.argv[3])

elif sys.argv[4] == "scatterplot":
    alpha_average = calculation.average_alpha_diversities(alpha_diversities_list)
    distance_average = support.reduce_average_distances(average_distances)
    ploter.scatterplot(alpha_average, distance_average, sys.argv[3])

elif sys.argv[5] == "dendrogram":
    ploter.dendrogram(abundances_matrix, sys.argv[4])

elif sys.argv[5] == "heatmap":
    ploter.heatmap(abundances_matrix)

elif sys.argv[5].split('-')[0] == "clustermap":
    if sys.argv[5].split('-')[1] == "log":
        calculation.log_matrix(abundances_matrix)
        ploter.cluster_map(abundances_matrix, '')

    elif sys.argv[5].split('-')[1] == "discrete":
        calculation.discretize_matrix(abundances_matrix, 0.0001)
        ploter.cluster_map(abundances_matrix, 'viridis')
        nested_abundances_matrix, nested_randomized_matrix = calculation.nestedness_assessment(abundances_matrix)
        print(nested_abundances_matrix, nested_randomized_matrix)

    elif sys.argv[5].split('-')[1] == "fold":
        matrix_log_fold = calculation.generate_log_fold_matrix(abundances_matrix)
        ploter.cluster_map(matrix_log_fold, 'RdBu')

# matrix = np.array([[1, 1, 1, 0], [1, 1, 1, 1], [1, 0, 0, 0], [1, 1, 0, 0]])
# print(calculation.nestedness(matrix))
# random_matrix = np.zeros((matrix.shape[0], matrix.shape[1]), dtype=int)
# random_matrix.ravel()[np.random.choice(matrix.shape[0] * matrix.shape[1],
#                                       calculation.count_ones_binary_matrix(matrix), replace=False)] = 1
# print(calculation.nestedness(random_matrix))
