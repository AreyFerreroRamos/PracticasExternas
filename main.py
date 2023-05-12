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


def relative_abundances_bacterial_genus(df_vertebrates):
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
    calculation.normalize_relative_abundances(relative_abundances, num_individuals)

    return num_zeros / num_abundances, relative_abundances


def alpha_diversities(df_vertebrates, df_metadata):
    alpha_diversities_individual = {}

    for individual in df_vertebrates:
        specie, sample_type = support.get_metadata(individual, df_metadata)

        if specie not in alpha_diversities_individual:
            alpha_diversities_individual[specie] = {'Wild': [], 'Captivity': []}

        num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

        alpha_diversity = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            if relative_abundance != 0:
                alpha_diversity += relative_abundance * math.log(relative_abundance)

        alpha_diversities_individual[specie][sample_type].append(0 - alpha_diversity)

    return alpha_diversities_individual


def matrix_abundances_individuals(df_vertebrates):
    matrix_individuals = np.empty((df_vertebrates.shape[1], df_vertebrates.shape[0]))

    num_individuals = 0
    for individual in df_vertebrates:
        num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

        column_genus = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            matrix_individuals[num_individuals][column_genus] = relative_abundance
            column_genus += 1

        num_individuals += 1

    return matrix_individuals


def matrix_abundances_vertebrates(df_vertebrates):
    matrix_vertebrates = np.empty((0, 0))

    num_genus = num_species = num_wild = num_captivity = 0
    specie_ant = ""
    for individual in df_vertebrates:
        specie, sample_type = support.get_metadata(individual, df_metadata)

        num_bacterial_species_per_individual = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            num_bacterial_species_per_individual += num_bacterial_species_per_genus
            if df_vertebrates.columns.get_loc(individual) == 0:
                num_genus += 1

        if not specie_ant or specie_ant != specie:
            if matrix_vertebrates.size == 0:
                matrix_vertebrates.resize(2, num_genus)
            else:
                calculation.normalize_matrix_vertebrates(matrix_vertebrates, num_species, num_wild, num_captivity)
                num_species += 2
                num_wild = num_captivity = 0
                matrix_vertebrates.resize((matrix_vertebrates.shape[0] + 2, matrix_vertebrates.shape[1]))
            specie_ant = specie

        column_genus = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            matrix_vertebrates[num_species + support.offset(sample_type)][column_genus] += relative_abundance
            column_genus += 1

        if sample_type == "Wild":
            num_wild += 1
        else:
            num_captivity += 1
    calculation.normalize_matrix_vertebrates(matrix_vertebrates, num_species, num_wild, num_captivity)

    return matrix_vertebrates


df_vertebrates = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)

ploter = show.select_ploter('seaborn')

if sys.argv[4] == "histogram":
    zeros_ratio, relative_abundances = relative_abundances_bacterial_genus(df_vertebrates)
    print("Total zeros: " + str(round(zeros_ratio * 100, 2)) + "%.")
    ploter.histogram(support.to_array(relative_abundances))

elif sys.argv[4] == "alpha-diversities" or sys.argv[4] == "boxplot":
    alpha_diversities = alpha_diversities(df_vertebrates, df_metadata)
    if sys.argv[4] == "alpha-diversities":
        ploter.alpha_diversities(alpha_diversities)
    else:
        calculation.t_test(alpha_diversities)
        ploter.boxplot(alpha_diversities, sys.argv[3])

else:
    if sys.argv[4] == "individuals":
        matrix_abundances = matrix_abundances_individuals(df_vertebrates)
    elif sys.argv[4] == "vertebrates":
        matrix_abundances = matrix_abundances_vertebrates(df_vertebrates)

    if sys.argv[5] == "dendrogram":
        ploter.dendrogram(matrix_abundances, sys.argv[4])

    elif sys.argv[5] == "heatmap":
        ploter.heatmap(matrix_abundances)

    elif sys.argv[5].split('-')[0] == "clustermap":
        if sys.argv[5].split('-')[1] == "log":
            calculation.log_matrix(matrix_abundances)
            ploter.cluster_map(matrix_abundances, '')
        elif sys.argv[5].split('-')[1] == "discrete":
            calculation.discretize_matrix(matrix_abundances, 0.0001)
            ploter.cluster_map(matrix_abundances, 'viridis')
        elif sys.argv[5].split('-')[1] == "fold":
            matrix_log_fold = calculation.generate_log_fold_matrix(matrix_abundances)
            ploter.cluster_map(matrix_log_fold, 'RdBu')
