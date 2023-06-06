import calculation_functions as calculation
import support_functions as support
import numpy as np
import math


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


def abundances_vertebrates_matrix(df_vertebrates, df_metadata):
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


def vertebrates_abundances_dictionary(df_vertebrates, df_metadata):
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


def abundances_matrices_specie(vertebrate_specie, df_vertebrates, df_metadata):
    specie_matrix = np.empty((0, df_vertebrates.shape[0]))

    num_individuals = 0
    for individual in df_vertebrates:
        specie, _ = support.get_specie_sample_type(individual, df_metadata)

        if specie == vertebrate_specie:
            specie_matrix.resize((specie_matrix.shape[0] + 1, specie_matrix.shape[1]))

            num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

            column_genus = 0
            for num_bacterial_species_per_genus in df_vertebrates[individual]:
                relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
                specie_matrix[num_individuals][column_genus] = relative_abundance
                column_genus += 1

            num_individuals += 1

    return specie_matrix


def abundances_matrices_specie_sample_type(vertebrate_specie, vertebrate_sample_type, df_vertebrates, df_metadata):
    specie_matrix = np.empty((0, df_vertebrates.shape[0]))

    num_individuals = 0
    for individual in df_vertebrates:
        specie, sample_type = support.get_specie_sample_type(individual, df_metadata)

        if specie == vertebrate_specie and sample_type == vertebrate_sample_type:
            specie_matrix.resize((specie_matrix.shape[0] + 1, specie_matrix.shape[1]))

            num_bacterial_species_per_individual = support.get_num_species_per_individual(individual, df_vertebrates)

            column_genus = 0
            for num_bacterial_species_per_genus in df_vertebrates[individual]:
                relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
                specie_matrix[num_individuals][column_genus] = relative_abundance
                column_genus += 1

            num_individuals += 1

    return specie_matrix
