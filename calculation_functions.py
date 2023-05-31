from scipy import stats
import fastcluster
import numpy as np
import random
import math


def t_test(alpha_diversities):
    for specie in alpha_diversities:
        t_stat, p_value = stats.ttest_ind(alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity'],
                                          equal_var=False)
        alpha_diversities[specie]['t_stat'] = t_stat
        alpha_diversities[specie]['p_value'] = p_value


def hierarchical_clustering(matrix):
    return fastcluster.linkage(matrix, method='average', metric='euclidean')


def average_alpha_diversities(alpha_diversities):
    alpha_average = {}

    for specie in alpha_diversities:
        alpha_average[specie] = {}
        for sample_type in alpha_diversities[specie]:
            alpha_average[specie][sample_type] = sum(alpha_diversities[specie][sample_type]) /\
                                                 len(alpha_diversities[specie][sample_type])
    return alpha_average


def normalize_relative_abundances_dictionary(relative_abundances, num_individuals):
    for bacterial_genus in relative_abundances:
        relative_abundances[bacterial_genus] /= num_individuals


def normalize_matrix(matrix, rows, columns, num_individuals):
    if num_individuals != 0:
        matrix[rows][columns] = matrix[rows][columns] / num_individuals
    else:
        matrix[rows][columns] = 0


def normalize_matrix_vertebrates(matrix_vertebrate_genus, num_specie, num_wild, num_captivity):
    num_genus = 0
    while num_genus < matrix_vertebrate_genus.shape[1]:
        normalize_matrix(matrix_vertebrate_genus, num_specie, num_genus, num_wild)
        normalize_matrix(matrix_vertebrate_genus, num_specie + 1, num_genus, num_captivity)
        num_genus += 1


def normalize_relative_abundances_list(relative_abundances, num_individuals):
    num_abundance = 0
    while num_abundance < len(relative_abundances):
        relative_abundances[num_abundance] /= num_individuals
        num_abundance += 1


def normalize_vertebrate_abundances(vertebrates_relatives_abundances, specie, num_wild, num_captivity):
    normalize_relative_abundances_list(vertebrates_relatives_abundances[specie]['Wild'], num_wild)
    normalize_relative_abundances_list(vertebrates_relatives_abundances[specie]['Captivity'], num_captivity)


def log_matrix(matrix):
    rows = 0
    while rows < matrix.shape[0]:
        columns = 0
        while columns < matrix.shape[1]:
            if matrix[rows][columns] > 0:
                matrix[rows][columns] = math.log(matrix[rows][columns], 10)
            else:
                matrix[rows][columns] = -10
            columns += 1
        rows += 1


def discretize_matrix(matrix, threshold):
    rows = 0
    while rows < matrix.shape[0]:
        columns = 0
        while columns < matrix.shape[1]:
            if matrix[rows][columns] > threshold:
                matrix[rows][columns] = 1
            else:
                matrix[rows][columns] = 0
            columns += 1
        rows += 1


def generate_log_fold_matrix(matrix_vertebrates_genus_sample_type):
    matrix_vertebrates_genus = np.empty((int(matrix_vertebrates_genus_sample_type.shape[0] / 2),
                                         matrix_vertebrates_genus_sample_type.shape[1]))
    num_genus = 0
    while num_genus < matrix_vertebrates_genus_sample_type.shape[1]:
        num_specie = 0
        while num_specie < matrix_vertebrates_genus_sample_type.shape[0]:
            if matrix_vertebrates_genus_sample_type[num_specie][num_genus] > 0 and \
                    matrix_vertebrates_genus_sample_type[num_specie + 1][num_genus] > 0:
                matrix_vertebrates_genus[int(num_specie / 2)][num_genus] = math.log(
                    matrix_vertebrates_genus_sample_type[num_specie][num_genus] /
                    matrix_vertebrates_genus_sample_type[num_specie + 1][num_genus], 10)
            else:
                matrix_vertebrates_genus[int(num_specie / 2)][num_genus] = 0
            num_specie += 2
        num_genus += 1

    return matrix_vertebrates_genus


def intra_distances(relative_abundances):
    distances = []

    for first_individual in relative_abundances:
        for second_individual in relative_abundances:
            if first_individual < second_individual:
                distances.append(np.linalg.norm(np.array(relative_abundances[first_individual]) -
                                                np.array(relative_abundances[second_individual])))

    return distances, sum(distances) / len(distances)


def inter_distances(first_relative_abundances, second_relative_abundances):
    distances = []

    for first_individual in first_relative_abundances:
        for second_individual in second_relative_abundances:
            distances.append(np.linalg.norm(np.array(first_relative_abundances[first_individual]) -
                                            np.array(second_relative_abundances[second_individual])))

    return distances, sum(distances) / len(distances)


def random_pairs(vertebrate_relative_abundances, total_couples):
    distance = 0
    num_couples = 0

    while num_couples < total_couples:
        first_specie, second_specie = random.sample(vertebrate_relative_abundances.keys(), 2)
        first_sample_type, _ = random.sample(vertebrate_relative_abundances[first_specie].keys(), 2)
        second_sample_type, _ = random.sample(vertebrate_relative_abundances[second_specie].keys(), 2)
        first_individual, _ = random.sample(vertebrate_relative_abundances[first_specie][first_sample_type].keys(), 2)
        second_individual, _ = random.sample(vertebrate_relative_abundances[second_specie][second_sample_type].keys(), 2)
        distance += np.linalg.norm(
            np.array(vertebrate_relative_abundances[first_specie][first_sample_type][first_individual]) -
            np.array(vertebrate_relative_abundances[second_specie][second_sample_type][second_individual]))
        num_couples += 1

    return distance / num_couples


def nestedness_assessment(matrix):
    if matrix.size == 0:
        return 0

    sorted_matrix = sort_matrix(matrix)
    rows, columns = sorted_matrix.shape
    nestedness_value = nestedness(sorted_matrix, rows, columns)

    return nestedness_value


def sort_matrix(matrix):
    # Calculate the summary of rows and columns in the matrix.
    acum_rows = np.sum(matrix, axis=1)
    acum_cols = np.sum(matrix, axis=0)

    # Calculate the rows and columns order based in the summatory of them in the matrix.
    sorted_rows = np.argsort(acum_rows)[::-1]
    sorted_cols = np.argsort(acum_cols)[::-1]

    # Reorder and return the matrix.
    return matrix[np.ix_(sorted_rows, sorted_cols)]


def nestedness(matrix, rows, columns):
    first_isocline = 0
    second_isocline = 0
    third_isocline = 0
    fourth_isocline = 0

    # Calculate the sum of the number of shared interactions between rows.
    for first_row in range(rows):
        for second_row in range(rows):
            if first_row < second_row:
                for col in range(columns):
                    if matrix[first_row][col] == 1 and matrix[second_row][col] == 1:
                        first_isocline += 1

    # Calculate the sum of the number of shared interactions between columns.
    for first_col in range(columns):
        for second_col in range(columns):
            if first_col < second_col:
                for row in range(rows):
                    if matrix[row][first_col] == 1 and matrix[row][second_col] == 1:
                        second_isocline += 1

    # Calculate the sum of the number of interactions of rows.
    for first_row in range(rows):
        for second_row in range(rows):
            if first_row < second_row:
                first_acum = 0
                second_acum = 0
                for col in range(columns):
                    first_acum += matrix[first_row][col]
                    second_acum += matrix[second_row][col]
                third_isocline += min(first_acum, second_acum)


    # Calculate the sum of the number of interactions of columns.
    for first_col in range(columns):
        for second_col in range(columns):
            if first_col < second_col:
                first_acum = 0
                second_acum = 0
                for row in range(rows):
                    first_acum += matrix[row][first_col]
                    second_acum += matrix[row][second_col]
                fourth_isocline += min(first_acum, second_acum)

    # Calculate and return the nestedness of the matrix.
    return (first_isocline + second_isocline) / (third_isocline + fourth_isocline)
