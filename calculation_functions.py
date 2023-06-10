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


def pearson_correlation(average_alpha, average_distance):
    pearson_corr, p_value = stats.pearsonr(average_alpha, average_distance)
    return pearson_corr, p_value


def spearman_correlation(average_alpha, average_distance):
    spearman_corr, p_value = stats.pearsonr(average_alpha, average_distance)
    return spearman_corr, p_value


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


def sort_matrix(matrix):
    # Calculate the summary of rows and columns in the matrix.
    acum_rows = np.sum(matrix, axis=1)
    acum_cols = np.sum(matrix, axis=0)

    # Calculate the rows and columns order based in the summatory of them in the matrix.
    sorted_rows = np.argsort(acum_rows)[::-1]
    sorted_cols = np.argsort(acum_cols)[::-1]

    # Reorder and return the matrix.
    return matrix[np.ix_(sorted_rows, sorted_cols)]


def nestedness(matrix):
    first_isocline = second_isocline = third_isocline = fourth_isocline = 0
    sum_rows = []
    sum_cols = []

    # Calculate and save the num of interactions of every row.
    for row in range(matrix.shape[0]):
        sum_rows.append(sum(matrix[row, :]))

    # Calculate and save the num of interactions of every col.
    for col in range(matrix.shape[1]):
        sum_cols.append(sum(matrix[:, col]))

    # Calculate the sum of the number of shared interactions between rows
    # and the sum of the minimum of pairs of interactions of rows.
    for first_row in range(matrix.shape[0] - 1):
        for second_row in range(first_row + 1, matrix.shape[0]):
            first_isocline += sum([first * second for first, second in zip(matrix[first_row, :], matrix[second_row, :])])
            third_isocline += min(sum_rows[first_row], sum_rows[second_row])

    # Calculate the sum of the number of shared interactions between columns
    # and the sum of the minimum of pairs of the number of interactions of columns.
    for first_col in range(matrix.shape[1] - 1):
        for second_col in range(first_col + 1, matrix.shape[1]):
            second_isocline += sum([first * second for first, second in zip(matrix[:, first_col], matrix[:, second_col])])
            fourth_isocline += min(sum_cols[first_col], sum_cols[second_col])

    # Calculate and return the nestedness of the matrix.
    return (first_isocline + second_isocline) / (third_isocline + fourth_isocline)


def count_ones_binary_matrix(matrix):
    num_ones = 0

    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            if matrix[row][column] == 1:
                num_ones += 1

    return num_ones


def generate_nested_values_randomized(matrix, num_randomized_matrices):
    nested_values_randomized = []
    num_ones = count_ones_binary_matrix(matrix)

    for i in range(num_randomized_matrices):
        randomized_matrix = np.zeros((matrix.shape[0], matrix.shape[1]), dtype=int)
        randomized_matrix.ravel()[np.random.choice(matrix.shape[0] * matrix.shape[1], num_ones, replace=False)] = 1
        nested_values_randomized.append(nestedness(randomized_matrix))

    return nested_values_randomized


def nestedness_assessment(matrix, num_randomized_matrices):
    if matrix.size == 0:
        return 0, 0

    # Generate as many randomized matrices from the real matrix as it is specified by parameter
    # and calculate their nestedness value.
    # nested_values = generate_nested_values_randomized(matrix, num_randomized_matrices)

    # Calculate the nestedness value of the real matrix.
    nested_value = nestedness(matrix)
    # nested_values.append(nested_value)

    # Sort the list of nestedness values.
    # nested_values.sort()

    # Calculate the fraction of randomized matrices that have a nestedness value greater than that of the real matrix.
    # p_value = (num_randomized_matrices - nested_values.index(nested_value)) / (num_randomized_matrices + 1)

    return nested_value
