import support_functions as support
from scipy import stats
import fastcluster
import numpy as np
import math


def normalize_relative_abundances_dictionary(relative_abundances, num_individuals):
    for bacterial_genus in relative_abundances:
        relative_abundances[bacterial_genus] /= num_individuals


def t_test(alpha_diversities):
    for specie in alpha_diversities:
        t_stat, p_value = stats.ttest_ind(alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity'],
                                          equal_var=False)
        alpha_diversities[specie]['t_stat'] = t_stat
        alpha_diversities[specie]['p_value'] = p_value


def hierarchical_clustering(matrix):
    return fastcluster.linkage(matrix, method='average', metric='euclidean')


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


def distances(relative_abundances):
    distances = []

    for first_individual in relative_abundances:
        for second_individual in relative_abundances:
            if first_individual < second_individual:
                first_list_abundances, second_list_abundances = support.pad_list_zeros(
                    relative_abundances[first_individual], relative_abundances[second_individual])
                distances.append(np.linalg.norm(np.array(first_list_abundances) - np.array(second_list_abundances)))

    return distances, sum(distances) / len(distances)


def nestedness(matrix):
    nestedness = 0

    row_sums = np.sum(matrix, axis=1)
    col_sums = np.sum(matrix, axis=0)
    print(row_sums)
    print(col_sums)

    mean_row_sums = np.mean(row_sums)
    mean_col_sums = np.mean(col_sums)
    print(mean_row_sums)
    print(mean_col_sums)

    for rows in range(matrix.shape[0]):
        for columns in range(matrix.shape[1]):
            if matrix[rows][columns] == 1:
                print(matrix[rows][columns])

    print(matrix)
