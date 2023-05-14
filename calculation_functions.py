from scipy import stats
import fastcluster
import numpy as np
import math


def normalize_relative_abundances(relative_abundances, num_individuals):
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


def normalize_vertebrate_abundances(vertebrates_relatives_abundances, specie, num_wild, num_captivity):
    num_abundance = 0
    while num_abundance < len(vertebrates_relatives_abundances[specie]['Wild']):
        vertebrates_relatives_abundances[specie]['Wild'][num_abundance] /= num_wild
        num_abundance += 1

    num_abundance = 0
    while num_abundance < len(vertebrates_relatives_abundances[specie]['Captivity']):
        vertebrates_relatives_abundances[specie]['Captivity'][num_abundance] /= num_captivity
        num_abundance += 1


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
