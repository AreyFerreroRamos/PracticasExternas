from scipy import stats
import fastcluster
import math


def t_test(alpha_diversities):
    for specie in alpha_diversities:
        t_stat, p_value = stats.ttest_ind(alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity'],
                                          equal_var=False)
        alpha_diversities[specie]['t_stat'] = t_stat
        alpha_diversities[specie]['p_value'] = p_value


def hierarchical_clustering(matrix):
    return fastcluster.linkage(matrix)


def normalize_matrix_vertebrates_genus(matrix_vertebrate_genus, num_specie, num_individuals_wild, num_individuals_captivity):
    num_genus = 0
    while num_genus < matrix_vertebrate_genus.shape[1]:
        if num_individuals_wild != 0:
            matrix_vertebrate_genus[num_specie][num_genus] /= num_individuals_wild
        else:
            matrix_vertebrate_genus[num_specie][num_genus] = 0
        if num_individuals_captivity != 0:
            matrix_vertebrate_genus[num_specie + 1][num_genus] /= num_individuals_captivity
        else:
            matrix_vertebrate_genus[num_specie + 1][num_genus] = 0
        num_genus += 1


def normalize_relative_abundances(relative_abundances, num_individuals):
    for bacterial_genus in relative_abundances:    
        relative_abundances[bacterial_genus] /= num_individuals


def alpha_diversity_specie_sample_type(alpha_diversities_specie, alpha_diversities_individual, sample_type):
    for specie in alpha_diversities_individual:
        alpha_diversity_specie = 0
        for alpha_diversity_individual in alpha_diversities_individual[specie][sample_type]:
            if alpha_diversity_individual != 0:
                alpha_diversity_specie += (alpha_diversity_individual / sum(alpha_diversities_individual[specie][sample_type])) * math.log(alpha_diversity_individual / sum(alpha_diversities_individual[specie][sample_type]))
        alpha_diversities_specie[specie][sample_type] = round(0 - alpha_diversity_specie, 4)


def alpha_diversity_specie(alpha_diversities_specie, alpha_diversities_individual):
    alpha_diversity_specie_sample_type(alpha_diversities_specie, alpha_diversities_individual, 'Wild')
    alpha_diversity_specie_sample_type(alpha_diversities_specie, alpha_diversities_individual, 'Captivity')
