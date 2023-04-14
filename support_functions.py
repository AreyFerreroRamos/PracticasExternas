import math


def cp_dictionary_to_array(relative_abundances_dictionary):
    relative_abundances_array = []
    
    for relative_abundance in relative_abundances_dictionary:
        relative_abundances_array.append(relative_abundances_dictionary[relative_abundance])
    return relative_abundances_array


def normalize_relative_abundances(relative_abundances_array, num_individuals):
    i = 0
    while i < len(relative_abundances_array):
        relative_abundances_array[i] /= num_individuals
        i += 1


def calculate_alpha_diversity_specie_sample_type(alpha_diversities_specie, alpha_diversities_individual, sample_type):
    for specie in alpha_diversities_individual:
        alpha_diversity_specie = 0
        for alpha_diversity_individual in alpha_diversities_individual[specie][sample_type]:
            if alpha_diversity_individual != 0:
                alpha_diversity_specie += (alpha_diversity_individual / sum(alpha_diversities_individual[specie][sample_type])) * math.log(alpha_diversity_individual / sum(alpha_diversities_individual[specie][sample_type]))
        alpha_diversities_specie[specie][sample_type] = round(0 - alpha_diversity_specie, 4)


def calculate_alpha_diversity_specie(alpha_diversities_specie, alpha_diversities_individual):
    calculate_alpha_diversity_specie_sample_type(alpha_diversities_specie, alpha_diversities_individual, 'Wild')
    calculate_alpha_diversity_specie_sample_type(alpha_diversities_specie, alpha_diversities_individual, 'Captivity')
    