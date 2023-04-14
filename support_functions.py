import math


def to_array(dictionary):
    array = []
    
    for key in dictionary:
        array.append(dictionary[key])
    return array


def normalize_relative_abundances(relative_abundances, num_individuals):
    for bacterial_genus in relative_abundances:    
        relative_abundances[bacterial_genus] /= num_individuals


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
    