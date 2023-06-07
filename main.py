import data_structures_functions as data_structures
import calculation_functions as calculation
import support_functions as support
import show_functions as show
import pandas as pd
import numpy as np
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


df_vertebrates = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)

if sys.argv[4] == "alpha-diversities" or sys.argv[4] == "boxplot-grid" or sys.argv[4] == "scatterplot":
    alpha_diversities_list = data_structures.alpha_diversities_dictionary(df_vertebrates, df_metadata)

if sys.argv[4] == "distances" or sys.argv[4] == "scatterplot":
    vertebrates_relatives_abundances = data_structures.vertebrates_abundances_dictionary(df_vertebrates, df_metadata)
    vertebrates_distances, average_distances = data_structures.vertebrates_distances_list(vertebrates_relatives_abundances)

elif sys.argv[4] == "individuals":
    abundances_matrix = data_structures.abundances_individuals_matrix(df_vertebrates)

elif sys.argv[4] == "vertebrates":
    abundances_matrix = data_structures.abundances_vertebrates_matrix(df_vertebrates, df_metadata)

elif len(sys.argv[4].split()) >= 2:
    if len(sys.argv[4].split()) == 2:
        abundances_matrix = data_structures.abundances_matrices_specie(support.get_code_specie(
            sys.argv[4], sys.argv[3]), df_vertebrates, df_metadata)
    elif sys.argv[4].split()[2] == "Wild":
        abundances_matrix = data_structures.abundances_matrices_specie_sample_type(support.get_code_specie(
            sys.argv[4].split()[0] + ' ' + sys.argv[4].split()[1], sys.argv[3]), 'Wild', df_vertebrates, df_metadata)
    elif sys.argv[4].split()[2] == "Captive":
        abundances_matrix = data_structures.abundances_matrices_specie_sample_type(support.get_code_specie(
            sys.argv[4].split()[0] + ' ' + sys.argv[4].split()[1], sys.argv[3]), 'Captivity', df_vertebrates, df_metadata)


ploter = show.select_ploter('seaborn')

if sys.argv[4] == "alpha-diversities":
    ploter.alpha_diversities(alpha_diversities_list)

elif sys.argv[4] == "histogram":
    zeros_ratio, global_relative_abundances = data_structures.relative_abundances_bacterial_genus_dictionary(df_vertebrates)
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
    show.select_ploter('pyplot').scatterplot(alpha_average, distance_average, sys.argv[3])

elif sys.argv[5] == "dendrogram":
    ploter.dendrogram(abundances_matrix, sys.argv[4])

elif sys.argv[5] == "heatmap":
    ploter.heatmap(abundances_matrix)

elif sys.argv[5].split(' ')[0] == "clustermap":
    if sys.argv[5].split()[1] == "log":
        calculation.log_matrix(abundances_matrix)
        ploter.cluster_map(abundances_matrix, '')

    elif sys.argv[5].split()[1] == "discrete":
        calculation.discretize_matrix(abundances_matrix, 0.0001)
        ploter.cluster_map(abundances_matrix, 'viridis')
        nested_abundances_matrix, p_value = calculation.nestedness_assessment(abundances_matrix, 1000)
        print(nested_abundances_matrix, p_value)

    elif sys.argv[5].split()[1] == "fold":
        matrix_log_fold = calculation.generate_log_fold_matrix(abundances_matrix)
        ploter.cluster_map(matrix_log_fold, 'RdBu')

# matrix = np.array([[1, 1, 1, 0], [1, 1, 1, 1], [1, 0, 0, 0], [1, 1, 0, 0]])
# nested_matrix, p_value = calculation.nestedness_assessment(matrix, 1000)
# print(nested_matrix, p_value)
