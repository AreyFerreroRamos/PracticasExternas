def get_num_species_per_individual(individual, df_vertebrates):
    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus

    return num_bacterial_species_per_individual


def get_specie_sample_type(individual, df_metadata):
    specie = sample_type = ""

    row = 1
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            specie = df_metadata.loc[row, df_metadata.columns[2]]
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
        else:
            row += 1

    return specie, sample_type


def get_name_specie(specie, name_file_codes_vertebrates):
    f_codes_vertebrates = open(name_file_codes_vertebrates, 'r')
    name_specie = ""

    for vertebrate_specie in f_codes_vertebrates:
        if specie == vertebrate_specie.split()[0]:
            name_specie = vertebrate_specie.split()[1].replace('_', ' ', 1)

    f_codes_vertebrates.close()
    return name_specie


def offset(sample_type):
    if sample_type == 'Wild':
        return 0
    else:
        return 1


def to_array(dictionary):
    array = []
    
    for key in dictionary:
        array.append(dictionary[key])
    return array


def generate_plot_data(distances, name_file_code_vertebrates):
    data = []
    labels = []

    for specie in distances:
        name_specie = get_name_specie(specie, name_file_code_vertebrates)
        name_specie = name_specie[0] + '. ' + name_specie.split(' ')[1]

        data.append(distances[specie]['Wild'])
        labels.append(name_specie + ' (W)')

        data.append(distances[specie]['Captivity'])
        labels.append(name_specie + ' (C)')

        data.append(distances[specie]['Wild-Captivity'])
        labels.append(name_specie + ' (W-C)')

    return data, labels


def generate_scatterplot_data(alpha_average, distance_average, name_file_code_vertebrates):
    average_alpha = []
    average_distance = []
    labels = []

    for specie in alpha_average:
        average_alpha.append(alpha_average[specie]['Wild'])
        average_distance.append(distance_average[specie]['Wild'])

        average_alpha.append(alpha_average[specie]['Captivity'])
        average_distance.append(distance_average[specie]['Captivity'])

        name_specie = get_name_specie(specie, name_file_code_vertebrates)
        labels.append(name_specie[0] + '. ' + name_specie.split(' ')[1])

    return average_alpha, average_distance, labels


def reduce_average_distances(average_distances):
    distances_average = {}

    for specie in average_distances:
        distances_average[specie] = {}
        distances_average[specie]['Wild'] = average_distances[specie]['Wild']
        distances_average[specie]['Captivity'] = average_distances[specie]['Captivity']

    return distances_average


def pad_list_average(first_list, second_list):
    if len(first_list) < len(second_list):
        first_list += [sum(first_list) / len(first_list)] * (len(second_list) - len(first_list))
    else:
        second_list += [sum(second_list) / len(second_list)] * (len(first_list) - len(second_list))

    return first_list, second_list


def create_colors():
    return [(1.0, 0.0, 0.0, 1.0),     # Red
            (0.0, 1.0, 0.0, 1.0),     # Green
            (0.0, 0.0, 1.0, 1.0),     # Blue
            (1.0, 1.0, 0.0, 1.0),     # Yellow
            (1.0, 0.0, 1.0, 1.0),     # Magenta
            (0.0, 1.0, 1.0, 1.0),     # Cyan
            (1.0, 0.5, 0.0, 1.0),     # Orange
            (0.5, 0.0, 1.0, 1.0),     # Purple
            (0.0, 1.0, 0.5, 1.0),     # Lime
            (0.0, 0.0, 0.0, 1.0),     # Black
            (0.0, 0.5, 1.0, 1.0),     # Azure
            (1.0, 0.0, 0.5, 1.0),     # Crimson
            (0.5, 0.0, 0.0, 1.0),     # Maroon
            (0.0, 0.5, 0.0, 1.0),     # Olive
            (0.0, 0.0, 0.5, 1.0),     # Navy
            (0.5, 0.5, 0.5, 1.0),     # Gray
            (1.0, 0.5, 0.5, 1.0),     # Pink
            (0.75, 1.0, 0.75, 1.0),   # Light green
            (0.5, 0.5, 1.0, 1.0),     # Light blue
            (1.0, 0.5, 1.0, 1.0),     # Fuchsia
            (1.0, 1.0, 0.75, 1.0),    # Light yellow
            (0.75, 1.0, 1.0, 1.0),    # Light cyan
            (0.5, 1.0, 0.0, 1.0),     # Chartreuse
            (0.75, 0.25, 0.75, 1.0),
            (0.25, 0.75, 0.25, 1.0)]
