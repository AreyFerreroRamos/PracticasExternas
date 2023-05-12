def get_num_species_per_individual(individual, df_vertebrates):
    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus

    return num_bacterial_species_per_individual


def get_specie_sample_type(individual, df_metadata):
    row = 1
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            specie = df_metadata.loc[row, df_metadata.columns[2]]
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
        else:
            row += 1

    return specie, sample_type


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


def pad_array(first_array, second_array):
    if len(first_array) < len(second_array):
        i = len(first_array)
        while i < len(second_array):
            first_array.append(sum(first_array) / len(first_array))
            i += 1
    elif len(first_array) > len(second_array):
        i = len(second_array)
        while i < len(first_array):
            second_array.append(sum(second_array) / len(second_array))
            i += 1
