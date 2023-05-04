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
