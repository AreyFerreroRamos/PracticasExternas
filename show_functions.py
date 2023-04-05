from matplotlib import pyplot as plt
from matplotlib import gridspec


def print_specie(specie, f_codes_vertebrates):
    for vertebrate_specie in f_codes_vertebrates:
        if specie == vertebrate_specie.split()[0]:
            print(specie+'\t'+vertebrate_specie.split()[1].replace('_', ' ', 1))
            break


def print_alpha_diversities_sample_type(alpha_diversities, sample_type):
    print('\n'+sample_type)
    for specie in alpha_diversities:
        print(specie, alpha_diversities[specie][sample_type])


def print_alpha_diversities(alpha_diversities):
    print_alpha_diversities_sample_type(alpha_diversities, 'Wild')
    print_alpha_diversities_sample_type(alpha_diversities, 'Captivity')


def show_plot(plot_type, alpha_diversities, sample_type):
    figure = plt.figure()
    spec = gridspec.GridSpec(nrows=5, ncols=5, figure=figure)

    row = column = 0
    for specie in alpha_diversities:
        ax_box = figure.add_subplot(spec[row, column])
        if plot_type == 'Boxplot':
            ax_box.boxplot([alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity']], labels=['Wild', 'Captivity'])
        else:
            ax_box.hist(alpha_diversities[specie][sample_type], bins=[0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])

        plt.title(specie)
        
        if row == int(spec.nrows / 2) and column == 0:
            plt.ylabel("Alpha diversity")
        if row == (spec.nrows - 1) and column == int(spec.ncols / 2):
            plt.xlabel("Sample type")
        
        column += 1
        if column >= 5:
            column = 0
            row += 1

    plt.suptitle("Bacterial diversity in "+sample_type+" vertebrate species")
    plt.show()
