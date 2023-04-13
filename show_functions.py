from scipy import stats
import seaborn as sns
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


def t_test(alpha_diversities):
    print("\nStudent's t-test")
    for specie in alpha_diversities:
        # Show gaussians.
        t_statistic, p_value = stats.ttest_ind(alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity'], equal_var=False)
        print(specie+":\tt statistic = "+str(round(t_statistic, 10))+"\tp-value = "+str(round(p_value, 10)))


def show_boxplot(alpha_diversities):
    #figure = plt.figure()
    #spec = gridspec.GridSpec(nrows=5, ncols=5, figure=figure)
    fig, ax_box = plt.subplots(5, 5)

    row = column = 0
    for specie in alpha_diversities:
        #ax_box = figure.add_subplot(spec[row, column])
        ax_box[row * 5 + column].boxplot([alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity']], labels=['Wild', 'Captivity'])
        #sns.boxplot(ax=ax_box[row * 5 + column])

        ax_box[row * 5 + column].set_title(specie)

        if row == 2 and column == 0:
            ax_box.set_ylabel("Alpha diversity")
        if row == 4 and column == 2:
            ax_box.set_xlabel("Sample type")
        
        if column >= 4:
            column = 0
            row += 1
        else:
            column += 1

    plt.suptitle("Bacterial diversity in vertebrate species")
    plt.show()


def show_histogram(relative_abundances, individual):
    plt.hist(relative_abundances[individual], bins=[0, 0.0001, 0.001, 0.01, 0.1, 1])
    plt.xscale('log')
    
    plt.title("Relative abundances of bacterial genus in "+individual)
    plt.ylabel("Num bacterial genus.")
    plt.xlabel("Discretized relative abundances")
    
    plt.show()
