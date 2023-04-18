from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns


def name_specie(specie, name_file_codes_vertebrates):
    f_codes_vertebrates = open(name_file_codes_vertebrates, 'r')
    
    for vertebrate_specie in f_codes_vertebrates:
        if specie == vertebrate_specie.split()[0]:
            return vertebrate_specie.split()[1].replace('_', ' ', 1)
    f_codes_vertebrates.close()


def alpha_diversities_sample_type(alpha_diversities, sample_type):
    print(sample_type)
    for specie in alpha_diversities:
        print(specie, alpha_diversities[specie][sample_type])
    print()


def alpha_diversities(alpha_diversities):
    alpha_diversities_sample_type(alpha_diversities, 'Wild')
    alpha_diversities_sample_type(alpha_diversities, 'Captivity')


def significance_conversion(p_value):
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "n.s."


def t_test(t_tests):
    print("\nStudent's t-test")
    for specie in t_tests:
        print(specie+":\tt statistic = "+str(round(t_tests[specie][0], 10))+"\tp-value = "+str(round(t_tests[specie][1], 10))+"\t"+significance_conversion(t_tests[specie][1]))


def histogram(relative_abundances):
    plt.hist(relative_abundances, bins=[0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
    
    plt.xscale('log')
    
    plt.title("Relative diversities of bacterial genus")
    plt.ylabel("Num bacterial genus")
    plt.xlabel("Relative diversities")
    
    plt.show()


def pyplot_boxplot(alpha_diversities, name_file_codes_vertebrates):
    figure = plt.figure()
    spec = gridspec.GridSpec(nrows=5, ncols=5, figure=figure)
    spec.update(hspace=0.5)

    row = column = 0
    for specie in alpha_diversities:
        ax_box = figure.add_subplot(spec[row, column])
        
        bp = ax_box.boxplot([alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity']], labels=['Wild ('+str(len(alpha_diversities[specie]['Wild']))+')', 'Captive ('+str(len(alpha_diversities[specie]['Captivity']))+')']) 
        ax_box.set_ylim(0.0, 5.1)

        xl = (bp['caps'][1].get_xdata()[0] + bp['caps'][1].get_xdata()[1]) / 2
        xr = (bp['caps'][3].get_xdata()[0] + bp['caps'][3].get_xdata()[1]) / 2
        yrange = (ax_box.get_ylim()[1] - ax_box.get_ylim()[0]) * 0.04
        yd = max(bp['caps'][1].get_ydata()[0], bp['caps'][3].get_ydata()[0]) + yrange
        yu = yd + yrange
        ax_box.plot([xl, xl, xr, xr], [yd, yu, yu, yd], lw=1, c='k')
        
        significance = significance_conversion(alpha_diversities[specie]['p_value'])
        if significance == 'n.s.':
            y = yu + yrange / 2
        else:
            y = yu - yrange / 2
        ax_box.text(x=(xl + xr) / 2, y=y, s=significance, fontsize=7)

        ax_box.tick_params(axis='x', labelsize=8)
        ax_box.tick_params(axis='y', labelsize=8)

        ax_box.set_title(name_specie(specie, name_file_codes_vertebrates), fontsize=9, y=0.95)
        
        if row == int(spec.nrows / 2) and column == 0:
            ax_box.set_ylabel("Alpha diversity", fontsize=11, labelpad=10)
        if row == (spec.nrows - 1) and column == int(spec.ncols / 2):
            ax_box.set_xlabel("Sample type", fontsize=11, labelpad=10)
        
        if column >= 4:
            column = 0
            row += 1
        else:
            column += 1

    plt.suptitle("Bacterial genus diversity in vertebrate species")
    plt.show()


def seaborn_boxplot(alpha_diversities, f_codes_vertebrates):
    figure, axes = plt.subplots(nrows=5, ncols=5)
    plt.subplots_adjust(hspace=0.5)

    row = column = 0
    for specie in alpha_diversities:
        ax_box = axes[row, column]

        #sns.boxplot()

        ax_box.tick_params(axis='x', labelsize=8)
        ax_box.tick_params(axis='y', labelsize=8)

        ax_box.set_title(name_specie(specie, f_codes_vertebrates), fontsize=9, y=0.95)

        if row == int(axes.shape[0] / 2) and column == 0:
            ax_box.set_ylabel("Alpha diversity")
        if row == (axes.shape[0] - 1) and column == int(axes.shape[1] / 2):
            ax_box.set_xlabel("Sample type")

        if column >= 4:
            column = 0
            row += 1
        else:
            column += 1

    figure.suptitle("Bacterial genus diversity in vertebrate species")
    plt.show()
