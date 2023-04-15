from matplotlib import pyplot as plt
from matplotlib import gridspec


def name_specie(specie, f_codes_vertebrates):
    for vertebrate_specie in f_codes_vertebrates:
        if specie == vertebrate_specie.split()[0]:
            print(specie+'\t'+vertebrate_specie.split()[1].replace('_', ' ', 1))
            break


def alpha_diversities_sample_type(alpha_diversities, sample_type):
    print('\n'+sample_type)
    for specie in alpha_diversities:
        print(specie, alpha_diversities[specie][sample_type])


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
        return "n.s"


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


def boxplot(alpha_diversities, t_tests):
    figure = plt.figure()
    spec = gridspec.GridSpec(nrows=5, ncols=5, figure=figure)

    row = column = 0
    for specie in alpha_diversities:
        ax_box = figure.add_subplot(spec[row, column])
        
        bp = ax_box.boxplot([alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity']], labels=['Wild ('+str(len(alpha_diversities[specie]['Wild']))+')', 'Captive ('+str(len(alpha_diversities[specie]['Captivity']))+')'])
        
        xl = (bp['caps'][1].get_xdata()[0] + bp['caps'][1].get_xdata()[1]) / 2
        xr = (bp['caps'][3].get_xdata()[0] + bp['caps'][3].get_xdata()[1]) / 2
        yd = min(bp['caps'][1].get_ydata()[0], bp['caps'][3].get_ydata()[0]) + 0.1
        yu = yd + 0.1
        plt.plot([xl, xl, xr, xr], [yd, yu, yu, yd], lw=1)

        y = min(bp['caps'][1].get_ydata()[0], bp['caps'][3].get_ydata()[0]) + 0.2
        ax_box.text(x=(xl + xr) / 2, y=y, s=significance_conversion(t_tests[specie][1]), ha='center', va='bottom')
        #ax_box.annotate(significance_conversion(t_tests[specie][1]), ((xl + xr) / 2, y))

        ax_box.set_title(specie)

        if row == int(spec.nrows / 2) and column == 0:
            ax_box.set_ylabel("Alpha diversity")
        if row == (spec.nrows - 1) and column == int(spec.ncols / 2):
            ax_box.set_xlabel("Sample type")
        
        if column >= 4:
            column = 0
            row += 1
        else:
            column += 1

    plt.suptitle("Bacterial genus diversity in vertebrate species")
    plt.show()
