import support_functions as support
from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import pandas as pd
import statannot
import abc


def alpha_diversities_sample_type(alpha_diversities, sample_type):
    print(sample_type)
    for specie in alpha_diversities:
        print(specie, alpha_diversities[specie][sample_type])
    print()


def alpha_diversities(alpha_diversities):
    alpha_diversities_sample_type(alpha_diversities, 'Wild')
    alpha_diversities_sample_type(alpha_diversities, 'Captivity')


class Ploter(abc.ABC):
    def get_name_specie(self, specie, name_file_codes_vertebrates):
        f_codes_vertebrates = open(name_file_codes_vertebrates, 'r')
        
        for vertebrate_specie in f_codes_vertebrates:
            if specie == vertebrate_specie.split()[0]:
                return vertebrate_specie.split()[1].replace('_', ' ', 1)
        f_codes_vertebrates.close()

    def significance_conversion(self, p_value):
        if p_value < 0.001:
            return "***"
        elif p_value < 0.01:
            return "**"
        elif p_value < 0.05:
            return "*"
        else:
            return "n.s."
    
    @abc.abstractmethod
    def histogram(self, relative_abundances):
        pass

    @abc.abstractmethod
    def boxplot(self, alpha_diversities, name_file_codes_vertebrates):
        pass


class PyplotPloter(Ploter):
    def histogram(self, relative_abundances):
        plt.hist(x=relative_abundances, bins=[0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])

        plt.xscale('log')
        
        plt.title("Relative diversities of bacterial genus")
        
        plt.ylabel("Num bacterial genus")
        plt.xlabel("Relative diversities")
        
        plt.show()

    def boxplot(self, alpha_diversities, name_file_codes_vertebrates):
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
            
            significance = self.significance_conversion(alpha_diversities[specie]['p_value'])
            if significance == 'n.s.':
                y = yu + yrange / 2
            else:
                y = yu - yrange / 2
            ax_box.text(x=(xl + xr) / 2, y=y, s=significance, fontsize=7)

            ax_box.tick_params(axis='x', labelsize=8)
            ax_box.tick_params(axis='y', labelsize=8)

            ax_box.set_title(self.get_name_specie(specie, name_file_codes_vertebrates), fontsize=9, y=0.95)
            
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


class SeabornPloter(Ploter):
    def histogram(self, relative_abundances):
        ax_hist = sns.histplot(data=relative_abundances, bins=[0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
        
        ax_hist.set_xscale('log')

        ax_hist.set_title("Relative diversities of bacterial genus")
        
        ax_hist.set_ylabel("Num bacterial genus")
        ax_hist.set_xlabel("Relative diversities")
        
        plt.show()    

    def boxplot(self, alpha_diversities, name_file_codes_vertebrates, option='manual'):
        figure, axes = plt.subplots(nrows=5, ncols=5)
        plt.subplots_adjust(hspace=0.5)

        row = column = 0
        for specie in alpha_diversities:
            ax_box = axes[row, column]

            wild = 'Wild ('+str(len(alpha_diversities[specie]['Wild']))+')'
            captive = 'Captive ('+str(len(alpha_diversities[specie]['Captivity']))+')'
            
            data = {wild: alpha_diversities[specie]['Wild'], captive: alpha_diversities[specie]['Captivity']}
            support.pad_array(data[wild], data[captive])
            data_df = pd.DataFrame(data)
            
            sns.boxplot(data=data_df, ax=ax_box, width=0.25)
            ax_box.set_ylim(0.0, 5.1)

            if (option == 'manual'):
                xl = ax_box.lines[3].get_xdata().mean()
                xr = ax_box.lines[9].get_xdata().mean()
                yrange = (ax_box.get_ylim()[1] - ax_box.get_ylim()[0]) * 0.04
                yd = max(ax_box.lines[3].get_ydata()[0], ax_box.lines[9].get_ydata()[0]) + yrange
                yu = yd + yrange
                ax_box.plot([xl, xl, xr, xr], [yd, yu, yu, yd], lw=1, c='k')

                significance = self.significance_conversion(alpha_diversities[specie]['p_value'])
                if significance == 'n.s.':
                    y = yu + yrange / 2
                else:
                    y = yu - yrange / 2
                ax_box.text(x=(xl + xr) / 2, y=y, s=significance, fontsize=7)
            else:
                statannot.add_stat_annotation(ax=ax_box, data=data_df, box_pairs=[(wild, captive)], test='t-test_ind', text_format='star')
            
            ax_box.tick_params(axis='x', labelsize=8)
            ax_box.tick_params(axis='y', labelsize=8)

            ax_box.set_title(self.get_name_specie(specie, name_file_codes_vertebrates), fontsize=9, y=0.95)

            if row == int(axes.shape[0] / 2) and column == 0:
                ax_box.set_ylabel("Alpha diversity", fontsize=11, labelpad=10)
            if row == (axes.shape[0] - 1) and column == int(axes.shape[1] / 2):
                ax_box.set_xlabel("Sample type", fontsize=11, labelpad=10)

            if column >= 4:
                column = 0
                row += 1
            else:
                column += 1

        figure.suptitle("Bacterial genus diversity in vertebrate species")
        plt.show()
