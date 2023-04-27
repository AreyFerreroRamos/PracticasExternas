import calculation_functions as calculation
import support_functions as support
from matplotlib import pyplot as plt
from matplotlib import gridspec
import scipy.cluster.hierarchy as sch
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


def select_ploter(library):
    if library == "matplotlib" or library == "pyplot":
        return PyplotPloter()
    return SeabornPloter()


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

    def histogram(self, relative_abundances):
        ax_hist = self.set_histogram(relative_abundances)
        
        ax_hist.set_xscale('log')
        
        ax_hist.set_title("Relative diversities of bacterial genus")
        
        ax_hist.set_ylabel("Num bacterial genus")
        ax_hist.set_xlabel("Relative diversities")
        
        plt.show()

    @abc.abstractmethod
    def set_histogram(self, relative_abundances):
        pass

    def boxplot(self, alpha_diversities, name_file_codes_vertebrates, mechanism='manual'):
        self.initialize_grid()
        
        row = column = 0
        for specie in alpha_diversities:
            ax_box = self.initialize_plot(row, column)

            self.wild = 'Wild ('+str(len(alpha_diversities[specie]['Wild']))+')'
            self.captive = 'Captive ('+str(len(alpha_diversities[specie]['Captivity']))+')'

            self.set_boxplot(ax_box, alpha_diversities, specie, self.wild, self.captive)
            ax_box.set_ylim(0.0, 5.1)
    
            self.mechanism = mechanism
            self.select_mechanism(ax_box, alpha_diversities, specie)

            ax_box.tick_params(axis='x', labelsize=8)
            ax_box.tick_params(axis='y', labelsize=8)

            ax_box.set_title(self.get_name_specie(specie, name_file_codes_vertebrates), fontsize=9, y=0.95)

            if row == int(self.get_nrows() / 2) and column == 0:
                ax_box.set_ylabel("Alpha diversity", fontsize=11, labelpad=10)
            if row == (self.get_nrows() - 1) and column == int(self.get_ncols() / 2):
                ax_box.set_xlabel("Sample type", fontsize=11, labelpad=10)

            if column >= 4:
                column = 0
                row += 1
            else:
                column += 1
        
        self.set_suptitle()
        plt.show()

    @abc.abstractmethod
    def initialize_grid(self):     
        pass 

    @abc.abstractmethod
    def initialize_plot(self, row, column):
        pass

    @abc.abstractmethod
    def set_boxplot(self, ax_box, alpha_diversities, specie, wild, captive):
        pass

    @abc.abstractmethod
    def select_mechanism(self, ax_box, alpha_diversities, specie):
        pass

    def set_significance(self, ax_box, alpha_diversities, specie):
        yrange = (ax_box.get_ylim()[1] - ax_box.get_ylim()[0]) * 0.04
        xl, xr, yd, yu = self.set_line_significance(ax_box, yrange)

        significance = self.significance_conversion(alpha_diversities[specie]['p_value'])
        if significance == 'n.s.':
            y = yu + yrange / 2
        else:
            y = yu - yrange / 2

        ax_box.plot([xl, xl, xr, xr], [yd, yu, yu, yd], lw=1, c='k')
        ax_box.text(x=(xl + xr) / 2, y=y, s=significance, fontsize=7)

    @abc.abstractmethod
    def set_line_significance(self, ax_box, yrange):
        pass

    @abc.abstractmethod
    def get_nrows(self):
        pass

    @abc.abstractmethod
    def get_ncols(self):
        pass

    @abc.abstractmethod
    def set_suptitle(self):
        pass

    def dendrogram(self, matrix, x_label):
        sch.dendrogram(calculation.hierarchical_clustering(matrix))

        plt.title('Dendrogram')
        plt.ylabel('Distance')
        plt.xlabel(x_label)

        plt.show()

    def heatmap(self, matrix):
        self.set_heatmap(matrix)
        plt.show()

    @abc.abstractmethod
    def set_heatmap(self, matrix):
        pass

    def cluster_map(self, matrix):
        self.set_cluster_map(matrix)
        plt.show()

    @abc.abstractmethod
    def set_cluster_map(self, matrix):
        pass


class PyplotPloter(Ploter):
    def set_histogram(self, relative_abundances):
        ax_hist = plt.figure().add_subplot()
        ax_hist.hist(x=relative_abundances, bins=[0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
        return ax_hist

    def initialize_grid(self):
        self.figure = plt.figure()
        self.spec = gridspec.GridSpec(nrows=5, ncols=5, figure=self.figure)
        self.spec.update(hspace=0.5)

    def initialize_plot(self, row, column):
        return self.figure.add_subplot(self.spec[row, column])

    def set_boxplot(self, ax_box, alpha_diversities, specie, wild, captive):
        self.bp = ax_box.boxplot([alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity']], labels=[wild, captive])

    def select_mechanism(self, ax_box, alpha_diversities, specie):
        self.set_significance(ax_box, alpha_diversities, specie)

    def set_line_significance(self, ax_box, yrange):
        xl = (self.bp['caps'][1].get_xdata()[0] + self.bp['caps'][1].get_xdata()[1]) / 2
        xr = (self.bp['caps'][3].get_xdata()[0] + self.bp['caps'][3].get_xdata()[1]) / 2
        yd = max(self.bp['caps'][1].get_ydata()[0], self.bp['caps'][3].get_ydata()[0]) + yrange
        yu = yd + yrange

        return xl, xr, yd, yu

    def get_nrows(self):
        return self.spec.nrows

    def get_ncols(self):
        return self.spec.ncols

    def set_suptitle(self):
        plt.suptitle("Bacterial genus diversity in vertebrate species")

    def set_heatmap(self, matrix):
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()

    def set_cluster_map(self, matrix):
        matrix = calculation.hierarchical_clustering(matrix)

        figure = plt.figure()
        spec = gridspec.GridSpec(nrows=1, ncols=2, figure=figure)

        ax_dendrogram = figure.add_subplot(spec[0, 0])
        sch.dendrogram(matrix, orientation='left')
        ax_dendrogram.axis('off')

        ax_heatmap = figure.add_subplot(spec[0, 1])
        heatmap = ax_heatmap.imshow(matrix, cmap='viridis')
        plt.colorbar(heatmap)


class SeabornPloter(Ploter):
    def set_histogram(self, relative_abundances):
        ax_hist = sns.histplot(data=relative_abundances, bins=[0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
        return ax_hist

    def initialize_grid(self):
        self.figure, self.axes = plt.subplots(nrows=5, ncols=5)
        plt.subplots_adjust(hspace=0.5)

    def initialize_plot(self, row, column):
        return self.axes[row, column]

    def set_boxplot(self, ax_box, alpha_diversities, specie, wild, captive):
        data = {wild: alpha_diversities[specie]['Wild'], captive: alpha_diversities[specie]['Captivity']}
        support.pad_array(data[wild], data[captive])
        self.data_df = pd.DataFrame(data)

        sns.boxplot(data=self.data_df, ax=ax_box, width=0.25)
    
    def select_mechanism(self, ax_box, alpha_diversities, specie):
        if self.mechanism == 'automatic':
            statannot.add_stat_annotation(ax=ax_box, data=self.data_df, box_pairs=[(self.wild, self.captive)], test='t-test_ind', text_format='star')
        else:
            self.set_significance(ax_box, alpha_diversities, specie)

    def set_line_significance(self, ax_box, yrange):
        xl = ax_box.lines[3].get_xdata().mean()
        xr = ax_box.lines[9].get_xdata().mean()
        yd = max(ax_box.lines[3].get_ydata()[0], ax_box.lines[9].get_ydata()[0]) + yrange
        yu = yd + yrange

        return xl, xr, yd, yu

    def get_nrows(self):
        return self.axes.shape[0]

    def get_ncols(self):
        return self.axes.shape[1]

    def set_suptitle(self):
        self.figure.suptitle("Bacterial genus diversity in vertebrate species")

    def set_heatmap(self, matrix):
        row_order = sch.leaves_list(calculation.hierarchical_clustering(matrix)) - 1
        column_order = sch.leaves_list(calculation.hierarchical_clustering(matrix[:, row_order].T)) - 1

        reordered_matrix = matrix[row_order, :]
        reordered_matrix = reordered_matrix[:, column_order]

        sns.heatmap(reordered_matrix, cmap='viridis')

    def set_cluster_map(self, matrix):
        sns.clustermap(matrix, cmap='viridis', method='average', metric='euclidean')
