import calculation_functions as calculation
import support_functions as support
from matplotlib import pyplot as plt
from matplotlib import gridspec
import scipy.cluster.hierarchy as sch
import seaborn as sns
import pandas as pd
import statannot
import abc


def select_ploter(library):
    if library == "matplotlib" or library == "pyplot":
        return PyplotPloter()
    return SeabornPloter()


class Ploter(abc.ABC):
    def alpha_diversities_sample_type(self, alpha_diversities_individual, sample_type):
        print(sample_type)
        for specie in alpha_diversities_individual:
            print(specie, [round(alpha_diversity, 2)
                           for alpha_diversity in alpha_diversities_individual[specie][sample_type]])

    def alpha_diversities(self, alpha_diversities_individual):
        self.alpha_diversities_sample_type(alpha_diversities_individual, 'Wild')
        print()
        self.alpha_diversities_sample_type(alpha_diversities_individual, 'Captivity')

    def stemplot(self, average_distances, name_file_code_vertebrates):
        data, labels = support.generate_plot_data(average_distances, name_file_code_vertebrates)

        figure = plt.figure(figsize=(11.5, 8.5))
        figure.subplots_adjust(bottom=0.2)

        ax_stem = figure.add_subplot()

        ax_stem.stem(data)

        ax_stem.set_xticks(range(len(labels)))
        ax_stem.set_xticklabels(labels)
        ax_stem.set_xticklabels(ax_stem.get_xticklabels(), rotation=90)

        ax_stem.tick_params(axis='x', labelsize=8)
        ax_stem.tick_params(axis='y', labelsize=8)

        ax_stem.set_ylabel('Average distances', fontsize=13, labelpad=12)
        ax_stem.set_xlabel('Vertebrate species', fontsize=13, labelpad=12)

        plt.suptitle('Average distance between wild, captive, and wild-captive individuals in every vertebrate specie',
                     fontsize=14)
        plt.show()

    def scatterplot(self, alpha_average, distance_average, name_file_code_vertebrates):
        average_alpha, average_distance, labels = support.generate_scatterplot_data(alpha_average, distance_average,
                                                                                    name_file_code_vertebrates)
        self.initialize_figure()
        self.figure.subplots_adjust(bottom=0.2)

        ax_scatter = self.initialize_plot()

        self.set_scatterplot(ax_scatter, average_alpha, average_distance, labels)
        plt.setp(self.legend.get_texts(), fontsize=8)

        pearson_corr, p_value = calculation.pearson_correlation(average_alpha, average_distance)
        spearman_corr, p_value = calculation.spearman_correlation(average_alpha, average_distance)
        ax_scatter.text(x=0.25, y=0.1, s=f'Pearson correlation: {pearson_corr:.3f}\np-value: {p_value:.2e}\n\n'
                                         f'Spearman correlation: {spearman_corr:.3f}\np-value: {p_value:.2e}')

        ax_scatter.tick_params(axis='x', labelsize=8)
        ax_scatter.tick_params(axis='y', labelsize=8)

        ax_scatter.set_ylabel('Average distances', fontsize=11, labelpad=10)
        ax_scatter.set_xlabel('Average alpha diversities', fontsize=11, labelpad=10)

        self.set_suptitle('Correlation between distances and alpha diversities in vertebrate species')
        plt.show()

    @abc.abstractmethod
    def set_scatterplot(self, ax_scatter, average_alpha, average_distance, labels):
        pass

    def histogram(self, relative_abundances):
        ax_hist = self.set_histogram(relative_abundances)
        
        ax_hist.set_xscale('log')

        ax_hist.set_ylabel("Num bacterial genus", fontsize=13, labelpad=12)
        ax_hist.set_xlabel("Relative diversities", fontsize=13, labelpad=12)

        self.set_suptitle("Relative diversities of bacterial genus in individuals of vertebrates species")
        plt.show()

    @abc.abstractmethod
    def set_histogram(self, relative_abundances):
        pass

    def boxplot(self, vertebrates_distances, base_line, name_file_code_vertebrates):
        data, labels = support.generate_plot_data(vertebrates_distances, name_file_code_vertebrates)

        self.initialize_figure()
        self.figure.subplots_adjust(bottom=0.2)

        ax_box = self.initialize_plot()

        self.set_boxplot(ax_box, data, labels)

        ax_box.axhline(y=base_line, color='black')

        ax_box.tick_params(axis='x', labelsize=8)
        ax_box.tick_params(axis='y', labelsize=8)

        ax_box.set_xticklabels(ax_box.get_xticklabels(), rotation='vertical')

        ax_box.set_ylabel('Distances', fontsize=13, labelpad=12)
        ax_box.set_xlabel('Vertebrate species', fontsize=13, labelpad=12)

        self.set_suptitle('Distances between wild, captive, and wild-captive individuals in every vertebrate specie')
        plt.show()

    @abc.abstractmethod
    def initialize_figure(self):
        pass

    @abc.abstractmethod
    def initialize_plot(self):
        pass

    @abc.abstractmethod
    def set_boxplot(self, ax_box, data, labels):
        pass

    def boxplot_grid(self, alpha_diversities, name_file_codes_vertebrates, mechanism='manual'):
        self.initialize_grid()
        
        row = column = 0
        for specie in alpha_diversities:
            ax_box = self.initialize_subplot(row, column)

            self.wild = 'Wild ('+str(len(alpha_diversities[specie]['Wild']))+')'
            self.captive = 'Captive ('+str(len(alpha_diversities[specie]['Captivity']))+')'

            self.set_boxplot_grid(ax_box, alpha_diversities, specie, self.wild, self.captive)
            ax_box.set_ylim(0.0, 5.1)
    
            self.mechanism = mechanism
            self.select_mechanism(ax_box, alpha_diversities, specie)

            ax_box.tick_params(axis='x', labelsize=8)
            ax_box.tick_params(axis='y', labelsize=8)

            ax_box.set_title(support.get_name_specie(specie, name_file_codes_vertebrates), fontsize=9, y=0.95)

            if row == int(self.get_nrows() / 2) and column == 0:
                ax_box.set_ylabel("Alpha diversity", fontsize=13, labelpad=12)
            if row == (self.get_nrows() - 1) and column == int(self.get_ncols() / 2):
                ax_box.set_xlabel("Sample type", fontsize=13, labelpad=12)

            if column >= 4:
                column = 0
                row += 1
            else:
                column += 1
        
        self.set_suptitle("Alpha diversities of bacterial genus in vertebrate species")
        plt.show()

    @abc.abstractmethod
    def initialize_grid(self):     
        pass 

    @abc.abstractmethod
    def initialize_subplot(self, row, column):
        pass

    @abc.abstractmethod
    def set_boxplot_grid(self, ax_box, alpha_diversities, specie, wild, captive):
        pass

    @abc.abstractmethod
    def select_mechanism(self, ax_box, alpha_diversities, specie):
        pass

    def significance_conversion(self, p_value):
        if p_value < 0.001:
            return "***"
        elif p_value < 0.01:
            return "**"
        elif p_value < 0.05:
            return "*"
        else:
            return "n.s."

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
    def set_suptitle(self, title):
        pass

    def dendrogram(self, matrix, label):
        sch.dendrogram(calculation.hierarchical_clustering(matrix))

        plt.ylabel('Distance between ' + label, fontsize=13, labelpad=12)
        plt.xlabel(label[0].upper() + label[1:], fontsize=13, labelpad=12)

        plt.suptitle('Dendrogram of the ' + label + ' matrix', fontsize=14)
        plt.show()

    def heatmap(self, matrix):
        row_order = sch.leaves_list(calculation.hierarchical_clustering(matrix)) - 1
        column_order = sch.leaves_list(calculation.hierarchical_clustering(matrix[:, row_order].T)) - 1

        reordered_matrix = matrix[row_order, :]
        reordered_matrix = reordered_matrix[:, column_order]

        self.set_heatmap(reordered_matrix)
        plt.show()

    @abc.abstractmethod
    def set_heatmap(self, reordered_matrix):
        pass

    def cluster_map(self, matrix, colour_map):
        self.set_cluster_map(matrix, colour_map)
        plt.show()

    @abc.abstractmethod
    def set_cluster_map(self, matrix, colour_map):
        pass


class PyplotPloter(Ploter):
    def set_scatterplot(self, ax_scatter, average_alpha, average_distance, labels):
        colors = support.create_colors()
        sp = []

        i = 0
        while i < len(average_alpha):
            sp.append(ax_scatter.scatter(x=average_alpha[i], y=average_distance[i], marker='o', s=20,
                                         color=colors[int(i / 2)]))
            sp.append(ax_scatter.scatter(x=average_alpha[i + 1], y=average_distance[i + 1], marker='s', s=15,
                                         color=colors[int(i / 2)]))
            i += 2

        self.legend = ax_scatter.legend(handles=sp[::2], labels=labels, loc='upper right', ncol=3)

    def set_histogram(self, relative_abundances):
        ax_hist = plt.figure(figsize=(11, 8.5)).add_subplot()
        ax_hist.hist(x=relative_abundances, bins=[0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
        return ax_hist

    def initialize_figure(self):
        self.figure = plt.figure(figsize=(11.5, 8.5))

    def initialize_plot(self):
        return self.figure.add_subplot()

    def set_boxplot(self, ax_box, data, labels):
        ax_box.boxplot(data, labels=labels)

    def initialize_grid(self):
        self.figure = plt.figure(figsize=(11.5, 8.5))
        self.spec = gridspec.GridSpec(nrows=5, ncols=5, figure=self.figure)
        self.spec.update(hspace=0.5)

    def initialize_subplot(self, row, column):
        return self.figure.add_subplot(self.spec[row, column])

    def set_boxplot_grid(self, ax_box, alpha_diversities, specie, wild, captive):
        self.bp = ax_box.boxplot([alpha_diversities[specie]['Wild'], alpha_diversities[specie]['Captivity']],
                                 labels=[wild, captive])

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

    def set_suptitle(self, title):
        plt.suptitle(title, fontsize=14)

    def set_heatmap(self, reordered_matrix):
        plt.imshow(reordered_matrix)
        plt.colorbar()

    def set_cluster_map(self, matrix, colour_map):
        matrix = calculation.hierarchical_clustering(matrix)

        figure = plt.figure()
        spec = gridspec.GridSpec(nrows=1, ncols=2, figure=figure)

        ax_dendrogram = figure.add_subplot(spec[0, 0])
        sch.dendrogram(matrix, orientation='left')
        ax_dendrogram.axis('off')

        ax_heatmap = figure.add_subplot(spec[0, 1])
        heatmap = ax_heatmap.imshow(matrix)
        plt.colorbar(heatmap)


class SeabornPloter(Ploter):
    def set_scatterplot(self, ax_scatter, average_alpha, average_distance, labels):
        sns.scatterplot(x=average_alpha[::2], y=average_distance[::2], ax=ax_scatter, hue=labels,
                        marker='o', s=30, palette='husl')
        sns.scatterplot(x=average_alpha[1::2], y=average_distance[1::2], ax=ax_scatter, hue=labels,
                        marker='s', s=25, palette='husl')

        self.legend = ax_scatter.legend(labels=labels, loc='upper right', ncol=3, markerscale=0.65)

    def set_histogram(self, relative_abundances):
        self.figure = plt.figure(figsize=(11, 8.5))
        ax_hist = sns.histplot(data=relative_abundances,
                               bins=[0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
        return ax_hist

    def initialize_figure(self):
        self.figure, self.axes = plt.subplots(nrows=1, ncols=1, figsize=(11.5, 8.5))

    def initialize_plot(self):
        return self.axes

    def set_boxplot(self, ax_box, data, labels):
        data_df = []

        num_data = 0
        while num_data < len(data):
            data_df.append(pd.DataFrame({labels[num_data]: data[num_data]}))
            num_data += 1

        sns.boxplot(data=pd.concat(data_df), ax=ax_box)

    def initialize_grid(self):
        self.figure, self.axes = plt.subplots(nrows=5, ncols=5, figsize=(11.5, 8.5))
        plt.subplots_adjust(hspace=0.5)

    def initialize_subplot(self, row, column):
        return self.axes[row, column]

    def set_boxplot_grid(self, ax_box, alpha_diversities, specie, wild, captive):
        data = {wild: (support.pad_list_average(alpha_diversities[specie]['Wild'].copy(),
                                                alpha_diversities[specie]['Captivity'].copy()))[0],
                captive: (support.pad_list_average(alpha_diversities[specie]['Wild'].copy(),
                                                   alpha_diversities[specie]['Captivity'].copy()))[1]}
        self.data_df = pd.DataFrame(data)

        sns.boxplot(data=self.data_df, ax=ax_box, width=0.25)
    
    def select_mechanism(self, ax_box, alpha_diversities, specie):
        if self.mechanism == 'automatic':
            statannot.add_stat_annotation(ax=ax_box, data=self.data_df, box_pairs=[(self.wild, self.captive)],
                                          test='t-test_ind', text_format='star')
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

    def set_suptitle(self, title):
        self.figure.suptitle(title, fontsize=14)

    def set_heatmap(self, reordered_matrix):
        sns.heatmap(reordered_matrix)

    def set_cluster_map(self, matrix, colour_map):
        if not colour_map:
            sns.clustermap(matrix, method='average', metric='euclidean')
        else:
            sns.clustermap(matrix, cmap=colour_map, method='average', metric='euclidean', center=0)
