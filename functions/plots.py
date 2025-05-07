import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write
import sys

sys.path.append("/functions")
from functions.metric import calculate_rmse, calculate_r2


class PhononPlotter:
    def __init__(
        self,
        distances_set: list,
        frequencies_set: list,
        x_labels: list,
        connections: list,
        colors=None,
        legend_labels=None,
        linestyles=None,
        linewidths=None,
        figsize=(10, 6),
        commensurate_points=None,
        meterial_name=None,
    ):
        """
        initializing the phonon plotter
        Parameters:
        distances_set: list of lists of distances for each phonon band structure
        frequencies_set: list of lists of frequencies for each phonon band structure
        x_labels: list of labels for the x-axis
        connections: list of booleans indicating if the x-path is connected or not
        colors: list of colors for each phonon band structure
        legend_labels: list of labels for the legend
        linestyles: list of linestyles for each phonon band structure
        linewidths: list of linewidths for each phonon band structure
        figsize: tuple of figure size
        commensurate_points: list of commensurate points for the x-axis
        meterial_name: name of the material
        """
        self.distances_set = distances_set
        self.frequencies_set = frequencies_set
        self.x_labels = x_labels
        self.connections = connections
        self.colors = (
            colors
            if colors is not None
            else ["black", "red", "blue", "green", "orange"]
        )
        self.legend_labels = (
            legend_labels
            if legend_labels is not None
            else ["DFT", "MP-0", "DFT+MP0", "ft MACE", "DFT+ft MACE"]
        )
        self.linestyles = (
            linestyles if linestyles is not None else ["-", "-", "-", ":", ":"]
        )
        self.linewidths = (
            linewidths if linewidths is not None else [1] * len(self.colors)
        )
        self.figsize = figsize
        self.commensurate_points = commensurate_points
        self.material_name = meterial_name

    def _create_xtick_labels(self):
        """create xtick labels for the phonon band structure plot
        Depending on wehter the x-path is connected or not, the labels have to repeat or change
        An example:
        x_labels = [Gamma,X,K,L, Gamma]
        connection = [True, True, False]
        '--------'---------'-----------'
        G        X        K| L         G

        x_labels: list of labels for example: [Gamma,X,K, Gamma]
        connection: list with True and False label if x-path is connected or not. for example: [True, True, False, True]
        """
        xtick_labels = []
        if False not in self.connections:
            return self.x_labels
        else:
            xtick_labels.append(self.x_labels[0])
            count = 1

            for connection in self.connections:
                if count >= len(self.x_labels) - 1:
                    xtick_labels.append(self.x_labels[-1])
                    break
                if connection:
                    xtick_labels.append(self.x_labels[count])
                    count += 1
                else:
                    xtick_labels.append(
                        str(self.x_labels[count]) + " | " + self.x_labels[count + 1]
                    )
                    count += 2
        return xtick_labels

    def _plot_phonon_bands(
        self, ax, distances, frequencies, color, linestyle, linewidth
    ):
        """plots phonon band structrue"""
        num_paths, num_kpoints, num_bands = frequencies.shape
        for path in range(num_paths):
            for band in range(num_bands):
                ax.plot(
                    distances[path],
                    frequencies[path, :, band],
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                )

    def _add_vertical_lines_and_commensurate_points(self, ax, distances, comm_points):
        """add vertical lines and commensurate points to the plot"""
        xticks = []
        num_paths = distances.shape[0]
        for path in range(num_paths):
            start, end = distances[path][0], distances[path][-1]
            ax.axvline(x=start, color="black", linestyle="--", linewidth=0.5)
            ax.axvline(x=end, color="black", linestyle="--", linewidth=0.5)
            xticks.extend([start, end])
            if comm_points is not None and len(comm_points) > 0:
                for point in comm_points[path]:
                    ax.plot(distances[path][point], 0, color="red", marker="D")
        return sorted(set(xticks))

    def beautiful_phonon_plotter(self):
        """creates a beautiful phonon plot
        gives back:
        fig, ax: figure and axis of the plot
        """
        fig, ax = plt.subplots(figsize=self.figsize)

        # plot all data
        for i in range(len(self.distances_set)):
            self._plot_phonon_bands(
                ax,
                np.array(self.distances_set)[i],
                np.array(self.frequencies_set)[i],
                self.colors[i],
                self.linestyles[i],
                self.linewidths[i],
            )
            # empty plot for legend
            ax.plot(
                [],
                [],
                color=self.colors[i],
                linestyle=self.linestyles[i],
                label=self.legend_labels[i],
            )

        # add vertcal lines and commensurate points
        default_comm_points = (
            self.commensurate_points
            if self.commensurate_points is not None
            else [[]] * np.array(self.distances_set)[0].shape[0]
        )
        xticks = self._add_vertical_lines_and_commensurate_points(
            ax, np.array(self.distances_set)[0], default_comm_points
        )

        # set xtick labels
        xtick_labels = self._create_xtick_labels()
        if len(xticks) == len(xtick_labels):
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels)

        # additional plot settings
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)
        ax.tick_params(axis="x", direction="in")
        ax.tick_params(axis="y", direction="in")
        ax.set_ylabel("Frequency [THz]")
        ax.set_xlabel("Wave vector")
        ax.set_xlim(0, np.array(self.distances_set)[0][-1, -1])
        ax.legend(loc="lower right")

        if self.material_name:
            plt.suptitle(self.material_name, y=0.95)
        plt.tight_layout()
        return fig, ax
