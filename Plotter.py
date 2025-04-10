import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
import numpy as np
from typing import List, Dict, Any
import matplotlib.colors as mcolors
from matplotlib.ticker import (FixedLocator, FixedFormatter)
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mplhep.plot import ErrorBarArtists

from Hist import Hist
import math


class Plotter:
    def __init__(self, experiment, base_output_dir, **kwargs):

        self.experiment = experiment
        hep.style.use(self.experiment)
        self.year = ''
        self.lumi = ''
        plt.ioff()  # interactive mode off; not to show figures on creation

        plt.rcParams['axes.linewidth'] = 2.0
        plt.rcParams['hatch.linewidth'] = 0.5

        # plt.rcParams.update({
        #     "text.usetex": True,
        #     # "font.family": "Helvetica",
        #     # "font.family": "Helvetica",
        #     # 'text.latex.preamble': r'\usepackage{newtxmath}',
        #     # https://stackoverflow.com/questions/64405910/matplotlib-font-common-font-with-latex
        #     #'text.latex.preamble':
        #     #    r"\usepackage{libertine} \usepackage[libertine]{newtxmath}",
        #     'font.size': 23
        # })

        self.rows = 1
        self.cols = 1

        self.fig = None
        self.axs = None
        self.current_axis = None

        # make directory to save plots
        self.base_output_dir = base_output_dir
        self.out_dir = ''

        # histogram
        self.hist_list = []
        self.hist_loc = []
        self.hist_as_stack = []
        self.hist_kwargs: List[Dict[str, Any]] = []  # kwargs for matplotlib

        # errorbar
        self.errorbar_list = []  # [(x_data, y_data)]
        self.errorbar_loc = []
        self.errorbar_kwargs: List[Dict[str, Any]] = []

        self.legend_handles = []
        self.legend_labels = []

        self.y_minimum = 0

    def reset(self):
        self.rows = 1
        self.cols = 1

        del self.fig
        del self.axs
        self.current_axis = None

        # histogram
        self.hist_list.clear()
        self.hist_loc.clear()
        self.hist_as_stack.clear()
        self.hist_kwargs.clear()

        # errorbar
        self.errorbar_list.clear()
        self.errorbar_loc.clear()
        self.errorbar_kwargs.clear()

        self.legend_handles.clear()
        self.legend_labels.clear()

        self.y_minimum = 0

    def set_experiment_label(self, label="Preliminary", location=(0, 0),
                             **kwargs):
        self.set_current_axis(location=location)
        self.current_axis.draw(self.current_axis.figure.canvas.get_renderer())

        is_data = True
        if label == "Simulation":
            is_data = False
            label = ""
        hep.cms.label(label, data=is_data, fontsize=20,
                      ax=self.current_axis, loc=0, pad=.0, **kwargs)

    def create_subplots(self, rows=1, cols=1,
                        figsize=(8, 8), left=0.1, right=0.95, bottom=0.1, top=0.9,
                        hspace=0.2, wspace=0.2, **gridspec_kw):

        # update self.rows and self.cols
        self.rows = rows
        self.cols = cols

        self.fig, self.axs = plt.subplots(self.rows,
                                          self.cols,
                                          figsize=figsize,
                                          gridspec_kw=gridspec_kw)
        plt.tight_layout()
        plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                            hspace=hspace, wspace=wspace)

    def set_current_axis(self, location=(0, 0)):
        self.current_axis = self.get_axis(location=location)

    def get_axis(self, location=(0, 0)):
        row, col = location
        if self.rows == 1 and self.cols == 1:
            return self.axs
        elif self.rows == 1 or self.cols == 1:
            if self.rows == 1:
                index = col
            else:
                index = row
            return self.axs[index]
        else:
            return self.axs[row][col]

    def add_hist(self, hist, location=(0, 0), as_stack=False, **kwargs):
        self.hist_list.append(hist)
        self.hist_loc.append(location)
        self.hist_kwargs.append(kwargs)
        self.hist_as_stack.append(as_stack)
        return len(self.hist_list)-1

    def add_comparison_pair(self, nominator_hist, denominator_hist,
                            location, ratio_location,
                            nominator_args, denominator_args):
        nominator_index = self.add_hist(nominator_hist, location=location, **nominator_args)
        denominator_index = self.add_hist(denominator_hist, location=location, **denominator_args)
        self.add_ratio_hist(nominator_index=nominator_index, denominator_index=denominator_index,
                            location=ratio_location, **nominator_args)

    def add_text(self, text, location=(0, 0), do_magic=True, **kwargs):
        self.set_current_axis(location=location)
        plt.rcParams['text.usetex'] = True
        at = AnchoredText(text, prop=dict(size=20), **kwargs)
        self.current_axis.add_artist(at)
        if do_magic:
            hep.plot.mpl_magic(self.current_axis)
        plt.rcParams['text.usetex'] = False

    def add_custom_axis_tick_labels(self, locates, labels, location=(0, 0)):
        self.set_current_axis(location=location)
        self.current_axis.set_xticks([])
        self.current_axis.xaxis.set_minor_locator(FixedLocator(locates))
        self.current_axis.xaxis.set_minor_formatter(FixedFormatter(labels))
        self.current_axis.tick_params(axis='x', which='both', rotation=35, labelsize=5)
        for label in self.current_axis.xaxis.get_minorticklabels():
            label.set_horizontalalignment('right')

    def draw_vlines(self, vlines, location=(0, 0), **kwargs):
        for line in vlines:
            self.get_axis(location=location).axvline(x=line, linestyle='--', linewidth=1)

    def add_errorbar(self, errorbar_data, location=(0, 0), **kwargs):
        self.errorbar_list.append(errorbar_data)
        self.errorbar_loc.append(location)
        self.errorbar_kwargs.append(kwargs)

    def _add_hists(self, index_to_sum):
        added_hist = None
        for index in index_to_sum:
            added_hist = self.hist_list[index] + added_hist
        return added_hist

    def adjust_y_scale(self, location=(0, 0)):
        self.set_current_axis(location)
        self.current_axis.set_ylim(ymin=self.y_minimum * 1e-2)

    def add_ratio_hist(self, nominator_index, denominator_index, location=(0, 0), **kwargs):
        # for nominator
        if isinstance(nominator_index, int):
            ratio_hist = self.hist_list[nominator_index] + None
        elif isinstance(nominator_index, list):
            # sum all histograms
            ratio_hist = self._add_hists(nominator_index)
        else:
            pass

        # for denominator
        if isinstance(denominator_index, int):
            ratio_hist = ratio_hist.divide(self.hist_list[denominator_index])
        elif isinstance(denominator_index, list):
            ratio_hist = ratio_hist.divide(self._add_hists(denominator_index))
        else:
            pass

        if 'color' in self.hist_kwargs[nominator_index]:
            kwargs['color'] = self.hist_kwargs[nominator_index]['color']
        if 'histtype' in self.hist_kwargs[nominator_index]:
            kwargs['histtype'] = self.hist_kwargs[nominator_index]['histtype']

        self.add_hist(ratio_hist, location=location, **kwargs)

    def draw_hist(self):
        # TODO handle stack
        bottom = 0
        for index, hist in enumerate(self.hist_list):
            values, bins, errors = hist.to_numpy()
            self.set_current_axis(location=self.hist_loc[index])
            self.y_minimum = np.min(values)
            # if 'label' in self.hist_kwargs[index]:
            #    print(self.hist_kwargs[index]['label'])
            #errors=False

            if not self.hist_as_stack[index]:
                # NOTE when histtype is errorbar or yerr is set(i.e., when ErrorBarArtists is used),
                # somehow the legend order changed
                # histtype options: 'step', 'fill', 'errorbar', 'bar', 'barstep', 'band'
                if 'yerr' in self.hist_kwargs[index]:
                    if not self.hist_kwargs[index]['yerr']:
                        errors = False  # only allow to disable yerr
                        del self.hist_kwargs[index]['yerr']
                artists = hep.histplot((values, bins), ax=self.current_axis, yerr=errors,
                                       **self.hist_kwargs[index])

                # Note according to histtype, legend_artist could be None
                # Need to properly handle each histtype to get legend_artists
                # there are two cases: artists[0] or artists[0].stairs
                if isinstance(artists[0], ErrorBarArtists) or artists[0].legend_artist is not None:
                    handle = artists[0]
                else:
                    handle = artists[0].stairs
            else:
                # stack
                _, _, patches = self.current_axis.hist(bins[:-1], bins, weights=values, histtype='bar',
                                                       bottom=bottom, **self.hist_kwargs[index])
                bottom += values
                handle = patches[0]

            label = self.hist_kwargs[index].get("label", None)
            # deside whether to show legend
            if label:
                self.legend_handles.append(handle)
                self.legend_labels.append(label)

            #print(self.legend_labels)
            #handles, labels = self.current_axis.get_legend_handles_labels()
            #print(labels)
            self.current_axis.set_xlim(bins[0], bins[-1])
        # write legends in reverse order, i.e., force the order of legend always reverse of the drawn order
        self.legend_handles.reverse()
        self.legend_labels.reverse()
        #print(self.legend_labels)

    def draw_matrix(self, rm_np, variable_name, **kwargs):

        self.set_current_axis((0,0))
        # hep.hist2dplot(rm_np, norm=mcolors.LogNorm(), ax=self.current_axis, **kwargs)
        hep.hist2dplot(rm_np, ax=self.current_axis, **kwargs)

        if (rm_np[1][-1] - rm_np[1][0] >= 5e2 and
                rm_np[2][-1] - rm_np[2][0] >= 5e2):
            self.current_axis.set_xscale("log")
            self.current_axis.set_yscale("log")

            # self.current_axis.set_xlim(0.2, rm_np[1][-1])
            # self.current_axis.set_ylim(0.2, rm_np[2][-1])

        for x in range(len(rm_np[1]) - 1):
            for y in range(len(rm_np[2]) - 1):
                c = rm_np[0][x][y]
                if math.isnan(c):
                    continue
                if c < 0.1:
                    continue
                x_half_width = (rm_np[1][x+1] - rm_np[1][x])/2
                y_half_width = (rm_np[2][y+1] - rm_np[2][y])/2
                self.current_axis.text(rm_np[1][x] + x_half_width,
                                       rm_np[2][y] + y_half_width,
                                       f'{c:.2f}', va='center', ha='center', fontsize=10, color='red')

        self.current_axis.set_ylabel(variable_name + "(Reco) [GeV]", fontsize=30)
        self.current_axis.set_xlabel(variable_name + "(Gen) [GeV]", fontsize=30)

    def show_legend(self, location=(0, 0), **kwargs):
        plt.rcParams['text.usetex'] = True
        self.set_current_axis(location=location)
        self.current_axis.legend(self.legend_handles, self.legend_labels, loc='best', )
        hep.plot.yscale_legend(self.current_axis)
        plt.rcParams['text.usetex'] = False

    def draw_errorbar(self):
        for index, data in enumerate(self.errorbar_list):
            x_data = np.array(data[0])
            y_data = np.array(data[1])

            x_value = x_data.T[0]
            y_value = y_data.T[0]
            x_error = x_data.T[1]
            y_error = y_data.T[1]
            self.set_current_axis(location=self.errorbar_loc[index])
            self.current_axis.errorbar(x_value, y_value,  xerr=x_error, yerr=y_error,
                                       **self.errorbar_kwargs[index])

    def comparison_plot_cosmetics(self, x_variable_name, y_log_scale=False, x_log_scale=False,
                                  bin_width_norm=False, ratio_name='Data/MC'):

        if x_log_scale:
            self.get_axis(location=(0, 0)).set_yscale("log")
        if x_log_scale:
            # FIXME loop over all axis
            self.get_axis(location=(0, 0)).set_xscale("log")
            self.get_axis(location=(1, 0)).set_xscale("log")

        self.get_axis(location=(0, 0)).set_xticklabels([])
        if bin_width_norm:
            self.get_axis(location=(0, 0)).set_ylabel("Events/GeV")
        else:
            self.get_axis(location=(0, 0)).set_ylabel("Events/bin")
        self.show_legend(location=(0, 0))

        self.get_axis(location=(1, 0)).set_ylim(0.4, 1.6)
        self.get_axis(location=(1, 0)).axhline(y=1, linestyle='--', linewidth=1, color='black')
        self.get_axis(location=(1, 0)).set_ylabel(ratio_name)
        self.get_axis(location=(1, 0)).set_xlabel(x_variable_name + " [GeV]")

    def draw_error_boxes(self, default_value, bins, errors,
                         location=(0,0), sys_name='', **kwargs):  # TODO add options for hatch cosmetics

        # default style
        fill = kwargs.get('fill', False)
        edgecolor = kwargs.get('edgecolor', None)
        facecolor = kwargs.get('facecolor', None)
        hatch = kwargs.get('hatch', '///')
        alpha = kwargs.get('alpha', 0.5)

        center = bins[:-1] + np.diff(bins) / 2.
        x_width = np.expand_dims(np.diff(bins) / 2., axis=0)
        x_width = np.append(x_width, x_width, axis=0)

        # symmetric error
        error_boxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum(),
                                 fill=fill, edgecolor=edgecolor, facecolor=facecolor,alpha=alpha,)
                       for x, y, xe, ye in zip(center, default_value, x_width.T, errors.T)]
        pc = PatchCollection(error_boxes, match_original=True, hatch=hatch, linewidth=0.0, zorder=100)
        # for legend only
        legend_ = Rectangle((0, 0), 0.0, 0.0,
                            fill=fill, edgecolor=edgecolor, facecolor=facecolor, label=sys_name, hatch=hatch,
                            alpha=alpha,
                            linewidth=0.5)

        self.set_current_axis(location=location)
        self.current_axis.add_collection(pc)
        self.current_axis.add_patch(legend_)

    def save_fig(self, out_name=''):
        out_file_name = out_name

        # print(f"save plot... {out_file_name}")
        self.fig.savefig(self.base_output_dir + "/" + out_file_name + ".pdf")
        # self.reset()
        plt.close()

