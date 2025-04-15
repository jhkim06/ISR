import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
import numpy as np
from typing import List, Dict, Any
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mplhep.plot import ErrorBarArtists
import math
from Hist import Hist

# PlottableHist class
class PlottableHist(Hist):
    def __init__(self, hist, location=(0, 0), as_denominator=False, as_stack=False, drawn=False, **kwargs):

        super(PlottableHist, self).__init__(hist.raw_root_hist,
                                            hist.label, hist.channel, hist.year, hist.is_measurement)

        self.hist = hist
        self.location = location
        self.as_denominator = as_denominator
        self.as_stack = as_stack
        self.drawn = drawn
        self.plot_kwargs = kwargs

    # def divide(self, other):
    #     return self.hist.divide(other.hist if isinstance(other, PlottableHist) else other)

class Plotter:
    def __init__(self, experiment, base_output_dir, **kwargs):
        self.experiment = experiment
        hep.style.use(self.experiment)
        self.year = ''
        self.lumi = ''
        plt.ioff()
        plt.rcParams['axes.linewidth'] = 2.0
        plt.rcParams['hatch.linewidth'] = 0.5
        plt.rcParams['text.usetex'] = True

        self.rows = 1
        self.cols = 1
        self.fig = None
        self.axs = None
        self.current_axis = None
        self.base_output_dir = base_output_dir
        self.out_dir = ''

        self.plottable_hists: List[PlottableHist] = []

        self.errorbar_list = []
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
        self.plottable_hists.clear()
        self.errorbar_list.clear()
        self.errorbar_loc.clear()
        self.errorbar_kwargs.clear()
        self.legend_handles.clear()
        self.legend_labels.clear()
        self.y_minimum = 0

    def set_experiment_label(self, label="Preliminary", location=(0, 0), **kwargs):
        self.set_current_axis(location=location)
        self.current_axis.draw(self.current_axis.figure.canvas.get_renderer())
        is_data = label != "Simulation"
        label = "" if not is_data else label
        plt.rcParams['text.usetex'] = False
        hep.cms.label(label, data=is_data, fontsize=20, ax=self.current_axis, loc=0, pad=.0, **kwargs)
        plt.rcParams['text.usetex'] = True

    def create_subplots(self, rows=1, cols=1, figsize=(8, 8),
                        left=0.1, right=0.95, bottom=0.1, top=0.9, hspace=0.2, wspace=0.2, **gridspec_kw):
        self.rows = rows
        self.cols = cols
        self.fig, self.axs = plt.subplots(rows, cols, figsize=figsize, gridspec_kw=gridspec_kw)
        plt.tight_layout()
        plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)

    def set_current_axis(self, location=(0, 0)):
        self.current_axis = self.get_axis(location=location)

    def get_axis(self, location=(0, 0)):
        row, col = location
        if self.rows == 1 and self.cols == 1:
            return self.axs
        elif self.rows == 1 or self.cols == 1:
            return self.axs[col] if self.rows == 1 else self.axs[row]
        else:
            return self.axs[row][col]

    def add_hist(self, hist, location=(0, 0), as_denominator=False, as_stack=False, **kwargs):
        p_hist = PlottableHist(hist, location=location, as_stack=as_stack, as_denominator=as_denominator, **kwargs)
        self.plottable_hists.append(p_hist)
        return len(self.plottable_hists) - 1

    def set_ratio_hist(self):
        nominator_index = []
        denominator_index = []
        for i, p_hist in enumerate(self.plottable_hists):
            if p_hist.drawn:
                continue
            if p_hist.as_denominator:
                denominator_index.append(i)
            else:
                nominator_index.append(i)

        if len(denominator_index) == 1 and len(nominator_index) > 1:
            return nominator_index, denominator_index[0]
        elif len(nominator_index) == 1 and len(denominator_index) > 1:
            return nominator_index[0], denominator_index
        elif len(denominator_index) == 1 and len(nominator_index) == 1:
            return nominator_index[0], denominator_index[0]
        else:
            print('Something went wrong in ratio plot setting...')
            exit(1)

    def draw_hist(self):
        bottom = 0
        for p_hist in self.plottable_hists:
            if p_hist.drawn:
                continue

            values, bins, errors = p_hist.to_numpy()
            if p_hist.location == -999:
                continue

            self.set_current_axis(location=p_hist.location)
            # self.y_minimum = min(self.y_minimum, np.min(values)) if self.y_minimum else np.min(values)
            self.y_minimum = np.min(values)

            if not p_hist.as_stack:
                if p_hist.plot_kwargs.get('yerr') is False:
                    errors = False
                    del p_hist.plot_kwargs['yerr']
                xerr = p_hist.plot_kwargs.pop('xerr', False)
                artists = hep.histplot((values, bins), ax=self.current_axis,
                                       yerr=errors, xerr=xerr, **p_hist.plot_kwargs)
                handle = artists[0] if isinstance(artists[0], ErrorBarArtists) or artists[0].legend_artist else (
                    artists[0].stairs)
            else:
                _, _, patches = self.current_axis.hist(bins[:-1], bins, weights=values,
                                                       histtype='bar', bottom=bottom, **p_hist.plot_kwargs)
                bottom += values
                handle = patches[0]

            label = p_hist.plot_kwargs.get("label", None)
            if label:
                self.legend_handles.append(handle)
                self.legend_labels.append(label)

            self.current_axis.set_xlim(bins[0], bins[-1])
            p_hist.drawn = True

    def add_ratio_hists(self, location=(0, 0), **kwargs):
        nominator_index, denominator_index = self.set_ratio_hist()
        use_year_as_ratio_label = all(p.location == -999 for p in self.plottable_hists)

        def get_hist(index):
            if isinstance(index, int):
                return self.plottable_hists[index].hist
            return sum((self.plottable_hists[i].hist for i in index[1:]),
                       self.plottable_hists[index[0]].hist)
            #return sum((self.plottable_hists[i].hist for i in index))

        def is_stackable(indices):
            return all(self.plottable_hists[i].as_stack for i in indices)

        def compute_error_band(hist):
            return hist.divide(hist)

        error_band = None

        if isinstance(nominator_index, int):
            reference_index = nominator_index
            nom_hist = self.plottable_hists[nominator_index].hist
            denom_hist = get_hist(denominator_index)
            ratio_hist = nom_hist.divide(denom_hist)
            error_band = compute_error_band(denom_hist)
        elif isinstance(denominator_index, int):
            denom_hist = self.plottable_hists[denominator_index].hist
            error_band = compute_error_band(denom_hist)
            if is_stackable(nominator_index):
                reference_index = denominator_index
                nom_hist = sum((self.plottable_hists[i].hist for i in nominator_index))
                ratio_hist = nom_hist.divide(denom_hist)
            else:
                for idx in nominator_index:
                    ratio_hist = self.plottable_hists[idx].hist.divide(denom_hist)
                    self.add_ratio_hist(ratio_hist, idx, location=location,
                                        use_year_as_ratio_label=use_year_as_ratio_label, **kwargs)
                self.draw_error_boxes(error_band.to_numpy()[0],
                                      error_band.to_numpy()[1],
                                      error_band.total_sym_sys,
                                      location=location,
                                      **{"facecolor": 'none', "alpha": 0.5, "fill": True, 'hatch': '///'})
                return
        else:
            print("❌ Invalid nominator/denominator config")
            exit(1)

        self.add_ratio_hist(ratio_hist, reference_index, location=location,
                            use_year_as_ratio_label=use_year_as_ratio_label, **kwargs)

        self.draw_error_boxes(error_band.to_numpy()[0],
                              error_band.to_numpy()[1],
                              error_band.total_sym_sys,
                              location=location,
                              **{"facecolor": 'none', "alpha": 0.5, "fill": True, 'hatch': '///'})

    def add_ratio_hist(self, ratio_hist, reference_index, location=(0, 0), use_year_as_ratio_label=False, **kwargs):
        ref_hist = self.plottable_hists[reference_index]
        if use_year_as_ratio_label:
            kwargs['label'] = getattr(ratio_hist, 'year', '')
        if 'color' in ref_hist.plot_kwargs:
            kwargs['color'] = ref_hist.plot_kwargs['color']
        if 'histtype' in ref_hist.plot_kwargs:
            kwargs['histtype'] = ref_hist.plot_kwargs['histtype']

        self.add_hist(ratio_hist, location=location, **kwargs)

    def show_legend(self, location=(0, 0), **kwargs):
        self.set_current_axis(location=location)
        self.current_axis.legend(self.legend_handles, self.legend_labels, loc='best', fontsize=17)
        hep.plot.yscale_legend(self.current_axis)

    def draw_error_boxes(self, default_value, bins, errors,
                         location=(0,0), sys_name='', **kwargs):
        fill = kwargs.get('fill', False)
        edgecolor = kwargs.get('edgecolor', None)
        facecolor = kwargs.get('facecolor', None)
        hatch = kwargs.get('hatch', '///')
        alpha = kwargs.get('alpha', 0.5)

        center = bins[:-1] + np.diff(bins) / 2.
        x_width = np.expand_dims(np.diff(bins) / 2., axis=0)
        x_width = np.append(x_width, x_width, axis=0)

        error_boxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum(),
                                 fill=fill, edgecolor=edgecolor, facecolor=facecolor,alpha=alpha,)
                       for x, y, xe, ye in zip(center, default_value, x_width.T, errors.T)]
        pc = PatchCollection(error_boxes, match_original=True, hatch=hatch, linewidth=0.0, zorder=100)

        legend_ = Rectangle((0, 0), 0.0, 0.0,
                            fill=fill, edgecolor=edgecolor, facecolor=facecolor, label=sys_name, hatch=hatch,
                            alpha=alpha,
                            linewidth=0.5)

        self.set_current_axis(location=location)
        self.current_axis.add_collection(pc)
        self.current_axis.add_patch(legend_)

    def set_common_ratio_plot_cosmetics(self,
                                        x_variable_name,
                                        x_log_scale=False,
                                        y_axis_name='Data',
                                        y_min=0.4,
                                        y_max=1.6
                                        ):

        if x_log_scale:
            self.get_axis(location=(0, 0)).set_xscale("log")
        self.get_axis(location=(0, 0)).set_ylabel(y_axis_name)
        self.show_legend(location=(0, 0))
        # self.adjust_y_scale()

        self.get_axis(location=(0, 0)).axhline(y=1, linestyle='--', linewidth=1, color='black')
        self.get_axis(location=(0, 0)).set_xlabel(x_variable_name + " [GeV]")
        self.get_axis(location=(0, 0)).set_ylim(y_min, y_max)

    def set_common_comparison_plot_cosmetics(self,
                                             x_variable_name,
                                             y_log_scale=False,
                                             x_log_scale=False,
                                             bin_width_norm=False,
                                             ratio_name='Data/MC',
                                             ratio_min = 0.4, ratio_max= 1.6):
        # usual comparison plot cosmetic (2 rows and 1 colum)
        if y_log_scale:
            self.get_axis(location=(0, 0)).set_yscale("log")

        if x_log_scale:
            self.get_axis(location=(0, 0)).set_xscale("log")
            self.get_axis(location=(1, 0)).set_xscale("log")

        # remove default tick labels to avoid clipping
        self.get_axis(location=(0, 0)).set_xticklabels([])

        if bin_width_norm:
            self.get_axis(location=(0, 0)).set_ylabel("Events/GeV")
        else:
            self.get_axis(location=(0, 0)).set_ylabel("Events/bin")

        self.show_legend(location=(0, 0))
        self.adjust_y_scale()

        self.get_axis(location=(1, 0)).set_ylim(ratio_min, ratio_max)
        self.get_axis(location=(1, 0)).axhline(y=1, linestyle='--', linewidth=1, color='black')
        self.get_axis(location=(1, 0)).set_ylabel(ratio_name)
        self.get_axis(location=(1, 0)).set_xlabel(x_variable_name + " [GeV]")

    def adjust_y_scale(self, location=(0, 0)):
        self.set_current_axis(location)
        self.current_axis.set_ylim(ymin=self.y_minimum * 1e-2)

    def add_text(self, text, location=(0, 0), do_magic=True, **kwargs):
        self.set_current_axis(location=location)
        #plt.rcParams['text.usetex'] = True
        at = AnchoredText(text, prop=dict(size=20), **kwargs)
        self.current_axis.add_artist(at)
        if do_magic:
            hep.plot.mpl_magic(self.current_axis)
        #plt.rcParams['text.usetex'] = False

    def add_comparison_pair(self, nominator_hist, denominator_hist,
                            location, ratio_location,
                            nominator_args, denominator_args):
        self.add_hist(nominator_hist, location=location, as_denominator=False, **nominator_args)
        self.add_hist(denominator_hist, location=location, as_denominator=True, **denominator_args)
        self.add_ratio_hists(location=ratio_location)

    def draw_matrix(self, rm_np, variable_name, **kwargs):

        self.set_current_axis((0, 0))
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
                x_half_width = (rm_np[1][x + 1] - rm_np[1][x]) / 2
                y_half_width = (rm_np[2][y + 1] - rm_np[2][y]) / 2
                self.current_axis.text(rm_np[1][x] + x_half_width,
                                       rm_np[2][y] + y_half_width,
                                       f'{c:.2f}', va='center', ha='center', fontsize=15, color='red')

        self.current_axis.set_ylabel(variable_name + "(Reco) [GeV]", fontsize=30)
        self.current_axis.set_xlabel(variable_name + "(Gen) [GeV]", fontsize=30)

    def add_custom_axis_tick_labels(self, locates, labels, location=(0, 0)):
        self.set_current_axis(location=location)
        self.current_axis.set_xticks([])
        self.current_axis.xaxis.set_minor_locator(FixedLocator(locates))
        self.current_axis.xaxis.set_minor_formatter(FixedFormatter(labels))
        self.current_axis.tick_params(axis='x', which='both', rotation=35, labelsize=5)
        for label in self.current_axis.xaxis.get_minorticklabels():
            label.set_horizontalalignment('right')

    def add_errorbar(self, errorbar_data, location=(0, 0), **kwargs):
        self.errorbar_list.append(errorbar_data)
        self.errorbar_loc.append(location)
        self.errorbar_kwargs.append(kwargs)

    # currently used to draw mean pt vs mean mass plot
    def draw_errorbar(self):
        for index, data in enumerate(self.errorbar_list):
            x_data = np.array(data[0])
            y_data = np.array(data[1])

            x_value = x_data.T[0]
            y_value = y_data.T[0]
            x_error = x_data.T[1]
            y_error = y_data.T[1]
            self.set_current_axis(location=self.errorbar_loc[index])
            self.current_axis.errorbar(x_value, y_value, xerr=x_error, yerr=y_error,
                                       **self.errorbar_kwargs[index])

    def draw_vlines(self, vlines, location=(0, 0), **kwargs):
        for line in vlines:
            self.get_axis(location=location).axvline(x=line, linestyle='--', linewidth=1)

    def save_fig(self, out_name=''):
        out_file_name = out_name

        # print(f"save plot... {out_file_name}")
        self.fig.savefig(self.base_output_dir + "/" + out_file_name + ".pdf")
        # self.reset()
        plt.close()

