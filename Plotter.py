import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mplhep.plot import ErrorBarArtists, StairsArtists
import math
from Hist import Hist
import copy
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.container import ErrorbarContainer
import matplotlib.colors as mcolors
from array import array
import ROOT


def extract_color_from_handle(handle):
    """
    Extracts a color from various matplotlib / mplhep artist types
    and returns it as a hex string, e.g. '#578fe0'.
    """

    raw_color = None

    # Case 1: mplhep ErrorBarArtists
    if isinstance(handle, ErrorBarArtists):
        container = handle.errorbar
        if isinstance(container, ErrorbarContainer):
            # this is usually a Line2D
            raw_color = container.lines[0].get_color()

    # Case 2: mplhep StairsArtists
    elif isinstance(handle, StairsArtists):
        # returns an RGBA tuple
        raw_color = handle.stairs.get_edgecolor()

    # Case 3: Line2D (e.g., for errorbar or line plot)
    elif isinstance(handle, mlines.Line2D):
        raw_color = handle.get_color()

    # Case 4: Patch (e.g., from bar hist)
    elif isinstance(handle, mpatches.Patch):
        raw_color = handle.get_facecolor()

    else:
        print(f"Unknown handle type: {type(handle)}")
        return None

    # Convert whatever you got (hex str, name, or RGBA tuple) into pure '#rrggbb'
    try:
        return mcolors.to_hex(raw_color)
    except Exception:
        # fallback: just str it
        return str(raw_color)

def find_non_negative_min(arr):
    # Get non-negative values
    non_negative = arr[arr >= 0]
    # Find the minimum among them
    if non_negative.size > 0:
        minimum = non_negative.min()
        return minimum if minimum > 0 else 1e-2
    else:
        return 1e-2


class PlotItem(Hist):
    def __init__(self, hist, location=(0, 0), as_denominator=False, as_stack=False,
                 drawn=False, not_to_draw=False, use_for_ratio=True, yerr=True, sym_err_name='Total',
                 show_err_band=True,
                 **kwargs):

        if isinstance(hist, tuple) and len(hist) == 3:
            values, bins, errors = hist
            bins = array('f', bins)
            raw_hist = ROOT.TH1D("hist_to_draw", "hist_to_draw", len(bins)-1, bins)
            for index, value in enumerate(values):
                raw_hist.SetBinContent(index+1, value)
                if errors is not None:
                    raw_hist.SetBinError(index+1, errors[index])
                else:
                    raw_hist.SetBinError(index+1, 0)
            hist = Hist(raw_hist, )

        super(PlotItem, self).__init__(hist.raw_root_hist,
                                       label=hist.label,
                                       channel=hist.channel,
                                       year=hist.year,
                                       is_measurement=hist.is_measurement, )

        self.systematic_raw_root_hists = copy.deepcopy(hist.systematic_raw_root_hists)
        # self.systematics = copy.deepcopy(hist.systematics)
        self.compute_systematic_rss_per_sysname()

        self.hist = hist
        self.location = location
        self.as_denominator = as_denominator
        self.as_stack = as_stack
        self.drawn = drawn
        self.not_to_draw = not_to_draw
        self.use_for_ratio = use_for_ratio
        self.plot_kwargs = kwargs
        self.show_y_err = yerr
        self.sym_err_name = sym_err_name
        self.show_err_band = show_err_band

    def divide(self, other=None):
        divided_hist = super().divide(other, use_default_denominator=True)
        return PlotItem(divided_hist, location=self.location, as_stack=self.as_stack, as_denominator=self.as_denominator,
                        use_for_ratio=self.use_for_ratio, not_to_draw=self.not_to_draw,
                        yerr=self.show_y_err,
                        sym_err_name=self.sym_err_name, show_err_band=self.show_err_band,
                        **self.plot_kwargs)

    def __add__(self, other=None, c1=1):
        added_hist = super().__add__(other, c1)
        return PlotItem(added_hist, location=self.location, as_stack=self.as_stack, as_denominator=self.as_denominator,
                        use_for_ratio=self.use_for_ratio, not_to_draw=self.not_to_draw,
                        yerr=self.show_y_err,
                        sym_err_name=self.sym_err_name, show_err_band=self.show_err_band,
                        **self.plot_kwargs)


# Plotter get Hist and draw it
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
        self._stack_bottoms: Dict[Tuple[int, int], np.ndarray] = {}
        self.current_axis = None
        self.base_output_dir = base_output_dir
        self.out_dir = ''

        self.plot_items: List[PlotItem] = []

        self.errorbar_list = []
        self.errorbar_loc = []
        self.errorbar_kwargs: List[Dict[str, Any]] = []

        self.legend_handles: Dict[Tuple[int, int], list] = {}
        self.legend_labels: Dict[Tuple[int, int], list] = {}

        self.y_minimum = 999.

    def reset(self):
        self.rows = 1
        self.cols = 1
        del self.fig
        del self.axs
        self._stack_bottoms.clear()
        self.current_axis = None
        self.plot_items.clear()
        self.errorbar_list.clear()
        self.errorbar_loc.clear()
        self.errorbar_kwargs.clear()
        self.legend_handles.clear()
        self.legend_labels.clear()
        self.y_minimum = 999.

    def init_plotter(self, figsize=(8,8), rows=1, cols=1,
                     left=0.15, right=0.95, bottom=0.15, top=0.95,):
        if rows == 2 and cols == 1:
            self.create_subplots(rows, cols, figsize=figsize,
                                 left=left, right=right, hspace=0.0, bottom=bottom, height_ratios=[1, 0.3])
        elif rows == 1 and cols == 1:
            self.create_subplots(rows, cols, figsize=figsize,
                                 left=left, right=right, hspace=0.0, bottom=bottom)
        else:
            # FIXME
            self.create_subplots(rows, cols, figsize=figsize,)

        self._stack_bottoms.clear()
        for r in range(rows):
            for c in range(cols):
                self._stack_bottoms[(r,c)] = 0
                self.legend_handles[(r,c)] = []
                self.legend_labels[(r,c)] = []

    def create_subplots(self, rows=1, cols=1, figsize=(8, 8),
                        left=0.1, right=0.95, bottom=0.1, top=0.9, hspace=0.2, wspace=0.2, **gridspec_kw):
        self.rows = rows
        self.cols = cols
        self.fig, self.axs = plt.subplots(rows, cols, figsize=figsize, gridspec_kw=gridspec_kw)
        plt.tight_layout()
        plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)

    def save_and_reset_plotter(self, hist_name, postfix=''):
        self.save_fig(hist_name + postfix)
        self.reset()

    def set_experiment_label(self, label="Preliminary", location=(0, 0), **kwargs):
        self.set_current_axis(location=location)
        self.current_axis.draw(self.current_axis.figure.canvas.get_renderer())
        is_data = label != "Simulation"
        label = "" if not is_data else label
        plt.rcParams['text.usetex'] = False
        hep.cms.label(label, data=is_data, fontsize=20, ax=self.current_axis, loc=0, pad=.0, **kwargs)
        plt.rcParams['text.usetex'] = True

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

    # hist can be Hist or (values, bin, error)
    def add_hist(self, hist, location=(0, 0), as_denominator=False, as_stack=False,
                 not_to_draw=False,
                 use_for_ratio=True, yerr=True, sym_err_name='Total', show_err_band=True,
                 **kwargs):
        p_hist = PlotItem(hist, location=location, as_stack=as_stack, as_denominator=as_denominator,
                          use_for_ratio=use_for_ratio, not_to_draw=not_to_draw,
                          yerr=yerr,
                          sym_err_name=sym_err_name, show_err_band=show_err_band,
                          **kwargs)
        self.plot_items.append(p_hist)
        return len(self.plot_items) - 1

    def set_ratio_hist(self):
        nominator_index = []
        denominator_index = []
        for i, p_hist in enumerate(self.plot_items):
            if not p_hist.use_for_ratio:
                continue
            if p_hist.as_denominator:
                denominator_index.append(i)
            else:
                nominator_index.append(i)
            p_hist.use_for_ratio = False  # do not use this histogram for ratio plot anymore

        if len(denominator_index) == 1 and len(nominator_index) > 1:
            return nominator_index, denominator_index[0]
        elif len(nominator_index) == 1 and len(denominator_index) > 1:
            return nominator_index[0], denominator_index
        elif len(denominator_index) == 1 and len(nominator_index) == 1:
            return nominator_index[0], denominator_index[0]
        else:
            print('Something went wrong in ratio plot setting...')
            print(f'# denominator {len(denominator_index)}, # nominator {len(nominator_index)}')
            exit(1)

    def draw_hist(self):
        for item_index, p_hist in enumerate(self.plot_items):
            if p_hist.drawn:
                continue
            # TODO allow case without TH1
            values, bins, errors = p_hist.to_numpy()
            if p_hist.not_to_draw:  # FIXME
                continue

            self.set_current_axis(location=p_hist.location)
            # self.y_minimum = min(self.y_minimum, find_non_negative_min(values)) if self.y_minimum else find_non_negative_min(values)
            self.y_minimum = find_non_negative_min(values)

            if p_hist.as_stack:
                handle = self._draw_stack(p_hist, values, bins)
                this_color = handle.get_facecolor()
            else:
                handle = self._draw_single(p_hist, values, bins, errors)
                this_color = extract_color_from_handle(handle)
                if 'color' not in p_hist.plot_kwargs and this_color is not None:
                    p_hist.plot_kwargs['color'] = this_color
                if p_hist.show_err_band:
                    self.draw_error_box_for_plot_item(item_index, p_hist.location)
            # check if kwargs has color, if not, draw and then get the color and set the color in the kwargs

            label = p_hist.plot_kwargs.get("label", None)
            if label:
                self.legend_handles[p_hist.location].append(handle)
                self.legend_labels[p_hist.location].append(label)

            self.current_axis.set_xlim(bins[0], bins[-1])
            p_hist.drawn = True

    def _draw_stack(self, item, values, bins):
        bottom = self._stack_bottoms[item.location]
        _, _, patches = self.current_axis.hist(bins[:-1], bins, weights=values,
                                histtype='bar', bottom=bottom, **item.plot_kwargs)
        self._stack_bottoms[item.location] += values
        return patches[0]

    def _draw_single(self, item,
                     values, bins, errs):
        plot_kwargs = copy.deepcopy(item.plot_kwargs)
        # plot_kwargs |= {"yerr": False}
        if item.show_y_err:
            yerr = errs
        else:
            yerr = False

        artist = hep.histplot((values, bins), ax=self.current_axis, yerr=yerr, xerr=False,
                              **plot_kwargs)[0]
        # belows are required to show legend properly
        if isinstance(artist, ErrorBarArtists) or artist.legend_artist:
            return artist
        if isinstance(artist, StairsArtists):
            return artist.stairs
        return artist

    def get_hist(self, index):
        if isinstance(index, int):
            return self.plot_items[index].hist
        # FIXME type is Hist not PlotItem
        return sum((self.plot_items[i] for i in index[1:]),
                   self.plot_items[index[0]])

    def draw_ratio_hists(self, location=(0, 0), show_normalized_error_band=True):
        nominator_index, denominator_index = self.set_ratio_hist()

        def is_stackable(indices):
            return all(self.plot_items[i].as_stack for i in indices)

        # Case 1: Both nominator and denominator are index
        if isinstance(nominator_index, int) and isinstance(denominator_index, int):
            nom_hist = self.plot_items[nominator_index]
            denom_hist = self.plot_items[denominator_index]
            ratio_hist = nom_hist.divide(denom_hist)  # TODO define divide for PlotItem!
            self.add_ratio_hist(ratio_hist, location)

        # Case 2: Only nominator is index (denominator is list)
        elif isinstance(nominator_index, int):
            nom_hist = self.plot_items[nominator_index]
            denom_hist = self.get_hist(denominator_index)
            ratio_hist = nom_hist.divide(denom_hist)
            self.add_ratio_hist(ratio_hist, location)

        # Case 3: Only denominator is index (nominator is list)
        elif isinstance(denominator_index, int):
            denom_hist = self.plot_items[denominator_index]

            if is_stackable(nominator_index):
                nom_hist = self.get_hist(nominator_index)
                ratio_hist = nom_hist.divide(denom_hist)
                self.add_ratio_hist(ratio_hist, location)
            else:
                for idx in nominator_index:
                    ratio_hist = self.plot_items[idx].divide(denom_hist)
                    self.add_ratio_hist(ratio_hist, location)
                    self.draw_hist()
                self.draw_normalized_error_band(denom_hist, location=location)
                return
        else:
            print("Invalid nominator/denominator config")
            exit(1)

        if show_normalized_error_band:
            self.draw_normalized_error_band(denom_hist, location=location)
        self.draw_hist()
        # TODO add on/off option

    def draw_normalized_error_band(self, hist, location=(0, 0)):
        error_band = hist.create_normalized_error_band()
        if 'histtype' in hist.plot_kwargs and hist.plot_kwargs['histtype'] == 'errorbar':
            self.add_hist(error_band, location=location, use_for_ratio=False, yerr=False, show_err_band=False,
                      **hist.plot_kwargs)
        self.draw_error_boxes(error_band.to_numpy()[0],
                              error_band.to_numpy()[1],
                              error_band.total_sym_err_array,
                              location=location,
                              **{"facecolor": 'none', "alpha": 0.5, "fill": True, 'hatch': '///'})

    def draw_error_box_for_plot_item(self, index, location=(0, 0),):
        plot_item = self.plot_items[index]

        self.draw_error_boxes(plot_item.to_numpy()[0],
                              plot_item.to_numpy()[1],
                              plot_item.get_sym_sys_err_array(plot_item.sym_err_name),
                              location=location,
                              **{"facecolor": plot_item.plot_kwargs['color'],
                                 "alpha": 0.2, "fill": True, 'hatch': None})

    def add_ratio_hist(self, ratio_hist, location=(0, 0)):
        ratio_hist.use_for_ratio = False
        ratio_hist.show_y_err = False
        ratio_hist.location = location
        ratio_hist.not_to_draw = False
        self.plot_items.append(ratio_hist)
        # return self.add_hist(ratio_hist, location=location, use_for_ratio=False, yerr=False, **kwargs)

    def show_legend(self, location=(0, 0), **kwargs_):
        self.set_current_axis(location=location)
        kwargs = {"loc": 'best', 'fontsize': 17} | kwargs_
        #print(self.legend_handles[location])
        #print(self.legend_labels[location])
        self.current_axis.legend(self.legend_handles[location],
                                 self.legend_labels[location], **kwargs)
        try:
            hep.plot.yscale_legend(self.current_axis)
        except RuntimeError:
            pass

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

    def set_isr_plot_cosmetics(self, channel, ):

        plt.grid(True, which='both', axis='x', linestyle='--', linewidth=0.7)
        plt.grid(True, which='major', axis='y', linestyle='--', linewidth=0.7)

        self.get_axis(location=(0, 0)).set_xlim(54, 300)
        self.get_axis(location=(0, 0)).set_ylim(15, 29)
        self.get_axis(location=(0, 0)).set_xscale("log")

        self.get_axis(location=(0, 0)).set_ylabel(r"Mean $p_{T}^{"+channel+"}$ GeV")
        self.get_axis(location=(0, 0)).set_xlabel(r"Mean $m^{"+channel+"}$ GeV")

    def set_common_ratio_plot_cosmetics(self,
                                        x_variable_name,
                                        x_log_scale=False,
                                        y_axis_name='Data',
                                        y_min=0.4,
                                        y_max=1.6):

        if x_log_scale:
            self.get_axis(location=(0, 0)).set_xscale("log")
        self.get_axis(location=(0, 0)).set_ylabel(y_axis_name)
        self.show_legend(location=(0, 0),)
        # self.adjust_y_scale()

        self.get_axis(location=(0, 0)).axhline(y=1, linestyle='--', linewidth=1, color='black')
        self.get_axis(location=(0, 0)).set_xlabel(x_variable_name)
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
        self.get_axis(location=(1, 0)).set_xlabel(x_variable_name)

    def adjust_y_scale(self, location=(0, 0)):
        self.set_current_axis(location)
        self.current_axis.set_ylim(ymin=self.y_minimum * 1e-1)

    def add_text(self, text, location=(0, 0), do_magic=True, **kwargs):
        self.set_current_axis(location=location)
        #plt.rcParams['text.usetex'] = True
        color='black'
        if 'color' in kwargs:
            color = kwargs['color']
            del kwargs['color']
        at = AnchoredText(text, prop=dict(size=20, color=color), **kwargs)
        self.current_axis.add_artist(at)
        if do_magic:
            hep.plot.mpl_magic(self.current_axis)  # error without legend?
        #plt.rcParams['text.usetex'] = False

    def add_comparison_pair(self, nominator_hist, denominator_hist,
                            location, ratio_location,
                            nominator_args, denominator_args):
        self.add_hist(nominator_hist, location=location, as_denominator=False, **nominator_args)
        self.add_hist(denominator_hist, location=location, as_denominator=True, **denominator_args)
        self.draw_ratio_hists(location=ratio_location)

    def draw_matrix(self, rm_np, x_axis_label="", y_axis_label="", **kwargs):
        self.set_current_axis((0, 0))
        # hep.hist2dplot(rm_np, norm=mcolors.LogNorm(), ax=self.current_axis, **kwargs)
        hep.hist2dplot(rm_np, ax=self.current_axis, **kwargs)

        if (rm_np[1][-1] - rm_np[1][0] >= 5e2 and
                rm_np[2][-1] - rm_np[2][0] >= 5e2):
            self.current_axis.set_xscale("log")
            self.current_axis.set_yscale("log")

            # self.current_axis.set_xlim(0.2, rm_np[1][-1])
            # self.current_axis.set_ylim(0.2, rm_np[2][-1])

        show_number = False
        if show_number:
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

        self.current_axis.set_ylabel(y_axis_label, fontsize=30)
        self.current_axis.set_xlabel(x_axis_label, fontsize=30)

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

            # use the last bin (total error) as error bar
            x_error = x_data.T[-1]
            y_error = y_data.T[-1]

            self.set_current_axis(location=self.errorbar_loc[index])
            #print(len(x_value), len(y_value))
            self.current_axis.errorbar(x_value, y_value, xerr=x_error, yerr=y_error,
                                       **self.errorbar_kwargs[index])

            handles, labels = self.current_axis.get_legend_handles_labels()
            label = self.errorbar_kwargs[index].get("label", None)
            if label:

                self.legend_handles[self.errorbar_loc[index]].append(handles[-1])
                self.legend_labels[self.errorbar_loc[index]].append(labels[-1])

    def draw_vlines(self, vlines, location=(0, 0), **kwargs):
        for line in vlines:
            self.get_axis(location=location).axvline(x=line, **kwargs)

    def save_fig(self, out_name=''):
        out_file_name = out_name

        # print(f"save plot... {out_file_name}")
        self.fig.savefig(self.base_output_dir + "/" + out_file_name + ".pdf")
        # self.reset()
        plt.close()

