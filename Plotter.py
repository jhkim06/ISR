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
import gc


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
        #print("use this color!")
        #print(mcolors.to_hex(raw_color))
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
                 drawn=False, not_to_draw=False, use_for_ratio=True, yerr=True, show_x_err=False,
                 sym_err_name='Total',
                 show_err_band=True, err_band_alpha=0.8, err_band_fill=True, err_band_hatch=None,
                 err_band_sys_name='',
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
        self.show_x_err = show_x_err
        self.sym_err_name = sym_err_name
        self.show_err_band = show_err_band
        self.err_band_alpha = err_band_alpha
        self.err_band_fill = err_band_fill
        self.err_band_hatch = err_band_hatch
        self.err_band_sys_name = err_band_sys_name

    def divide(self, other=None):
        divided_hist = super().divide(other, use_default_denominator=True)
        return PlotItem(divided_hist, location=self.location, as_stack=self.as_stack, as_denominator=self.as_denominator,
                        use_for_ratio=self.use_for_ratio, not_to_draw=self.not_to_draw,
                        yerr=self.show_y_err, show_x_err=self.show_x_err,
                        sym_err_name=self.sym_err_name, show_err_band=self.show_err_band,
                        err_band_alpha=self.err_band_alpha, err_band_fill=self.err_band_fill,
                        err_band_hatch=self.err_band_hatch, err_band_sys_name=self.err_band_sys_name,
                        **self.plot_kwargs)

    def __add__(self, other=None, c1=1):
        added_hist = super().__add__(other, c1)
        return PlotItem(added_hist, location=self.location, as_stack=self.as_stack, as_denominator=self.as_denominator,
                        use_for_ratio=self.use_for_ratio, not_to_draw=self.not_to_draw,
                        yerr=self.show_y_err, show_x_err=self.show_x_err,
                        sym_err_name=self.sym_err_name, show_err_band=self.show_err_band,
                        err_band_alpha=self.err_band_alpha, err_band_fill=self.err_band_fill,
                        err_band_hatch=self.err_band_hatch, err_band_sys_name=self.err_band_sys_name,
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

        self.show_legend_for_normalized_error=False

        self.y_minimum = 999.

    def reset(self):
        self.rows = 1
        self.cols = 1
        self.fig.clf()
        #plt.close("all")
        #plt.close(self.fig)
        del self.fig
        self.fig = None
        self.axs = None
        self._stack_bottoms.clear()
        self.current_axis = None
        self.plot_items.clear()
        self.errorbar_list.clear()
        self.errorbar_loc.clear()
        self.errorbar_kwargs.clear()
        self.legend_handles.clear()
        self.legend_labels.clear()
        self.y_minimum = 999.
        #print(len(plt.get_fignums()), "figures still open")
        gc.collect()

    def init_plotter(self, figsize=(8,8), rows=1, cols=1,
                     left=0.15, right=0.95, bottom=0.15, top=0.95,):
        if rows == 3 and cols == 1:
            self.create_subplots(rows, cols, figsize=figsize,
                                 left=left, right=right, hspace=0.0, bottom=bottom, height_ratios=[1, 0.32, 0.27])
        elif rows == 2 and cols == 1:
            self.create_subplots(rows, cols, figsize=figsize,
                                 left=left, right=right, hspace=0.0, bottom=bottom, height_ratios=[1, 0.3])
        elif rows == 1 and cols == 1:
            self.create_subplots(rows, cols, figsize=figsize,
                                 left=left, right=right, hspace=0.0, bottom=bottom)
        else:
            # FIXME
            self.create_subplots(rows, cols, figsize=figsize,)

        self.set_current_axis()

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
        #self.fig, self.axs = plt.subplots(rows, cols, figsize=figsize, gridspec_kw=gridspec_kw)

        # https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
        self.fig = plt.figure(num=1, figsize=figsize)
        self.fig.set_size_inches(*figsize, forward=True)

        # Step 2: Build a GridSpec with all the kw
        gs_kw = dict(hspace=hspace, wspace=wspace, **gridspec_kw)
        gs = self.fig.add_gridspec(rows, cols, **gs_kw)

        # Step 3: Add each subplot cell-by-cell
        axes = []
        for irow in range(rows):
            for icol in range(cols):
                ax = self.fig.add_subplot(gs[irow, icol])
                axes.append(ax)

        # reshape to same layout as subplots() would have
        self.axs = np.array(axes).reshape(rows, cols) \
                   if (rows>1 and cols>1) \
                   else (axes[0] if (rows==1 and cols==1) else axes)

        #plt.tight_layout()
        plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top)

    def save_and_reset_plotter(self, hist_name, out_sub_dirs=''):
        self.save_fig(hist_name, out_sub_dirs)
        self.reset()

    def get_n_colors(self, n):
        # Get the original tab10 colormap
        tab10 = plt.get_cmap('tab10')
        colors = [tab10(i) for i in range(10)]

        # Generate extended colors by modifying brightness/saturation
        extended_colors = []
        for i in range(n):
            base_color = colors[i % 10]
            # Slightly adjust the brightness and saturation
            factor = 0.9 + 0.1 * (i // 10)  # Change factor for every new cycle
            new_color = tuple(min(1, c * factor) for c in base_color)
            extended_colors.append(new_color)

        return extended_colors

    # Recommended to use, after draw all component
    def set_experiment_label(self, label="Preliminary", location=(0, 0),
                             **kwargs):
        self.set_current_axis(location=location)
        is_data = label != "Simulation"
        label = "" if not is_data else label

        plt.rcParams['text.usetex'] = False
        hep.cms.label(label, data=is_data, fontsize=20, ax=self.current_axis, loc=0, pad=.0,
                      **kwargs)
        plt.rcParams['text.usetex'] = True

        self.fig.canvas.draw()
        ax = self.current_axis
        offset_text = ax.yaxis.get_offset_text()  # scientific exponent
        exponent_str = offset_text.get_text()

        if exponent_str and ax.get_yscale() != "log":
            offset_text.set_fontsize(14)
            self.fig.canvas.draw()
        #    renderer = self.fig.canvas.get_renderer()
        #    bbox = offset_text.get_window_extent(renderer)
        #    x1_disp, y1_disp = bbox.x1, bbox.y1
        #    x1_axes, y1_axes = ax.transAxes.inverted().transform((x1_disp, y1_disp))
        #    all_texts = ax.texts
        #    x, _ = all_texts[1].get_position()
        #    all_texts[1].set_x(x1_axes+x)  # 'CMS'
        #    x, _ = all_texts[2].get_position()  # 'Preliminary' for example
        #    all_texts[2].set_x(x1_axes+x)
        #else:
        #    print("No offset text present")

        # move the patch
        self.fig.canvas.draw()

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
                 use_for_ratio=True, yerr=True, show_x_err=False, sym_err_name='Total', show_err_band=True,
                 err_band_alpha=1.0, err_band_fill=True, err_band_hatch=None, err_band_sys_name="",
                 **kwargs):

        p_hist = PlotItem(hist, location=location, as_stack=as_stack, as_denominator=as_denominator,
                          use_for_ratio=use_for_ratio, not_to_draw=not_to_draw,
                          yerr=yerr, show_x_err=show_x_err,
                          sym_err_name=sym_err_name, show_err_band=show_err_band,
                          err_band_alpha=err_band_alpha, err_band_fill=err_band_fill,
                          err_band_hatch=err_band_hatch, err_band_sys_name=err_band_sys_name,
                          **kwargs)

        self.plot_items.append(p_hist)

        return len(self.plot_items)-1

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
            values, bins, errors = p_hist.to_numpy(stat=True)
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

        artist = hep.histplot((values, bins), ax=self.current_axis, yerr=yerr, xerr=item.show_x_err,
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

        return sum((self.plot_items[i] for i in index[-2::-1]),  # FIXME check length of items
                   self.plot_items[index[-1]])

    def draw_ratio_hists(self, location=(0, 0), show_y_error=True, show_error_band=True,
                         show_normalized_error_band=True, sym_err_name='Total',
                         show_legend_for_normalized_error=False):
        nominator_index, denominator_index = self.set_ratio_hist()
        self.show_legend_for_normalized_error=show_legend_for_normalized_error

        def is_stackable(indices):
            return all(self.plot_items[i].as_stack for i in indices)

        # Case 1: Both nominator and denominator are index
        if isinstance(nominator_index, int) and isinstance(denominator_index, int):
            nom_hist = self.plot_items[nominator_index]
            denom_hist = self.plot_items[denominator_index]
            ratio_hist = nom_hist.divide(denom_hist)
            self.add_ratio_hist(ratio_hist, location, show_y_error, show_error_band, sym_err_name)

        # Case 2: Only nominator is index (denominator is list)
        elif isinstance(nominator_index, int):
            nom_hist = self.plot_items[nominator_index]
            denom_hist = self.get_hist(denominator_index)
            ratio_hist = nom_hist.divide(denom_hist)
            self.add_ratio_hist(ratio_hist, location, show_y_error, show_error_band, sym_err_name)

        # Case 3: Only denominator is index (nominator is list)
        elif isinstance(denominator_index, int):
            denom_hist = self.plot_items[denominator_index]

            if is_stackable(nominator_index):
                nom_hist = self.get_hist(nominator_index)
                ratio_hist = nom_hist.divide(denom_hist)
                self.add_ratio_hist(ratio_hist, location, show_y_error, show_error_band, sym_err_name)
            else:
                for idx in nominator_index:
                    ratio_hist = self.plot_items[idx].divide(denom_hist)
                    self.add_ratio_hist(ratio_hist, location, show_y_error, show_error_band, sym_err_name)
                    #print(ratio_hist.plot_kwargs)
                    self.draw_hist()
                if show_normalized_error_band:
                    self.draw_normalized_error_band(denom_hist, location=location)
                self.draw_hist()
                return
        else:
            print("Invalid nominator/denominator config")
            exit(1)

        if show_normalized_error_band:
            self.draw_normalized_error_band(denom_hist, location=location)
        self.draw_hist()

    def draw_normalized_error_band(self, hist, location=(0, 0)):
        error_band = hist.create_normalized_error_band()
        if ('histtype' in hist.plot_kwargs and hist.plot_kwargs['histtype'] == 'errorbar' and
                hist.err_band_hatch is not None):
            kwargs=dict()
            if hist.plot_kwargs:
                kwargs = kwargs | hist.plot_kwargs
            kwargs["markersize"] = 3
            self.add_hist(error_band, location=location,
                          use_for_ratio=False, yerr=True, show_err_band=False,
                          **kwargs)
        sys_name = hist.err_band_sys_name
        sys_stat_name = 'Stat.'
        if not self.show_legend_for_normalized_error:
            sys_name = ''
            sys_stat_name = ''
        self.draw_error_boxes(error_band.to_numpy()[0],
                              error_band.to_numpy()[1],
                              error_band.total_sym_err_array,
                              location=location,
                              zorder=0,
                              sys_name=sys_name,
                              **{"facecolor": hist.plot_kwargs['color'],
                                 "alpha": 0.8,
                                 "fill": hist.err_band_fill,
                                 'hatch': hist.err_band_hatch})

        # stat
        self.draw_error_boxes(error_band.to_numpy()[0],
                              error_band.to_numpy()[1],
                              error_band.get_sym_sys_err_array('stat'),
                              location=location,
                              zorder=1,
                              sys_name=sys_stat_name,
                              **{"facecolor": 'mistyrose',  # TODO
                                 "alpha": 0.8,
                                 "fill": hist.err_band_fill,
                                 'hatch': hist.err_band_hatch})

    def draw_error_box_for_plot_item(self, index, location=(0, 0),):
        plot_item = self.plot_items[index]

        sys_name = plot_item.err_band_sys_name
        sys_stat_name = 'Stat.'
        if self.show_legend_for_normalized_error:
            sys_name = ''
            sys_stat_name = ''
        if plot_item.err_band_alpha != 0.0:
            zorder = 1
            if plot_item.is_measurement:
                zorder = 1001
            self.draw_error_boxes(plot_item.to_numpy()[0],
                                  plot_item.to_numpy()[1],
                                  plot_item.get_sym_sys_err_array(plot_item.sym_err_name),
                                  location=location,
                                  zorder=zorder,
                                  sys_name=sys_name,
                                  **{"facecolor": plot_item.plot_kwargs['color'],
                                     "alpha": plot_item.err_band_alpha,
                                     "fill": plot_item.err_band_fill,
                                     'hatch': plot_item.err_band_hatch})

            if plot_item.err_band_hatch is None:
                zorder = 1
                if plot_item.is_measurement:
                    zorder = 1002
                self.draw_error_boxes(plot_item.to_numpy()[0],
                                      plot_item.to_numpy()[1],
                                      plot_item.get_sym_sys_err_array('stat'),
                                      location=location,
                                      zorder=zorder,
                                      sys_name=sys_stat_name,
                                      **{"facecolor": 'mistyrose',  # TODO
                                         "alpha": 0.8,
                                         "fill": plot_item.err_band_fill,
                                         'hatch': plot_item.err_band_hatch})

    def add_ratio_hist(self, ratio_hist, location=(0, 0), show_y_err=True,
                       show_err_band=True, sys_err_name='Total'):
        ratio_hist.use_for_ratio = False

        if ratio_hist.err_band_hatch is not None:
            if show_y_err:
                ratio_hist.show_y_err = True
        else:
            # if hatch is not used
            ratio_hist.show_y_err = False
            ratio_hist.show_x_err = False
        #ratio_hist.show_y_err = show_y_err
        ratio_hist.location = location
        ratio_hist.not_to_draw = False
        ratio_hist.show_err_band = show_err_band
        ratio_hist.sym_err_name = sys_err_name
        if ratio_hist.err_band_fill == True and ratio_hist.err_band_alpha==0.0:  # assume it is suppressed temporarily
            ratio_hist.err_band_alpha = 0.8
        self.plot_items.append(ratio_hist)
        # return self.add_hist(ratio_hist, location=location, use_for_ratio=False, yerr=False, **kwargs)

    def show_legend(self, location=(0, 0), reverse=True, data_label_first=False,
                    show_only_sys_legends=False, font_size=17,
                    **kwargs_):
        self.set_current_axis(location=location)
        kwargs = {"loc": 'best', 'fontsize': font_size} | kwargs_
        if reverse:
            self.current_axis.legend(self.legend_handles[location][::-1],
                                     self.legend_labels[location][::-1], **kwargs)
        else:
            if data_label_first:
                data_index = -1
                for index, label in enumerate(self.legend_labels[location]):
                    if "Data" in label:
                        data_index = index
                        break
                if data_index != -1:
                    labels = copy.deepcopy(self.legend_labels[location])
                    handles = copy.deepcopy(self.legend_handles[location])

                    item_to_move = labels.pop(data_index)
                    labels.insert(0, item_to_move)

                    item_to_move = handles.pop(data_index)
                    handles.insert(0, item_to_move)

                    self.current_axis.legend(handles, labels, **kwargs)
            else:
                if show_only_sys_legends:
                    handles = []
                    labels = []
                    for index, label in enumerate(self.legend_labels[location]):
                        if ("Stat." in label or "Total uncertainty" in label or
                                r"Syst.(Theory $\oplus$ measurement)"  in label or
                                r"Scale $\oplus$ PDF $\oplus$ $\alpha_s$" in label):
                            handles.append(self.legend_handles[location][index])
                            labels.append(label)
                    if len(handles) > 0:
                        n_col = len(handles)
                        kwargs = kwargs | {"ncol": n_col}
                        self.current_axis.legend(handles, labels, **kwargs)
                else:
                    self.current_axis.legend(self.legend_handles[location],
                                             self.legend_labels[location], **kwargs)
        try:
            hep.plot.yscale_legend(self.current_axis)
        except (RuntimeError, StopIteration):
            pass

    def draw_error_boxes(self, default_value, bins, errors,
                         location=(0,0), sys_name='', zorder=10, **kwargs):
        fill = kwargs.get('fill', False)
        edgecolor = kwargs.get('edgecolor', None)
        facecolor = kwargs.get('facecolor', None)
        hatch = kwargs.get('hatch', '///')
        alpha = kwargs.get('alpha', 0.8)

        center = bins[:-1] + np.diff(bins) / 2.
        x_width = np.expand_dims(np.diff(bins) / 2., axis=0)
        x_width = np.append(x_width, x_width, axis=0)

        error_boxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum(),
                                 fill=fill, edgecolor=edgecolor, facecolor=facecolor,alpha=alpha,)
                       for x, y, xe, ye in zip(center, default_value, x_width.T, errors.T)]
        pc = PatchCollection(error_boxes, match_original=True, hatch=hatch, linewidth=0.0, zorder=zorder)

        legend_ = Rectangle((0, 0), 0.0, 0.0,
                            fill=fill, edgecolor=edgecolor, facecolor=facecolor, label=sys_name, hatch=hatch,
                            alpha=alpha,
                            linewidth=0.5)

        self.set_current_axis(location=location)
        # FIXME
        if sys_name:
            self.legend_handles[location].append(legend_)
            self.legend_labels[location].append(sys_name)

        self.current_axis.add_collection(pc)
        self.current_axis.add_patch(legend_)

    def set_isr_plot_cosmetics(self, channel, y_min=15, y_max=29):

        plt.grid(True, which='both', axis='x', linestyle='--', linewidth=0.7)
        plt.grid(True, which='major', axis='y', linestyle='--', linewidth=0.7)

        self.get_axis(location=(0, 0)).set_xlim(54, 300)
        self.get_axis(location=(0, 0)).set_ylim(y_min, y_max)
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
                                             font_size=17,
                                             y_min_scale=1.0,
                                             ratio_min = 0.6, ratio_max= 1.4,
                                             show_sys_legend_for_ratio=False,
                                             rows=2):
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

        self.show_legend(location=(0, 0), font_size=font_size)
        if show_sys_legend_for_ratio:
            self.show_legend(location=(1, 0), font_size=10, reverse=False,
                             show_only_sys_legends=show_sys_legend_for_ratio)
        if y_min_scale != 1.0:
            self.adjust_y_scale(scale=y_min_scale)

        self.get_axis(location=(1, 0)).set_ylim(ratio_min, ratio_max)
        self.get_axis(location=(1, 0)).axhline(y=1, linestyle='--', linewidth=1, color='black', zorder=1000)
        self.get_axis(location=(1, 0)).set_ylabel(ratio_name)
        if ratio_min == 0.6 and ratio_max == 1.4:
            self.get_axis(location=(1, 0)).set_yticks([0.8, 1.0, 1.2])
        elif ratio_min == 0.8 and ratio_max == 1.2:
            self.get_axis(location=(1, 0)).set_yticks([0.9, 1.0, 1.1])
        else:
            pass
        self.get_axis(location=(1, 0)).grid(True, which='major', axis='y', linestyle='-', linewidth=0.7)

        if rows > 2:
            self.get_axis(location=(1, 0)).set_xticklabels([])
            self.get_axis(location=(2, 0)).set_ylim(ratio_min, ratio_max)
            self.get_axis(location=(2, 0)).axhline(y=1, linestyle='--', linewidth=1, color='black', zorder=1000)
            if ratio_min == 0.6 and ratio_max == 1.4:
                self.get_axis(location=(2, 0)).set_yticks([0.8, 1.0, 1.2])
            else:
                self.get_axis(location=(2, 0)).set_yticks([0.9, 1.0, 1.1])
            self.get_axis(location=(2, 0)).set_xlabel(x_variable_name)
            self.get_axis(location=(2, 0)).grid(True, which='major', axis='y', linestyle='-', linewidth=0.7)
        else:
            self.get_axis(location=(1, 0)).set_xlabel(x_variable_name)

    def adjust_y_scale(self, location=(0, 0), scale=0.1):
        self.set_current_axis(location)
        self.current_axis.set_ylim(ymin=self.y_minimum * scale)

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
            #hep.plot.mpl_magic(self.current_axis)  # error without legend?
            try:
                hep.plot.mpl_magic(self.current_axis)
            except RuntimeError as e:
                # swallow only the “Could not fit legend” error,
                # let any other RuntimeError bubble up
                if "Could not fit legend" not in str(e):
                    raise
        #plt.rcParams['text.usetex'] = False

    def add_comparison_pair(self, nominator_hist, denominator_hist,
                            location, ratio_location,
                            nominator_args, denominator_args):
        self.add_hist(nominator_hist, location=location, as_denominator=False, **nominator_args)
        self.add_hist(denominator_hist, location=location, as_denominator=True, **denominator_args)
        self.draw_ratio_hists(location=ratio_location)

    def draw_matrix(self, rm_np, x_axis_label="", y_axis_label="", show_number=False,
                    number_fontsize=3,
                    **kwargs):
        self.set_current_axis((0, 0))
        # hep.hist2dplot(rm_np, norm=mcolors.LogNorm(), ax=self.current_axis, **kwargs)
        hep.hist2dplot(rm_np, ax=self.current_axis, **kwargs)

        if (rm_np[1][-1] - rm_np[1][0] >= 5e2 and
                rm_np[2][-1] - rm_np[2][0] >= 5e2):
            self.current_axis.set_xscale("log")
            self.current_axis.set_yscale("log")

            # self.current_axis.set_xlim(0.2, rm_np[1][-1])
            # self.current_axis.set_ylim(0.2, rm_np[2][-1])

        if show_number:
            for x in range(len(rm_np[1]) - 1):
                for y in range(len(rm_np[2]) - 1):
                    c = rm_np[0][x][y]
                    if math.isnan(c):
                        continue
                    if c < 0.01:
                        continue
                    x_half_width = (rm_np[1][x + 1] - rm_np[1][x]) / 2
                    y_half_width = (rm_np[2][y + 1] - rm_np[2][y]) / 2
                    self.current_axis.text(rm_np[1][x] + x_half_width,
                                           rm_np[2][y] + y_half_width,
                                           f'{c:.2f}', va='center', ha='center',
                                           fontsize=number_fontsize, color='red')

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
            label = self.errorbar_kwargs[index].get("label", None)
            if "fill_between_only" in self.errorbar_kwargs[index]:
                del self.errorbar_kwargs[index]['fill_between_only']
                temp_kwargs = {}
                if "color" in self.errorbar_kwargs[index]:
                    temp_kwargs['color'] = self.errorbar_kwargs[index]['color']
                if "label" in self.errorbar_kwargs[index]:
                    temp_kwargs['label'] = self.errorbar_kwargs[index]['label']
                    del self.errorbar_kwargs[index]['label']
                handle = self.current_axis.fill_between(x=x_value, y1=y_value-y_error, y2=y_value+y_error,
                                                        alpha=0.1,
                                                        **temp_kwargs)
                self.current_axis.errorbar(x_value, y_value, **self.errorbar_kwargs[index])
            else:
                self.current_axis.errorbar(x_value, y_value, xerr=x_error, yerr=y_error,
                                           **self.errorbar_kwargs[index])

                handles, labels = self.current_axis.get_legend_handles_labels()
                handle = handles[-1]
            if label:
                self.legend_handles[self.errorbar_loc[index]].append(handle)
                self.legend_labels[self.errorbar_loc[index]].append(label)

    def update_legend(self, location=(0, 0)):
        handles, labels = self.get_axis(location=location).get_legend_handles_labels()
        #label = self.errorbar_kwargs[index].get("label", None)
        self.legend_handles[location].append(handles[-1])
        self.legend_labels[location].append(labels[-1])

    def draw_vlines(self, vlines, location=(0, 0), **kwargs):
        for line in vlines:
            self.get_axis(location=location).axvline(x=line, **kwargs)

    def save_fig(self, out_name='', out_sub_dirs=''):
        if out_sub_dirs:
            out_file_name = self.base_output_dir + out_sub_dirs + out_name + ".pdf"
        else:
            out_file_name = self.base_output_dir + "/" + out_name + ".pdf"
        #print(f"save plot... {out_file_name}")
        self.fig.savefig(out_file_name)


