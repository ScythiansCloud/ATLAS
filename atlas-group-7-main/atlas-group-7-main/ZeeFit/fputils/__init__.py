import time
import concurrent
import numpy as np
import uproot as up
import matplotlib.pyplot as plt
import ROOT

from ipywidgets import IntProgress, HTML, VBox
from IPython.display import display

class NextTreeInterruption(Exception):
    pass

class Analysis:
    """
    An analysis object is used to keep track of samples, xsecs (normlization
    in general), luminosity, files, plotting labels and meta-data.
    """

    def __init__(self, luminosity_fb):
        """
        Creates a new Analysis object. The given luminosity (in 1/femtobarn)
        is used to normalize the MC samples.
        """
        self.luminosity = luminosity_fb
        self.names = []
        self.labels = []
        self.files = []
        self.expected_events = []  # Contains None for Data samples
        self.meta_data = []
        self.histogram_ids = []
        self.next_histogram_id = 0

    def add_sample(self, name, label, files, xsec_fb=None, filter_efficiency=1,
                   branching_ratio=1, meta_data=None):
        if name in self.names:
            sample_id = self.names.index(name)
            existing_label = self.labels[sample_id]
            if label != existing_label:
                raise ValueError(f"Cannot add sample with name {name!r}, "
                                 f"label conflict: {label!r} should be "
                                 f"{existing_label!r}")
            self.histogram_ids.append(self.histogram_ids[sample_id])
        else:
            self.histogram_ids.append(self.next_histogram_id)
            self.next_histogram_id += 1

        self.names.append(name)
        self.labels.append(label)
        self.files.append(files)
        self.meta_data.append(meta_data)

        n_events = self._calc_n_events(xsec_fb, filter_efficiency,
                                       branching_ratio)
        self.expected_events.append(n_events)


    def _calc_n_events(self, xsec, filter_efficiency, branching_ratio):
        if xsec is None:
            return None

        return self.luminosity * xsec * filter_efficiency * branching_ratio

    def _inverted_paths(self):
        """
        Returns a dict with unique paths as keys. The items are lists of
        sample indices which need include this path.
        """
        inverse = {}
        for i, paths in enumerate(self.files):
            for path in paths:
                index_list = inverse.setdefault(path, [])
                index_list.append(i)
        return inverse

    def TH1F_stack(self, nbins, xmin, xmax):
        """
        Returns a list of identical histogram, one histogram for each unique
        sample name.
        """
        
        stack = [ROOT.TH1F("", "", nbins, xmin, xmax)
                 for i in range(self.next_histogram_id)]
        for histogram_id, label in zip(self.histogram_ids, self.labels):
            stack[histogram_id].SetTitle(label)

        return stack

    def dataframe_iterate(self, tree_name="eventTree", branches=None, fraction=1.0):
        """
        Iterates over all given root files, reads the trees in a chunked manner
        and yields a Pandas DataFrame for each chunk. For each chunk the
        method yields:
            data_frame, normalization, meta_data, histogram_index

        If the variable fraction is provided (<=1), the fraction of events from each file will be read.
        """
        inv_paths = self._inverted_paths()
        trees = [up.open(path)[tree_name] for path in inv_paths]
        
        processed_events = 0
        total_events = int(sum([fraction*tree.num_entries*len(sample_ids) 
                            for tree, sample_ids
                            in zip(trees, inv_paths.values())]))

        progress = IntProgress(min=0, max=total_events, value=0)

        label = HTML()
        box = VBox(children=[label, progress])
        display(box)
        
        start = time.time()
        eta = None

        try:
            # Loop over files and samples
            for tree_i, (file_name, sample_ids) in enumerate(inv_paths.items()):
                tree = trees[tree_i]
                root_file = up.open(file_name)
                info_tree = root_file["infoTree"]
                sow_name = "nEventsProcessedStage0_MCEventWeightSum"
                if sow_name in info_tree:
                    sum_of_weights =sum(info_tree.arrays('nEventsProcessedStage0_MCEventWeightSum')['nEventsProcessedStage0_MCEventWeightSum'])
                else:
                    sum_of_weights = None

                try:
                    # Loop over chunks in file
                    for df_count, df in enumerate(tree.iterate(library='pd',filter_name=branches,entry_stop=fraction*tree.num_entries)):
                        for sample_id in sample_ids:
                            normalization = self.expected_events[sample_id]
                            if normalization is not None:
                                # Assume MC
                                normalization /= sum_of_weights
                            else:
                                # Assume Data
                                normalization = 1

                            yield (df, normalization, self.meta_data[sample_id],
                                   self.histogram_ids[sample_id])
                        
                            processed_events += len(df)
                            progress.value = processed_events
                            label.value = f"Progress: {processed_events / total_events * 100:.0f}% " \
                                   f"({processed_events:,d} / {total_events:,d} events)"

                            elapsed = time.time() - start
                            events_per_second = processed_events / elapsed
                            current_eta = (total_events - processed_events) / events_per_second

                            # Show ETA only after 10 seconds and after loading the
                            # second file
                            if elapsed > 10 and (tree_i > 1 or len(trees) ==1):
                                if eta is None:
                                    eta = current_eta

                                # Low pass filter of ETA
                                eta = (current_eta + 9 * eta) / 10
                                label.value += f", ETA: {eta:.0f}s"
                except NextTreeInterruption as e:
                    pass
                        
        except:
            progress.bar_style = 'danger'
            raise
        else:
            progress.bar_style = 'warning' if fraction<1.0 else 'success' 
            if fraction<1.0:
                progress.value = processed_events
                label.value = f"Progress: {processed_events:,d} events (debug mode)"
            else:
                progress.value = processed_events
                label.value = f"Progress: 100% ({processed_events:,d} events)"
        

def apply_calib(E, eta, alpha, eta_edges):
    """
    Applies the energy callibration by using the the appropriate values from the
    alpha variable depending on eta. Please note that the arguments and return
    value are arrays.
    """
    alpha_array = np.zeros(len(E))
    
    eta_idx = np.searchsorted(eta_edges, eta.abs()) - 1
    for idx, a in enumerate(alpha):
        # We apply different alpha correction factors, one at a time, only to
        # events in the appropriate eta bin
        alpha_array[eta_idx == idx] = a
        
    # The entries in alpha_array are the alpha factors to the corresponding
    # events in the E and eta array.
        
    return E / (1 + alpha_array)

def fill_histo(array, histogram, weights=None):
    """
    Fill all entries of the given array into the histogram in a performant
    way.

    >>> import ROOT
    >>> histo = ROOT.TH1F("h", "", 4, 0, 4)
    >>> data = np.array([0, 1, 3, 0, 3])
    >>> fill_histo(data, histo)
    >>> histo.GetBinContent(1)
    2.0
    >>> histo.GetBinContent(2)
    1.0
    >>> histo.GetBinContent(3)
    0.0
    >>> histo.GetBinContent(4)
    2.0
    >>> histo.GetEntries()
    5.0
    """
    bins = edges(histogram)
    
    bin_contents, _ = np.histogram(array, bins=bins, weights=weights)
    if weights is None:
        bin_errors2 = bin_contents
    else:
        bin_errors2, _ = np.histogram(array, bins=bins, weights=weights**2)

    entries = histogram.GetEntries()
    for i, (bin_content, bin_error2) in enumerate(zip(bin_contents, bin_errors2), 1):
        old_error = histogram.GetBinError(i)
        histogram.SetBinContent(i, histogram.GetBinContent(i) + bin_content)
        histogram.SetBinError(i, np.sqrt(bin_error2 + old_error**2))
        
    histogram.SetEntries(entries + len(array))

def edges(histogram):
    """
    Extracts the bin edges of a ROOT histogram and returns them as a numpy
    histogram. The length of the returned array is thus number of bins + 1.

    >>> import ROOT
    >>> histo = ROOT.TH1F("h", "", 4, 0, 5)
    >>> edges(histo)
    array([0.  , 1.25, 2.5 , 3.75, 5.  ])
    """
    x = histogram.GetXaxis()
    bins = [x.GetBinLowEdge(1)]
    bins += [x.GetBinUpEdge(i + 1) for i in range(x.GetNbins())]
    bins = np.array(bins)

    return bins

def dump_histo(histogram, uncertainty=True):
    """
    Extracts the bin contents and their uncertainties from the given histogram
    and returns them as a tuple of numpy arrays.

    >>> import ROOT
    >>> histo = ROOT.TH1F("h", "", 3, 0, 3)
    >>> histo.Fill(0.1)
    1
    >>> histo.Fill(2.1, 2);
    3
    >>> dump_histo(histo)
    (array([1., 0., 2.]), array([1., 0., 2.]))

    If `uncertainty` is `False` (default `True), the method returns only the
    bin contents as a numpy array.

    >>> dump_histo(histo, uncertainty=False)
    array([1., 0., 2.])
    """
    n_bins = histogram.GetNbinsX()
    bin_contents = [histogram.GetBinContent(i + 1) for i in range(n_bins)]

    if uncertainty:
        uncertainties = [histogram.GetBinError(i + 1) for i in range(n_bins)]
        return np.array(bin_contents), np.array(uncertainties)

    return np.array(bin_contents)

def hist_points(histogram, axes=None):
    """
    Plots the given histogram as a separate line histogram. The `axes`
    argument is used to plot if not omitted.

    Returns the axes.
    
    >>> import ROOT
    >>> first = ROOT.TH1F("", "First", 10, 0, 10)
    >>> second = ROOT.TH1F("", "Second", 10, 0, 10)
    >>> axes = hist_points(first)
    >>> axes = hist_points(second)
    """
    if axes is None:
        fig, axes = plt.subplots()

    bins = edges(histogram)

    bin_contents, uncertainty = dump_histo(histogram)
    bin_centers = (bins[1:] + bins[:-1]) / 2
    bin_widths = (bins[1:] - bins[:-1]) / 2
    axes.errorbar(bin_centers, bin_contents, uncertainty,
                  bin_widths, markersize=4,color='k', fmt='o',
                  label=histogram.GetTitle())

    return axes

def hist_line(*histograms, axes=None,**args):
    """
    Plots the given histograms as a separate line histograms. The `axes`
    argument is used to plot if not omitted.

    Returns the axes.
    
    >>> import ROOT
    >>> first = ROOT.TH1F("", "First", 10, 0, 10)
    >>> second = ROOT.TH1F("", "Second", 10, 0, 10)
    >>> axes = hist_list(first, second)
    """
    if axes is None:
        fig, axes = plt.subplots()

    if len(histograms) == 0:
        return axes

    for histogram in histograms:
        bins = edges(histogram)

        bin_contents, uncertainty = dump_histo(histogram)
        axes.hist(bins[:-1], bins=bins, weights=bin_contents,
                  histtype='step',label=histogram.GetTitle(),**args)

    return axes

def hist_stack(*histograms, axes=None):
    """
    Plots the given histograms as a stacked histogram. The `axes`
    argument is used to plot if not omitted. An exception is raised if the
    histograms have differing bin edges.

    Returns the axes.
    
    >>> import ROOT
    >>> first = ROOT.TH1F("", "First", 10, 0, 10)
    >>> second = ROOT.TH1F("", "Second", 10, 0, 10)
    >>> axes = hist_stack(first, second)
    """
    if axes is None:
        fig, axes = plt.subplots()

    if len(histograms) == 0:
        return axes

    bottom = 0
    bins = edges(histograms[0])

    total_uncertainty = np.zeros(len(bins) - 1)
    for i, histogram in enumerate(histograms):
        if (bins != edges(histogram)).any():
            raise Exception("Bins of {histogram!r} are not aligned.")

        label = histogram.GetTitle() if len(histograms) > 1 else None

        bin_contents, uncertainty = dump_histo(histogram)
        axes.hist(bins[:-1], bins=bins, weights=bin_contents, bottom=bottom,
                 label=label, histtype='stepfilled', zorder=-i)
        bottom += bin_contents
        total_uncertainty += uncertainty**2

    total_uncertainty = np.sqrt(total_uncertainty )
    axes.hist(bins[:-1], bins=bins,
              weights=2 * total_uncertainty,
              bottom=bottom - total_uncertainty,
              histtype='stepfilled', zorder=0,
              fill=False, hatch='/////',
              linewidth=0, edgecolor="#666666")

    return axes

def reverse_legend(axes=None):
    if axes is None:
        axes = plt.gca()

    handles, labels = axes.get_legend_handles_labels()
    
    #Plot Lines instead of rectangles for line_histos:
    import matplotlib.legend_handler
    handles=[matplotlib.legend_handler.Line2D([],[],color=handle.get_edgecolor())
             if isinstance(handle,matplotlib.patches.Polygon) and handle.get_edgecolor()!=(0,0,0,0)
             else handle
             for handle in handles]
    axes.legend(handles[::-1], labels[::-1], frameon=False, loc=1)

def dataframe_iterate(paths, tree_name="eventTree", name=False, fraction=1.0):
    """
    Iterates over all given root files, reads the trees in a chunked manner
    and yields a Pandas DataFrame for each chunk.

    If the `name` argument is `True`, the method returns a tuple:
        current file name, data frame
    
    If the variable fraction is provided (<=1), the fraction of events from each file will be read.
    """
    trees = [up.open(path)[tree_name] for path in paths]
    
    total_events = int(sum([fraction*tree.num_entries for tree in trees]))
    processed_events = 0
    eta = None
    
    progress = IntProgress(min=0, max=total_events, value=0)

    label = HTML()
    box = VBox(children=[label, progress])
    display(box)
    
    start = time.time()

    try:
        # Loop over file
        for i, (file_name, tree) in enumerate(zip(paths, trees), 1):
            # Loop over chunks in file
            for df in tree.iterate(library='pd',entry_stop=fraction*tree.num_entries):
                if name:
                    yield file_name, df
                else:
                    yield df
                processed_events += len(df)
                progress.value = processed_events
                label.value = f"Progress: {processed_events / total_events * 100:.0f}% " \
                       f"({processed_events:,d} / {total_events:,d} events)"

                elapsed = time.time() - start
                events_per_second = processed_events / elapsed
                current_eta = (total_events - processed_events) / events_per_second

                # Show ETA only after 10 seconds and after loading the
                # second file
                if elapsed > 10 and (i > 1 or len(trees) ==1):
                    if eta is None:
                        eta = current_eta

                    # Low pass filter of ETA
                    eta = (current_eta + 9 * eta) / 10
                    label.value += f", ETA: {eta:.0f}s"
                
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'warning' if fraction<1.0 else 'success' 
        if fraction<1.0:
            progress.value = processed_events
            label.value = f"Progress: {processed_events:,d} events (debug mode)"
        else:
            progress.value = processed_events
            label.value = f"Progress: 100% ({processed_events:,d} events)"


class ZeeModel:
    """
    Fit model for Zee calibration. The model describes a combination of signal
    and background events. The background is modelled as an exponentally
    falling spectrum. The signal is a convuluation of a breit wigner and the
    crystal ball function.
    """

    def __init__(self, histo, xMin, xMax):
        """
        Creates a new fit model to be fitted to the given histogram. THe fit
        range is limited by xMin and xMax.
        """
        self.x = ROOT.RooRealVar("x", "Mll", xMin, xMax)
        self.data = ROOT.RooDataHist("data", "data", ROOT.RooArgList(self.x), histo)

        # Signal function
        self.meanBW  = ROOT.RooRealVar("meanBW", "meanBW",
                                       self.data.mean(self.x),
                                       self.data.mean(self.x) - 3,
                                       self.data.mean(self.x) + 3)
        self.widthBW = ROOT.RooRealVar("widthBW","widthBW", 2.5, 0, 10)
        self.breitWigner = ROOT.RooBreitWigner("BreitWigner", "BreitWigner",
                                               self.x,
                                               self.meanBW,
                                               self.widthBW)

        self.meanCB  = ROOT.RooRealVar("meanCB" , "meanCB" ,   0.6840,   0,     10)
        self.widthCB = ROOT.RooRealVar("widthCB", "WidthCB",   2,        0.01,  10)
        self.nCB     = ROOT.RooRealVar("nCB",     "nCB"    , 115.,      10,    400)
        self.alphaCB = ROOT.RooRealVar("alphaCB", "alphaCB",   0.542,    0,     20)

        self.CBshape  = ROOT.RooCBShape("CBshape", "CBshape",
                                        self.x,
                                        self.meanCB,
                                        self.widthCB,
                                        self.alphaCB,
                                        self.nCB)
        self.ZeeModel = ROOT.RooFFTConvPdf("ZeeModel", "BreitWigner (X) crystalball",
                                           self.x,
                                           self.breitWigner,
                                           self.CBshape)
        self.ZeeModel.setBufferFraction(0.25)


        # Background function
        self.etau = ROOT.RooRealVar("etau", "etau", -0.1, -1, 0)
        self.expo = ROOT.RooExponential("expoBkg", "expoBkg", self.x, self.etau)

        self.N = self.data.sumEntries();
        self.nsig = ROOT.RooRealVar("nsig", "number of signal events", self.N, 0, 2 * self.N)
        self.nbkg = ROOT.RooRealVar("nbkg", "number of background events", 0, 0, 2 * self.N)
        self.eSig = ROOT.RooExtendPdf("eSig", "eSig", self.ZeeModel, self.nsig)
        self.eBkg = ROOT.RooExtendPdf("eBkg", "eBkg", self.expo, self.nbkg)

        self.model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(self.eSig, self.eBkg))
        
    def fit_to_data(self):
        """
        Performs the actual fit and returns the fitted mean value of the breit
        wigner.
        """
        if (self.data.sumEntries()>0.):
            self.result = self.model.fitTo(self.data,
                                 ROOT.RooFit.SumW2Error(True),
                                 ROOT.RooFit.Save(),
                                 ROOT.RooFit.PrintLevel(-1),
                                 ROOT.RooFit.PrintEvalErrors(-1),
                                 ROOT.RooFit.Warnings(0))

            # self.plot = x.frame()
            # self.data.plotOn(self.plot, ROOT.RooFit.MarkerSize(1.2))
            # self.model.plotOn(self.plot,
            #              ROOT.RooFit.LineColor(2),
            #              ROOT.RooFit.LineStyle(2),
            #              ROOT.RooFit.LineWidth(4),
            #              ROOT.RooFit.DrawOption("L"))
            #
            # self.result.plotOn(self.plot, "meanBW", "widthBW")
            # self.plot.Draw()
            #
            # self.result.Print()

            self.mean = self.result.floatParsFinal().find("meanBW")
            return self.mean.getVal()
        else:
            return 0.
