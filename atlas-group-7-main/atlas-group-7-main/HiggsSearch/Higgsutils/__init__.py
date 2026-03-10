import time
import concurrent
import numpy as np
import uproot as up
import matplotlib.pyplot as plt
import ROOT
import collections
import awkward as ak
import numba as nb
from ipywidgets import IntProgress, HTML, VBox
from IPython.display import display
import atlasify as atl
atl.ATLAS = "FP2 - ATLAS"  # Do not remove

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
        Returns an ordered dictionary of identical histogram, one histogram for each unique
        sample name.
        """
        
        ordDict = collections.OrderedDict([(name,ROOT.TH1F("", "", nbins, xmin, xmax))
                                         for i,name in enumerate(self.names)])
        for name, label in zip(self.names, self.labels):
            ordDict[name].SetTitle(label)

        return ordDict

    def akward_iterate(self, tree_name="eventTree", branches=None, fraction=1.0):
        """
        Iterates over all given root files, reads the trees in a chunked manner
        and yields an Akward Array for each chunk. For each chunk the
        method yields:
            akward_array, meta_data, histogram_index

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

                try:
                    # Loop over chunks in file
                    for df_count, df in enumerate(tree.iterate(filter_name=branches,entry_stop=fraction*tree.num_entries)):
                        for sample_id in sample_ids:

                            yield (df, self.meta_data[sample_id],
                                   self.names[sample_id])
                        
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

    def ScanHistograms(self,hist_stack, isUp=False,x_label='',y_label=''):
        background_names=[]
        background_hist=[]
        signal_names=[]
        signal_hist=[]
        for i,name in enumerate(self.names):
            if(self.meta_data[i]['mc']):
                if(self.meta_data[i]['signal']):
                    signal_names.append(name)
                    signal_hist.append(hist_stack[name])
                else:
                    background_names.append(name)
                    background_hist.append(hist_stack[name])


        signal_combined=sum_histograms(signal_hist)
        background_combined=sum_histograms(background_hist)

        hist_significance=background_combined.Clone("significance")
        significance=[]

        nbins=background_combined.GetNbinsX()
        for i in range(1,nbins+1):
            s = 1.;
            b = 1.;
            if(isUp):
                s = signal_combined.Integral(0,i)
                b = background_combined.Integral(0,i) 

            else:
                s = signal_combined.Integral(i,nbins) 
                b = background_combined.Integral(i,nbins);

            hist_significance.SetBinError(i,0)
            if(b<0.0001):
                hist_significance.SetBinContent(i,0);
            else:
                z = s/pow(b,0.5);
                hist_significance.SetBinContent(i,z);
        if(signal_combined.Integral()>0):
            signal_combined.Scale(1./signal_combined.Integral())
        if(background_combined.Integral()):
            background_combined.Scale(1./background_combined.Integral())
        signal_combined.SetTitle('signal')
        background_combined.SetTitle('background')    

        ax=hist_line(signal_combined,color='r')
        bins = edges(hist_significance)
        bin_centers = (bins[1:] + bins[:-1]) / 2
        sigLine,=plt.plot(bin_centers,dump_histo(hist_significance, uncertainty=False),color='g',marker='x',ls='None')
        print(sigLine)
        legend = plt.legend([sigLine], [r"s/$\sqrt{b}$"+(' (upper)' if isUp else ' (lower)')], loc=(0.4,0.82))
        plt.gca().add_artist(legend)

        hist_line(background_combined,color='k',axes=ax)
        atl.atlasify(subtext='Higgs search')
        ax.legend()
        ax.set_xlabel(x_label, ha='right', x=0.95)
        ax.set_ylabel(y_label, ha='right', x=0.95)
        return ax

    def PValue(self,hist_stack,observed=False,scan_window=10):
        background_names=[]
        background_hist=[]
        signal_masses=[]
        signal_names=[]
        signal_hist=[]
        for i,name in enumerate(self.names):
            if(self.meta_data[i]['mc']):
                if(self.meta_data[i]['signal']):
                    signal_names.append(name)
                    signal_hist.append(hist_stack[name])
                    signal_masses.append(float(self.meta_data[i]['HiggsMass']))
                else:
                    background_names.append(name)
                    background_hist.append(hist_stack[name])

        background_combined=sum_histograms(background_hist)
        hdata=hist_stack['data']
        scan_window=5

        pvalues=[]
        for i,name in enumerate(signal_names):
            mass=signal_masses[i]
            nbkg = background_combined.Integral(background_combined.FindBin(mass-(scan_window/2)),background_combined.FindBin(mass+(scan_window/2)));
            ndata=0
            if(observed):
                ndata = hdata.Integral(hdata.FindBin(mass-(scan_window/2.)),hdata.FindBin(mass+(scan_window/2.)));
            else:
                Expected_combined=sum_histograms(background_hist+[signal_hist[i]])
                ndata = Expected_combined.Integral(Expected_combined.FindBin(mass-(scan_window/2.)),Expected_combined.FindBin(mass+(scan_window/2.)));

            Poisson = ROOT.TF1("f","TMath::Poisson(x,[0])",0,200);
            Poisson.SetParameter(0,nbkg);
            sampling_grid_size=1000
            xval=np.array(sampling_grid_size*[0.]);
            wval=np.array(sampling_grid_size*[0.]);
            Poisson.CalcGaussLegendreSamplingPoints(sampling_grid_size,xval,wval,1e-15);
            pval = Poisson.IntegralFast(sampling_grid_size,xval,wval,ndata,10000);
            pvalues.append(pval)

            Gaussian = ROOT.TF1("f2","TMath::Gaus(x)",-100,100);
            xq=np.array([1-pval])
            yq=np.array([0.])
            Gaussian.GetQuantiles(1, yq, xq);

            print("---------------------------------------------------------------------------")
            print(f"- Mass value {mass:.2f}, p-value {pval:.2e}; n_bkg {nbkg:.1f}; {'n_observed' if observed else 'n_bkg+sig:'} {ndata:.1f}")
            print(f"     Quantile of {xq[0]:.2f} is {yq[0]:.2f} Sigma ")
            print("---------------------------------------------------------------------------")
        fig, ax = plt.subplots()
        plt.plot(signal_masses,pvalues)
        ax.set_xlabel(r'$m_H$ [GeV]', ha='right', x=0.95)
        ax.set_ylabel('p-value'+(' (observed)' if observed else ' (expected)'), ha='left', x=0.95)
        plt.plot([min(signal_masses),max(signal_masses)],[3e-7,3e-7],color='r')
        atl.atlasify(subtext="Higgs search",axes=ax)

        return pvalues,ax

    def FractionFit(self,min,max,hist_stack):
        background_names=[]
        background_hist=[]
        signal_names=[]
        signal_hist=[]
        for i,name in enumerate(self.names):
            if(self.meta_data[i]['mc']):
                if(self.meta_data[i]['signal']):
                    signal_names.append(name)
                    signal_hist.append(hist_stack[name])
                else:
                    background_names.append(name)
                    background_hist.append(hist_stack[name])

        background_combined=sum_histograms(background_hist)
        print(signal_names)
        signal_combined=sum_histograms(signal_hist)
        hdata=hist_stack['data']
        mc_samples = ROOT.TObjArray(2)
        mc_samples.Add(background_combined)
        mc_samples.Add(signal_combined)

        FractionFitter=ROOT.TFractionFitter(hdata, mc_samples, "V");
        FractionFitter.Constrain(0,0.,2.);             
        FractionFitter.Constrain(0,0.,2.);             
        FractionFitter.SetRangeX(hdata.FindBin(min),hdata.FindBin(max));
        status = FractionFitter.Fit()
        print(f"Fit status: {status}")
        value1,error1=np.array([0.]),np.array([0.])
        FractionFitter.GetResult(0, value1, error1);
        value2,error2=np.array([0.]),np.array([0.])
        FractionFitter.GetResult(1, value2, error2);
        if(status==0):

            fig, ax = plt.subplots()
            result = FractionFitter.GetPlot();
            ndata = hdata.Integral(hdata.FindBin(min),hdata.FindBin(max));
            nbg   = background_combined.Integral(hdata.FindBin(min),hdata.FindBin(max));
            nsig  = signal_combined.Integral(hdata.FindBin(min),hdata.FindBin(max));
            nfit = result.Integral(-1,-1);
            bkgSF=float((value1*nfit)/nbg)
            sigSF=float((value2*nfit)/nsig)

            hist_points(hdata,axes=ax);
            label=f"{float(value1):.1f} bkg+{float(value2):.1f} sig"
            result.SetTitle('Fit (bkg+sig)')
            hist_line(result,axes=ax,lw=2)
            fitted_signal=signal_combined.Clone()
            fitted_signal.Scale(sigSF)
            fitted_signal.SetTitle(f"Higgs sig * {sigSF:.2f}")
            hist_line(fitted_signal,axes=ax,lw=2,color='r')
            fitted_background=background_combined.Clone()
            fitted_background.Scale(bkgSF)
            fitted_background.SetTitle(f"backgrounds * {bkgSF:.2f}")
            hist_line(fitted_background,axes=ax,lw=2,color='g',alpha=0.5)
            ax.set_xlabel(r'$m_H$ [GeV]', ha='right', x=0.95)
            ax.set_ylabel('Events', ha='left', x=0.95)
            atl.atlasify(subtext="Higgs search",axes=ax)

            print("-------------------------------------------------------------------")
            print(f"FractionFit result: bg {float(value1):.2f}+/-{float(error1):.2f} sig {float(value2):.2f}+/-{float(error2):.2f}")
            print(f"ndata {ndata:.1f} nbg {nbg:.1f} nsig {nsig:.1f}")
            print(f"bgSF {bkgSF:.2f}; sigSF {sigSF:.2f}")
            print("-------------------------------------------------------------------")
            return ax
def sum_histograms(histograms):
    """
    Returns ROOT.TH1 histogram with the sum of the input ROOT.TH1 histograms.
    Errors are added in quadratures.
    """
    result = histograms[0].Clone()
    result.Reset()
    result.SetTitle("")
    for histogram in histograms:
        result.Add(histogram)
    
    return result

def fill_histo(array, histogram, weights=None):
    """
    Fill all entries of the given an array/akward_array into the histogram in a performant
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
    if(isinstance(array,ak.Array)):
        array=np.array(array)
    if(isinstance(weights,ak.Array)):
        weights=np.array(weights)

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

def hist_points(histogram, axes=None,**args):
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
    if('color' not in args):
        args['color']='k'
    if('markersize' not in args):
        args['markersize']=4
    if('label' not in args):
        args['label']=histogram.GetTitle()
    if('fmt' not in args):
        args['fmt']='o'

    if axes is None:
        fig, axes = plt.subplots()

    bins = edges(histogram)

    bin_contents, uncertainty = dump_histo(histogram)
    bin_centers = (bins[1:] + bins[:-1]) / 2
    bin_widths = (bins[1:] - bins[:-1]) / 2
    axes.errorbar(bin_centers, bin_contents, uncertainty,
                  bin_widths,**args)

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

    total_uncertainty = np.sqrt(total_uncertainty)

    axes.hist(bins[:-1], bins=bins,
              weights=2 * total_uncertainty,
              bottom=bottom - total_uncertainty,
              histtype='stepfilled', zorder=0,
              fill=False, hatch='/////',
              linewidth=0, edgecolor="#666666")

    return axes

def reverse_legend(axes=None,**args):
    if axes is None:
        axes = plt.gca()

    handles, labels = axes.get_legend_handles_labels()
    
    #Plot Lines instead of rectangles for line_histos:
    import matplotlib.legend_handler
    handles=[matplotlib.legend_handler.Line2D([],[],color=handle.get_edgecolor())
             if isinstance(handle,matplotlib.patches.Polygon) and handle.get_edgecolor()!=(0,0,0,0)
             else handle
             for handle in handles]
    if('frameon' not in args):
        args['frameon']=False
    if('loc' not in args):
        args['loc']=1
    axes.legend(handles[::-1], labels[::-1],**args)
    
def akward_iterate(paths, tree_name="eventTree", name=False, fraction=1.0):
    """
    Iterates over all given root files, reads the trees in a chunked manner
    and yields an Akward Array for each chunk.

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
            for df in tree.iterate(entry_stop=fraction*tree.num_entries):
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




@nb.jit(forceobj=False,nopython=True)
def calc_invmass(pts,etas,phis):
    # To students: No need to modify this function!
    #
    # This function takes numpy arrays of pts, etas and phis and calculates the invariant mass of the provided particles
    # Needed to code it explicitely to be able to run it within numba (in a performat way)
    theta=2*np.arctan(np.exp(-etas))
    p=pts/np.sin(theta)
    pz=p*np.cos(theta)
    px=p*np.sin(theta)*np.cos(phis)
    py=p*np.sin(theta)*np.sin(phis)
    sumpz=np.sum(pz)
    sumpx=np.sum(px)
    sumpy=np.sum(py)
    sumE=np.sum(p)
    minv=sumE**2 - sumpz**2 - sumpx**2 - sumpy**2
    return np.sqrt(minv)

@nb.jit(forceobj=False,nopython=True)
def invariant_masses(events):
    # To students: No need to modify this function!
    #
    # This function accepts akward array and calculates for each row (event) the following quantities:
    # -> m12 -- invariant mass of the oposite charge same flavour lepton pair with the closest value to Z mass in the event
    # -> m34 -- invariant mass of the same-flavour opposite charge leptons fornming next most likely Z-candidate (leptons from the first pair are excluded)
    # -> mllll -- invariant mass of the four leptons forming the two most likely Z-candidates
    # -> i12, j12 -- indices of the leptons forming the leading (by mass-agreement) Z-candidate
    # -> i34, j34 -- indices of the leptons forming the sub-leading Z-candidate
    #
    # These results are returned as numpy arrays, with one row corresponding to one event
    # This function is numba-compatible improving the computation time
    # 
    zmass=91.18
    out12 = np.empty(len(events), np.float64)
    out34 = np.empty(len(events), np.float64)
    out4l = np.empty(len(events), np.float64)
    i12 = np.empty(len(events), np.int32)
    j12 = np.empty(len(events), np.int32)
    i34 = np.empty(len(events), np.int32)
    j34 = np.empty(len(events), np.int32)
    #leptIndicesOut = np.empty(len(events), np.array)
    selectedLeptonIndices = np.array([-1,-1,-1,-1])
    
    num_pairs = 0
    for event in events:
        minMass12=140000;
        minMass34=140000;
        out12[num_pairs] = -1.0;
        out34[num_pairs] = -1.0;
        out4l[num_pairs] = -1.0;
        # Let's loop over all lepton pairs and find the one with invariant mass cosest to Z mass
        for i in range(len(event.lep_pt)):
            for j in range(i + 1, len(event.lep_pt)):
                if((event.lep_charge[i])*(event.lep_charge[j])<0. and abs((event.lep_id[i])-(event.lep_id[j]))<1.):
                    mass=np.sqrt(2*event.lep_pt[i]*event.lep_pt[j]*(np.cosh(event.lep_eta[i] - event.lep_eta[j]) - np.cos(event.lep_phi[i] - event.lep_phi[j])))
                    if(abs(mass-zmass)<minMass12):
                        minMass12=abs(mass-zmass)
                        out12[num_pairs] = mass
                        selectedLeptonIndices[0]=i # i12
                        selectedLeptonIndices[1]=j # j12
        # Let's loop over all lepton pairs, excluding ones from the previous step and find the one with invariant mass cosest to Z mass
        for i in range(len(event.lep_pt)):
            for j in range(i + 1, len(event.lep_pt)):
                if((i not in selectedLeptonIndices[:2]) and (j not in selectedLeptonIndices[:2])):
                    if((event.lep_charge[i])*(event.lep_charge[j])<0. and abs((event.lep_id[i])-(event.lep_id[j]))<1.):
                        mass=np.sqrt(2*event.lep_pt[i]*event.lep_pt[j]*(np.cosh(event.lep_eta[i] - event.lep_eta[j]) - np.cos(event.lep_phi[i] - event.lep_phi[j])))
                        if(abs(mass-zmass)<minMass34):
                            minMass34=abs(mass-zmass)
                            out34[num_pairs] = mass
                            selectedLeptonIndices[2]=i # i34
                            selectedLeptonIndices[3]=j # j34
        i12[num_pairs]=selectedLeptonIndices[0]
        j12[num_pairs]=selectedLeptonIndices[1]
        i34[num_pairs]=selectedLeptonIndices[2]
        j34[num_pairs]=selectedLeptonIndices[3]
        # TUTORS
        """
        In the ROOT macro, m4l is calculated as inv mass of the leading 4-leptons.
        Calculating it from Z-decay product candidates seems to be a better choice.
        """ 
        # END TUTORS
        if(out12[num_pairs]>0 and out34[num_pairs]>0):
            out4l[num_pairs]=calc_invmass(np.array(event.lep_pt)[selectedLeptonIndices],
                                np.array(event.lep_eta)[selectedLeptonIndices],
                               np.array(event.lep_phi)[selectedLeptonIndices])
        num_pairs += 1
    return (out12,out34,out4l,i12,j12,i34,j34)

def add_invariant_masses(akw):
    # To students: No need to modify this function!
    #
    # This function dresses akward array with the results of 'invariant_masses' function. 
    # This function itself is not running in a performant way, but calls invariant_masses which does (with numba).
    #
    m12,m34,mllll,i12,j12,i34,j34=invariant_masses(akw)
    akw['m12']=ak.Array(m12)
    akw['m34']=ak.Array(m34)
    akw['mllll']=ak.Array(mllll)
    akw['i12']=ak.Array(i12)
    akw['j12']=ak.Array(j12)
    akw['i34']=ak.Array(i34)
    akw['j34']=ak.Array(j34)
    return akw

def sort_leptons_by_pt(akw_array):
    # To students: No need to modify this function!
    #
    # Function accepts akward array with fields 'lep_pt', 'lep_eta','lep_phi','lep_E',`lep_ptiso`,`lep_etiso`,`lep_charge`,`lep_id`,`lep_d0`,`lep_d0sig`,`lep_z0` and `lep_z0sig`
    # sorting them by lep_pt in descending order
    sorted_pt=ak.argsort(-akw_array.lep_pt)
    if(ak.count_nonzero(sorted_pt)>0): # Ignore if no entry to change
        # akw_array.lep_pt=akw_array.lep_pt[sorted_pt]
        # akw_array.lep_eta=akw_array.lep_eta[sorted_pt]
        # akw_array.lep_phi=akw_array.lep_phi[sorted_pt]
        # akw_array.lep_E=akw_array.lep_E[sorted_pt]
        # akw_array.lep_ptiso=akw_array.lep_ptiso[sorted_pt]
        # akw_array.lep_etiso=akw_array.lep_etiso[sorted_pt]
        # akw_array.lep_charge=akw_array.lep_charge[sorted_pt]
        # akw_array.lep_id=akw_array.lep_id[sorted_pt]
        # akw_array.lep_d0=akw_array.lep_d0[sorted_pt]
        # akw_array.lep_d0sig=akw_array.lep_d0sig[sorted_pt]
        # akw_array.lep_z0=akw_array.lep_z0[sorted_pt]
        # akw_array.lep_z0sig=akw_array.lep_z0sig[sorted_pt]
        akw_array["lep_pt"]=akw_array.lep_pt[sorted_pt]
        akw_array["lep_eta"]=akw_array.lep_eta[sorted_pt]
        akw_array["lep_phi"]=akw_array.lep_phi[sorted_pt]
        akw_array["lep_E"]=akw_array.lep_E[sorted_pt]
        akw_array["lep_ptiso"]=akw_array.lep_ptiso[sorted_pt]
        akw_array["lep_etiso"]=akw_array.lep_etiso[sorted_pt]
        akw_array["lep_charge"]=akw_array.lep_charge[sorted_pt]
        akw_array["lep_id"]=akw_array.lep_id[sorted_pt]
        akw_array["lep_d0"]=akw_array.lep_d0[sorted_pt]
        akw_array["lep_d0sig"]=akw_array.lep_d0sig[sorted_pt]
        akw_array["lep_z0"]=akw_array.lep_z0[sorted_pt]
        akw_array["lep_z0sig"]=akw_array.lep_z0sig[sorted_pt]
    return akw_array

def sort_jets_by_pt(akw_array):
    # To students: No need to modify this function!
    #
    # Function accepts akward array with fields 'jet_pt', 'jet_eta' and 'jet_phi' sorting them by jet_pt in descending order
    sorted_pt=ak.argsort(-akw_array.jet_pt)
    if(ak.count_nonzero(sorted_pt)>0): # Ignore if no entry to change
        akw_array["jet_pt"]=akw_array.jet_pt[sorted_pt]
        akw_array["jet_eta"]=akw_array.jet_eta[sorted_pt]
        akw_array["jet_phi"]=akw_array.jet_phi[sorted_pt]
    return akw_array

def preselection_cuts(akw):
    # To students: No need to modify this function!
    #
    # Trigger & overlap removal cuts
    
    
    # Necessary treatment to remove Zbb overlap in samples.
    # Do not change that! Just take it as it is.
    # If you are really interested in this, ask your tutor. 
    hfor = akw.hfor_type
    channel = akw.m_mc_channel_number
    IsOverlapping_part1=ak.all([ak.any([ak.all([channel>=147105,channel <=147110],axis=0),
                ak.all([channel>=178354,channel<=178358],axis=0),
                ak.all([channel>=178369,channel<=178373],axis=0)],axis=0),
            akw.isData==False,
            akw.hfor_type==0],axis=0)
    IsOverlapping_part2=ak.all([akw.isData==False,
                                ak.any([ak.all([channel>=147113,channel <=147118],axis=0),
                                        ak.all([channel>=178374,channel<=178378],axis=0),
                                        ak.all([channel>=178359,channel<=178363],axis=0)],axis=0),
                                akw.hfor_type==0],axis=0)
    IsOverlapping=ak.any([IsOverlapping_part1,IsOverlapping_part2],axis=0)
    akw=akw[ak.where(IsOverlapping,False,True)]
    
    #Apply Trigger
    PassTrigger=ak.any([akw.EF_e24vhi_medium1,akw.EF_e60_medium1,akw.EF_2e12Tvh_loose1, 
          akw.EF_mu24i_tight,akw.EF_mu36_tight,akw.EF_2mu13,
          akw.EF_mu18_tight_mu8_EFFS,akw.EF_e12Tvh_medium1_mu8,
          akw.EF_e24vhi_loose1_mu8],axis=0)
    akw=akw[ak.where(PassTrigger,True,False)]
    
    return akw

def blind_data(histogram_stack_ref,data_name='data',background_names=['ggH125', 'VBFH125', 'WH125', 'ZH125', 'ttH125', 'ZZ', 'WZ', 'Zee', 'Zeebb', 'Zmumu', 'Zmumubb', 'Ztautau', 'Top']):
    histogram_stack=histogram_stack_ref.copy()
    backgrounds = [histogram_stack[name] for name in background_names]
    hBkg=sum_histograms(backgrounds)
    
    histogram_stack[data_name]=hBkg.Clone()
    histogram_stack[data_name].SetTitle("Data (blind)")
    for binIndex in range(histogram_stack[data_name].GetNbinsX()+1):
        histogram_stack[data_name].SetBinContent(binIndex,0.)
        histogram_stack[data_name].SetBinError(binIndex,0.)
    return histogram_stack
