import os
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.colors as mcolors
import awkward as ak
import numpy as np
import ROOT as rt

class Plotter():

    def __init__(self):
        self.label_upper_right = r'$\it{MAIA}$ Detector Concept'
        self.data_label = 'Simulation, no BIB'
        self.lattice_label = r'Lattice v04'
        self.com_tev = 10
        self.outdir = None

    def SetUpperRightLabel(self,val):
        self.label_upper_right = val

    def SetDataLabel(self,val):
        self.data_label = val

    def SetLatticeLabel(self,val):
        self.lattice_label = val

    def SetCOMTev(self,val):
        self.com_tev = val

    def SetOutputDirectory(self,val):
        self.outdir = val
        if(self.outdir is not None):
            os.makedirs(self.outdir,exist_ok=True)

    @staticmethod
    def _theta_to_eta(theta):
        return -np.log(np.tan(theta / 2.0))

    @staticmethod
    def _eta_to_theta(eta):
        return 2.0 * np.arctan(np.exp(-eta))

    @staticmethod
    def _add_eta_axis(ax, eta_ticks=None, fontsize=18, position=-0.17, top=False):
        """Add a secondary x-axis showing eta corresponding to theta [rad].
        top=True places it above the plot; top=False (default) places it below.
        """
        if eta_ticks is None:
            eta_ticks = np.array([-2, -1, 0, 1, 2])
        if top:
            secax = ax.secondary_xaxis('top', functions=(Plotter._theta_to_eta, Plotter._eta_to_theta))
            secax.set_xlabel(r'Muon $\eta$', fontsize=fontsize, labelpad=6, loc='right')
        else:
            ax.xaxis.labelpad = 2
            secax = ax.secondary_xaxis(position, functions=(Plotter._theta_to_eta, Plotter._eta_to_theta))
            secax.set_xlabel('')
            # Use ax.transAxes so the label sits at exactly x=1.0 — same reference as the
            # main xlabel (loc='right'), guaranteeing horizontal alignment.
            ax.text(1.0, position - 0.05, r'Muon $\eta$',
                    transform=ax.transAxes, ha='right', va='top', fontsize=fontsize)
        secax.set_xticks(eta_ticks)
        secax.tick_params(which='major', length=12, width=1.6)
        secax.tick_params(which='minor', length=6,  width=1.0)
        secax.minorticks_on()
        return secax

    @staticmethod
    def _add_theta_axis(ax, theta_ticks=None, fontsize=18, position=-0.17, top=False):
        """Add a secondary x-axis showing theta [rad] when the main axis is eta.
        top=True places it above; top=False (default) places it below.
        """
        if top:
            secax = ax.secondary_xaxis('top', functions=(Plotter._eta_to_theta, Plotter._theta_to_eta))
            secax.set_xlabel(r'True Muon $\theta$ [rad]', fontsize=fontsize, loc='right')
            secax.tick_params(which='major', length=12, width=1.6, top=True)
            secax.tick_params(which='minor', length=6,  width=1.0, top=True)
            if theta_ticks is not None:
                secax.set_xticks(theta_ticks)
        else:
            ax.xaxis.labelpad = 2
            secax = ax.secondary_xaxis(position, functions=(Plotter._eta_to_theta, Plotter._theta_to_eta))
            secax.set_xlabel('')
            ax.text(1.0, position - 0.07, r'True Muon $\theta$ [rad]',
                    transform=ax.transAxes, ha='right', va='top', fontsize=fontsize)
            secax.tick_params(which='major', length=12, width=1.6)
            secax.minorticks_off()
            # Theta ticks at η = ±2, ±1.5, ±1, ±0.5, 0 (labeled, symmetric about π/2)
            if theta_ticks is None:
                eta_ref = np.array([-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.])
                theta_ticks = np.sort(2.0 * np.arctan(np.exp(-eta_ref)))
            secax.set_xticks(theta_ticks)
            secax.set_xticklabels([f'{t:.2f}' for t in theta_ticks])
            secax.tick_params(which='major', labelsize=fontsize - 6)
            # ±0.5, ±1.5 minor ticks on the main eta axis (tick mark only, no label)
            ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([-1.5, -0.5, 0.5, 1.5]))
            ax.tick_params(axis='x', which='minor', length=6, width=1.0)
        return secax

    def plot_processed_data(self,processed_results, labels=None, xlabel='', ylabel='', title='', misctext=None, fontsize=20, log=False, xlog=False, ylim=None, xlim=None, legend_loc=(0.4,0.71), label_block_y_up=0.92, misctext_y_up=0.77, savename='', eta_axis=False, theta_axis=False):
        """
        Plot the processed RMS data.

        Parameters:
            processed_results (list): List of dictionaries containing processed data.
            labels (list of str, optional): Labels for the data series.
            xlabel (str, optional): X-axis label.
            ylabel (str, optional): Y-axis label.
            title (str, optional): Plot title.
            fontsize (int, optional): Font size for the plot labels.
            log (bool, optional): Use a logarithmic y-axis scale.
            ylim (tuple, optional): Y-axis limits.
            xlim (tuple, optional): X-axis limits.
            savepdf (bool, optional): Save the plot as a PDF file.
            savename (str, optional): Name of the saved PDF file.
        """
        fig, ax = plt.subplots(figsize=(9,6))

        # Loop over the processed results (dictionary of bin centers, rms values, sem values, and theta width) and create errorbar graphs
        for i, hist in enumerate(processed_results):
            label = labels[i] if labels else ''
            color = None if labels else plt.rcParams['axes.prop_cycle'].by_key()['color'][0]
            nbins = hist.GetNbinsX()
            print("i:", i)
            print("nbins:", nbins)
            bin_centers = [hist.GetBinCenter(i+1) for i in range(nbins)]
            rms_values = [hist.GetBinContent(i+1) for i in range(nbins)]
            errors = [hist.GetBinError(i+1) for i in range(nbins)]
            bin_widths = np.array([hist.GetBinWidth(i+1) for i in range(nbins)])
            ax.errorbar(
                bin_centers,
                rms_values,
                xerr=bin_widths/2.,
                yerr=errors,
                fmt='o',
                markersize=3, #6
                label=label,
                markerfacecolor=color,
                markeredgecolor=color,
                ecolor=color
            )

        # Custom Muon Collider text
        if len(title) > 0:
            ax.set_title(title, fontsize=fontsize)
        else:
            ax.set_title(self.label_upper_right, fontsize=fontsize, loc = 'right')

        plt.text(0.04, label_block_y_up, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
        plt.text(0.04, label_block_y_up-0.075, self.data_label, transform=plt.gca().transAxes)
        com_label = r'$\sqrt{s}$ = ' +r'{}'.format(self.com_tev) +   r' TeV'
        combined_label = '{}, {}'.format(self.lattice_label,com_label)
        plt.text(0.04, label_block_y_up-0.15, combined_label, transform=plt.gca().transAxes)
        # plt.text(0.04, 0.71, com_label, transform=plt.gca().transAxes)

        if(misctext is not None):
            if(type(misctext) != list):
                misctext = [misctext]
            for i,line in enumerate(misctext):
                plt.text(0.04, misctext_y_up - 0.075 * (i+1), line, fontsize=fontsize-3, transform=plt.gca().transAxes)

        # handle axes
        ax.set_xlabel(xlabel, loc='right', fontsize=fontsize+5)
        ax.set_ylabel(ylabel, loc='top'  , fontsize=fontsize+5)
        if log : ax.set_yscale('log')
        if xlog: ax.set_xscale('log')
        if ylim is not None: ax.set_ylim(ylim)
        if xlim is not None: ax.set_xlim(xlim)
        if(labels):
            ax.legend(fontsize=fontsize-3,loc=legend_loc)
        ax.tick_params(labelsize=fontsize)
        #ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        if eta_axis:
            self._add_eta_axis(ax, fontsize=fontsize+5)
        if theta_axis:
            self._add_theta_axis(ax, fontsize=fontsize+5)

        if len(savename)>0:
            if(self.outdir is not None):
                savename = '{}/{}'.format(self.outdir,savename)
            plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

        plt.show()

    def plot_efficiencies(self,results, min_value, max_value, xlabel=None, labels="", misctext='', bottom_label=None, savename='',xlim=None,ylim=None, label_block_y_up=None, upTo10GeVBool=False, eta_axis=False, theta_axis=False):
        plt.figure(figsize=(8, 6))

        for i, hist in enumerate(results):
            # depending on version of code, might be handling TH1 or TEfficiency
            if(type(hist) in [rt.TH1, rt.TH1D, rt.TH1F]):
                nbins = hist.GetNbinsX()
                print(f"\nHistogram {i} (TH1):")
                for ibin in range(1, nbins + 1):
                    n_entries = hist.GetBinContent(ibin)
                    low_edge  = hist.GetBinLowEdge(ibin)
                    high_edge = low_edge + hist.GetBinWidth(ibin)

                    print(f"  Bin {ibin}: [{low_edge:.3f}, {high_edge:.3f})  entries = {n_entries}")
                bin_centers = [hist.GetBinCenter(i+1) for i in range(nbins)]
                efficiencies = [hist.GetBinContent(i+1) for i in range(nbins)]
                errors_up = [hist.GetBinError(i+1)/2. for i in range(nbins)]
                errors_down = [hist.GetBinError(i+1)/2. for i in range(nbins)]
                bin_widths = np.array([hist.GetBinWidth(i+1) for i in range(nbins)])
            else:
                h_total = hist.GetTotalHistogram()
                nbins = h_total.GetNbinsX()
                h_pass  = hist.GetPassedHistogram()

                print(f"\nHistogram {i} (TEfficiency):")
                for ibin in range(1, nbins + 1):
                    n_total = h_total.GetBinContent(ibin)
                    n_pass  = h_pass.GetBinContent(ibin)
                    eff     = hist.GetEfficiency(ibin)

                    low_edge  = h_total.GetBinLowEdge(ibin)
                    high_edge = low_edge + h_total.GetBinWidth(ibin)

                    print(
                        f"  Bin {ibin}: [{low_edge:.3f}, {high_edge:.3f})  "
                        f"passed = {n_pass}, total = {n_total}, eff = {eff:.3f}"
                    )
                bin_centers = [h_total.GetBinCenter(i+1) for i in range(nbins)]
                efficiencies = [hist.GetEfficiency(i+1) for i in range(nbins)]
                errors_up = [hist.GetEfficiencyErrorUp(i+1) for i in range(nbins)]
                errors_down = [hist.GetEfficiencyErrorLow(i+1) for i in range(nbins)]
                bin_widths = np.array([h_total.GetBinWidth(i+1) for i in range(nbins)])

            label = labels[i] if labels else None
            plt.errorbar(
                bin_centers,
                efficiencies,
                yerr=(errors_down,errors_up),
                xerr=bin_widths/2.,
                fmt='o',
                label=label,
                markersize=3
            )


        # Set x and y limits (with some leeway)

        #if "bib_eff_vs_pt" in savename: max_value=1000
        #if "nobib_eff_vs_pt_barrel" in savename: max_value=5000
        #if "nobib_eff_vs_pt_endcap" in savename: max_value=3000
        # if "vs_theta" in savename:
        #     minvalue = 0
        #     maxvalue = 180
        if(xlim is None):
            if "vs_pt" in savename :
                if(not upTo10GeVBool):
                    plt.xscale('log')
                
                x, y = [min_value, max_value], [1, 1]
                plt.xlim(min_value, max_value)
            else :
                plt.xlim(min_value-10, max_value+10)
                x, y = [min_value-10, max_value+10], [1, 1]
        else:
            x, y = xlim, [1, 1]
            plt.xlim(*xlim)
        plt.xlim(min_value, max_value)

        if(ylim is None): ylim = (0,1.35)
        plt.ylim(*ylim)

        plt.plot(x, y, linestyle="dashed",linewidth=1,color='xkcd:grey')

        if xlabel is not None: plt.xlabel(xlabel, loc='right', fontsize=20)
        plt.ylabel('Reconstruction Efficiency', loc='top', fontsize=20)
        plt.title(self.label_upper_right, fontsize=20, loc = 'right')

        # Custom Muon Collider text
        if(label_block_y_up is None):
            label_block_y_up = 0.92
        plt.text(0.04, label_block_y_up, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
        plt.text(0.04, label_block_y_up-0.075, self.data_label, transform=plt.gca().transAxes)
        combined_label = self.lattice_label + r', ' + r'$\sqrt{s}$ = ' +r'{}'.format(self.com_tev) +   r' TeV'
        plt.text(0.04, label_block_y_up-0.15 , combined_label, transform=plt.gca().transAxes)

        # bottom label
        if(bottom_label is not None):
            if(type(bottom_label) != list): bottom_label = [bottom_label]
            for i,line in enumerate(bottom_label):
                plt.text(0.42, 0.14 - i*0.07, line, transform=plt.gca().transAxes)

        if len(misctext): plt.text(0.76-len(misctext)*0.004, 0.92, misctext, transform=plt.gca().transAxes)

        plt.xticks(fontsize=20)
        if "nobib_eff_vs_pt_barrel" in savename:
            plt.xticks(ticks=[1,10,100,1000,5000],fontsize=20)
        plt.yticks(fontsize=20)
        # plt.legend(loc='lower left', fontsize=20)
        if(labels is not None):
            plt.legend(loc=(0.575, 0.725), fontsize=16)

        if eta_axis:
            self._add_eta_axis(plt.gca(), fontsize=20)
        if theta_axis:
            self._add_theta_axis(plt.gca(), fontsize=20)

        if len(savename)>0:
            if(self.outdir is not None):
                savename = '{}/{}'.format(self.outdir,savename)
            plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

        plt.show()

    # Histogram plotting function that is intended to compare fake and truth-matched data
    #def plot_fake_and_true_distributions(self,data, fake_key, truth_key, bins, x_label, y_label, x_range=None, x_scale ='linear', y_scale='linear', savename='', custom_data_func=None):
    def plot_fake_and_true_distributions(self,data, fake_key, truth_key, bins, x_label, y_label, x_range=None, y_scale='linear', savename='', custom_data_func=None):

        # Set up plotting parameters
        plt.style.use('seaborn-v0_8-colorblind')
        fontsize = 20
        plt.rcParams['font.size'] = fontsize

        plt.figure(figsize=(8, 6))

        if custom_data_func:
            fake_data, truth_data = custom_data_func(data)
        # Data in these jsons is stored as jagged arrays varying per event per track etc. so we need to flatten them
        else:
            fake_data = ak.flatten(data[fake_key])
            truth_data = ak.flatten(data[truth_key])

        # Bins input can either be a single itn or can be used as a np.linspace array
        if isinstance(bins, int):
            bins = (bins, bins)
        elif isinstance(bins, np.ndarray):
            bins = (bins, bins)

        # Create the histograms
        fake_weights  = np.ones_like(fake_data) / len(fake_data)
        truth_weights = np.ones_like(truth_data) / len(truth_data)
        plt.hist(fake_data, bins=bins[0], linewidth=1.5, histtype='step', weights=fake_weights, color='red', label='Fake')#density=True,
        plt.hist(truth_data, bins=bins[1], linewidth=1.5, histtype='step', weights=truth_weights, color='blue', label='Truth-Matched')

        if x_range:
            plt.xlim(x_range)

        plt.xlabel(x_label, loc='right')
        plt.ylabel(y_label, loc='top')
        #plt.xscale(x_scale)
        plt.yscale(y_scale)
        plt.title(r'$\it{MAIA}$ Detector Concept', fontsize=fontsize, loc='right')

        # Get current y-axis limits and add 10% padding to the top # TODO: This adds 40% for linear, no? Seems fine. -Jan
        y_min, y_max = plt.gca().get_ylim()
        plt.ylim(y_min, y_max * 1.4)
        if y_scale=="log": plt.ylim(y_min, y_max * 10)
        # ---- Adjust x-limits ----
        #x_min, x_max = plt.gca().get_xlim()
        #if x_scale == "log":
        #    plt.gca().set_xlim(0, 10e6)

        # Custom Muon Collider text
        plt.text(0.04, 0.9, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
        plt.text(0.04, 0.83, self.data_label, transform=plt.gca().transAxes)
        plt.text(0.04, 0.76, self.lattice_label, transform=plt.gca().transAxes)
        com_label = r'$\sqrt{s}$ = ' +r'{}'.format(self.com_tev) +   r' TeV'
        plt.text(0.04, 0.69, com_label, transform=plt.gca().transAxes)

        # Or use mplhep to add the text
        # hep.atlas.text("Muon Collider")

        plt.legend(frameon=False, loc = 'upper right', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.tight_layout()

        dist=fake_key.removeprefix("fake_")

        if(savename == ''):
            savename="{}_dist".format(dist)

        if len(savename)>0:
            if(self.outdir is not None):
                savename = '{}/{}'.format(self.outdir,savename)
            plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

        plt.show()

    def plot_distributions(self,data, labels, bins, x_label, y_label, x_range=None, y_scale='linear', savename=''):

        # Set up plotting parameters
        plt.style.use('seaborn-v0_8-colorblind')
        fontsize = 20
        plt.rcParams['font.size'] = fontsize

        plt.figure(figsize=(8, 6))

        # Bins input can either be a single itn or can be used as a np.linspace array
        if isinstance(bins, int):
            bins = (bins, bins)
        elif isinstance(bins, np.ndarray):
            bins = (bins, bins)

        # Create the histograms
        for i,array in enumerate(data):
            plt.hist(ak.flatten(array), bins=bins[0], linewidth=1.5, histtype='step', label=labels[i])

        if x_range:
            plt.xlim(x_range)

        plt.xlabel(x_label, loc='right')
        plt.ylabel(y_label, loc='top')
        plt.yscale(y_scale)
        plt.title(r'$\it{MAIA}$ Detector Concept', fontsize=fontsize, loc='right')

        # Get current y-axis limits and add 10% padding to the top # TODO: This adds 40% for linear, no? Seems fine. -Jan
        y_min, y_max = plt.gca().get_ylim()
        plt.ylim(y_min, y_max * 1.4)
        if y_scale=="log": plt.ylim(y_min, y_max * 10)

        # Custom Muon Collider text
        plt.text(0.04, 0.9, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
        plt.text(0.04, 0.83, self.data_label, transform=plt.gca().transAxes)
        plt.text(0.04, 0.76, self.lattice_label, transform=plt.gca().transAxes)
        com_label = r'$\sqrt{s}$ = ' +r'{}'.format(self.com_tev) +   r' TeV'
        plt.text(0.04, 0.69, com_label, transform=plt.gca().transAxes)

        # Or use mplhep to add the text
        # hep.atlas.text("Muon Collider")

        plt.legend(frameon=False, loc = 'upper right', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.tight_layout()

        if len(savename)>0:
            if(self.outdir is not None):
                savename = '{}/{}'.format(self.outdir,savename)
            plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

        plt.show()

    @staticmethod
    def _make_symlog_bins(lo, hi, n, linthresh):
        """
        Create bin edges matched to a symlog scale so that each displayed bin
        has roughly equal visual width.  Bins are log-spaced outside linthresh
        and linear inside it.
        """
        n_log  = max(4, int(n * 0.45))   # bins per log decade side
        n_lin  = max(3, n - 2 * n_log)   # bins in the linear region

        neg_bins = -np.geomspace(linthresh, abs(lo), n_log + 1)[::-1]   # lo … -linthresh
        lin_bins =  np.linspace(-linthresh, linthresh, n_lin + 1)        # -linthresh … +linthresh
        pos_bins =  np.geomspace(linthresh, hi, n_log + 1)               # linthresh … hi

        return np.unique(np.concatenate([neg_bins, lin_bins[1:-1], pos_bins]))

    def plot_2d_histogram(self, datax, datay, bins=100, xlabel='', ylabel='',
                          xlim=None, ylim=None, savename='', misctext=None,
                          log_color=True, log_y=False, square_y=False, fontsize=20):
        """
        Plot a 2D histogram (TH2F equivalent) as a heatmap.

        Parameters:
            datax: x-axis data (awkward or numpy array, will be flattened).
            datay: y-axis data (awkward or numpy array, will be flattened).
            bins (int or [int,int]): Number of bins in x and y.
            xlabel (str): X-axis label.
            ylabel (str): Y-axis label.
            xlim (tuple): X-axis range, also used to clip data.
            ylim (tuple): Y-axis range, also used to clip data.
            savename (str): Output filename stem (no extension).
            misctext (str or list): Extra annotation lines.
            log_color (bool): Use log color scale.
            log_y (bool): Use symlog scale on the y-axis.  Bins are automatically
                          log-spaced to match the visual scale.
            square_y (bool): Square the y values before plotting (e.g. to show
                             residual² rather than residual).  Implies log_y=True
                             with a plain log scale since the result is positive.
            fontsize (int): Font size for labels.
        """
        x = ak.to_numpy(ak.flatten(datax)).astype(float)
        y = ak.to_numpy(ak.flatten(datay)).astype(float)

        if square_y:
            y = y ** 2
            log_y = True   # always log when squared (all positive)

        # Remove NaN/Inf entries
        finite_mask = np.isfinite(x) & np.isfinite(y)
        x = x[finite_mask]
        y = y[finite_mask]

        fig, ax = plt.subplots(figsize=(10, 7))

        # --- build bin arrays -------------------------------------------------
        n_x = bins[0] if isinstance(bins, (list, tuple)) else bins
        n_y = bins[1] if isinstance(bins, (list, tuple)) else bins

        x_lo = float(xlim[0]) if xlim is not None else float(x.min())
        x_hi = float(xlim[1]) if xlim is not None else float(x.max())
        y_lo = float(ylim[0]) if ylim is not None else float(y.min())
        y_hi = float(ylim[1]) if ylim is not None else float(y.max())

        xbins = np.linspace(x_lo, x_hi, n_x + 1)

        if square_y:
            # All positive after squaring → plain log bins
            y_lo_pos = max(y_lo, float(np.percentile(y[y > 0], 0.1))) if np.any(y > 0) else 1e-10
            ybins = np.geomspace(y_lo_pos, y_hi, n_y + 1)
        elif log_y:
            linthresh = float(np.percentile(np.abs(y[y != 0]), 1)) if np.any(y != 0) else 1e-3
            ybins = self._make_symlog_bins(y_lo, y_hi, n_y, linthresh)
        else:
            ybins = np.linspace(y_lo, y_hi, n_y + 1)

        # Clip data to the chosen range so out-of-range entries don't appear
        in_range = (x >= x_lo) & (x <= x_hi) & (y >= ybins[0]) & (y <= ybins[-1])
        x = x[in_range]
        y = y[in_range]

        # --- draw -------------------------------------------------------------
        weights = np.ones_like(x) / len(x)
        norm = mcolors.LogNorm() if log_color else None
        _, _, _, img = ax.hist2d(x, y, bins=[xbins, ybins], weights=weights, norm=norm)
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label('Normalized tracks', fontsize=fontsize - 2)
        cbar.ax.tick_params(labelsize=fontsize - 4)

        if square_y:
            ax.set_yscale('log')
        elif log_y:
            ax.set_yscale('symlog', linthresh=linthresh)

        # --- labels -----------------------------------------------------------
        ax.set_title(self.label_upper_right, fontsize=fontsize, loc='right')
        ax.set_xlabel(xlabel, loc='right', fontsize=fontsize + 3)
        ax.set_ylabel(ylabel, loc='top',   fontsize=fontsize + 3)
        ax.tick_params(labelsize=fontsize - 2)

        label_block_y_up = 0.55
        ax.text(0.04, label_block_y_up,        "Muon Collider", fontweight='bold', style='italic', transform=ax.transAxes, fontsize=fontsize)
        ax.text(0.04, label_block_y_up - 0.075, self.data_label, transform=ax.transAxes, fontsize=fontsize - 2)
        com_label = r'$\sqrt{s}$ = ' + r'{}'.format(self.com_tev) + r' TeV'
        combined_label = '{}, {}'.format(self.lattice_label, com_label)
        ax.text(0.04, label_block_y_up - 0.15, combined_label, transform=ax.transAxes, fontsize=fontsize - 2)

        if misctext is not None:
            if not isinstance(misctext, list):
                misctext = [misctext]
            for i, line in enumerate(misctext):
                ax.text(0.04, label_block_y_up - 0.225 - 0.075 * i, line, fontsize=fontsize - 4, transform=ax.transAxes)

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(ybins[0], ybins[-1])

        if len(savename) > 0:
            out = '{}/{}'.format(self.outdir, savename) if self.outdir is not None else savename
            plt.savefig(out + ".pdf", format='pdf', bbox_inches='tight')

        plt.show()