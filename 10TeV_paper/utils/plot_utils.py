import os
import matplotlib.pyplot as plt
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

    def plot_processed_data(self,processed_results, labels=None, xlabel='', ylabel='', title='', misctext=None, fontsize=20, log=False, xlog=False, ylim=None, xlim=None, savename=''):
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
        fig, ax = plt.subplots(figsize=(8,6))

        # Loop over the processed results (dictionary of bin centers, rms values, sem values, and theta width) and create errorbar graphs
        for i, hist in enumerate(processed_results):
            label = labels[i] if labels else ''
            color = None if labels else plt.rcParams['axes.prop_cycle'].by_key()['color'][0]
            nbins = hist.GetNbinsX()
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

        plt.text(0.04, 0.92, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
        plt.text(0.04, 0.85, self.data_label, transform=plt.gca().transAxes)
        plt.text(0.04, 0.78, self.lattice_label, transform=plt.gca().transAxes)
        com_label = r'$\sqrt{s}$ = ' +r'{}'.format(self.com_tev) +   r' TeV'
        plt.text(0.04, 0.71, com_label, transform=plt.gca().transAxes)

        if(misctext is not None):
            if(type(misctext) != list):
                misctext = [misctext]
            for i,line in enumerate(misctext):
                plt.text(0.04, 0.71 - 0.07 * (i+1), line, fontsize=fontsize-3, transform=plt.gca().transAxes)

        # handle axes
        ax.set_xlabel(xlabel, loc='right', fontsize=fontsize+5)
        ax.set_ylabel(ylabel, loc='top'  , fontsize=fontsize+5)
        if log : ax.set_yscale('log')
        if xlog: ax.set_xscale('log')
        if ylim is not None: ax.set_ylim(ylim)
        if xlim is not None: ax.set_xlim(xlim)
        if labels:
            ax.legend(fontsize=fontsize-3,loc=(0.5,0.71))
        ax.tick_params(labelsize=fontsize)
        #ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        if len(savename)>0:
            if(self.outdir is not None):
                savename = '{}/{}'.format(self.outdir,savename)
            plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

        plt.show()

    def plot_efficiencies(self,results, min_value, max_value, xlabel=None, labels="", misctext='', bottom_label=None, savename='',xlim=None):
        plt.figure(figsize=(8, 6))

        for i, hist in enumerate(results):
            nbins = hist.GetNbinsX()
            bin_centers = [hist.GetBinCenter(i+1) for i in range(nbins)]
            efficiencies = [hist.GetBinContent(i+1) for i in range(nbins)]
            errors = [hist.GetBinError(i+1) for i in range(nbins)]
            bin_widths = np.array([hist.GetBinWidth(i+1) for i in range(nbins)])

            label = labels[i] if labels else None
            plt.errorbar(
                bin_centers,
                efficiencies,
                yerr=errors,
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
                plt.xscale('log')
                x, y = [min_value, max_value], [1, 1]
                plt.xlim(min_value, max_value)
            else :
                plt.xlim(min_value-10, max_value+10)
                x, y = [min_value-10, max_value+10], [1, 1]
        else:
            x, y = xlim, [1, 1]
            plt.xlim(*xlim)

        plt.ylim(0, 1.35)

        plt.plot(x, y, linestyle="dashed")

        if xlabel is not None: plt.xlabel(xlabel, loc='right', fontsize=20)
        plt.ylabel('Reconstruction Efficiency', loc='top', fontsize=20)
        plt.title(self.label_upper_right, fontsize=20, loc = 'right')

        # Custom Muon Collider text
        plt.text(0.04, 0.92, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
        plt.text(0.04, 0.85, self.data_label, transform=plt.gca().transAxes)
        combined_label = self.lattice_label + r', ' + r'$\sqrt{s}$ = ' +r'{}'.format(self.com_tev) +   r' TeV'
        plt.text(0.04, 0.77, combined_label, transform=plt.gca().transAxes)

        # bottom label
        if(bottom_label is not None):
            if(type(bottom_label) != list): bottom_label = [bottom_label]
            for i,line in enumerate(bottom_label):
                plt.text(0.04, 0.14 - i*0.07, line, transform=plt.gca().transAxes)

        if len(misctext): plt.text(0.8-len(misctext)*0.004, 0.92, misctext, transform=plt.gca().transAxes)

        plt.xticks(fontsize=20)
        if "nobib_eff_vs_pt_barrel" in savename:
            plt.xticks(ticks=[1,10,100,1000,5000],fontsize=20)
        plt.yticks(fontsize=20)
        # plt.legend(loc='lower left', fontsize=20)
        if(labels is not None):
            plt.legend(loc=(0.575, 0.725), fontsize=16)

        if len(savename)>0:
            if(self.outdir is not None):
                savename = '{}/{}'.format(self.outdir,savename)
            plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

        plt.show()

    # Histogram plotting function that is intended to compare fake and truth-matched data
    def PlotHistogram(self,data, fake_key, truth_key, bins, x_label, y_label, x_range=None, y_scale='linear', savename='', custom_data_func=None):

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
        plt.yscale(y_scale)
        plt.title(r'$\it{MAIA}$ Detector Concept', fontsize=fontsize, loc='right')

        # Get current y-axis limits and add 10% padding to the top
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

        dist=fake_key.removeprefix("fake_")

        if(savename == ''):
            savename="{}_dist".format(dist)

        if len(savename)>0:
            if(self.outdir is not None):
                savename = '{}/{}'.format(self.outdir,savename)
            plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

        plt.show()