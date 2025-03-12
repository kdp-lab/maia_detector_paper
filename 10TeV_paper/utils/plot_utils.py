
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np

# Histogram plotting function that is intended to compare fake and truth-matched data
def PlotHistogram(data, fake_key, truth_key, bins, x_label, y_label, x_range=None, y_scale='linear', file_name=None, custom_data_func=None):

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
    plt.text(0.04, 0.83, "Simulation, with BIB", transform=plt.gca().transAxes)
    plt.text(0.04, 0.76, r'Lattice v0.4', transform=plt.gca().transAxes)
    plt.text(0.04, 0.69, r'$\sqrt{s}$ = 10 TeV', transform=plt.gca().transAxes)

    # Or use mplhep to add the text
    # hep.atlas.text("Muon Collider")

    plt.legend(frameon=False, loc = 'upper right', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()

    dist=fake_key.removeprefix("fake_")
    file_name="bib_{}_dist.pdf".format(dist)
    plt.savefig(file_name)
    plt.show()


# In[5]:


# Function to choose which histogram to plot based on the plot_type argument
def choose_histogram(data, plot_type):

    if plot_type == 'pt':
        PlotHistogram(data, 'fake_pt', 'mcp_mu_pt', np.linspace(0, 1000, 100), 'Track $p_T$ [GeV]', 'Normalized Count', x_range=(0, 1000), y_scale='log')

    elif plot_type == 'phi':
        PlotHistogram(data, 'fake_phi', 'mcp_mu_phi', (100), r'Track $\phi$', 'Normalized Count', x_range=(-np.pi, np.pi), y_scale='linear')

    elif plot_type == 'eta':
        PlotHistogram(data, 'fake_eta', 'mcp_mu_eta', (20), r'Track $\eta$', 'Normalized Count', y_scale='linear')

    elif plot_type == 'chi2_ndf':
        PlotHistogram(data, 'fake_chi2', 'LC_chi2', (30), r'Track $\chi^2/n_{dof}$', 'Normalized Count', x_range=(0, 3), y_scale='linear', custom_data_func=lambda d: (ak.flatten(d['fake_chi2']) / ak.flatten(d['fake_ndf']), ak.flatten(d['LC_chi2']) / ak.flatten(d['LC_ndf'])))

    elif plot_type == 'd0':
        PlotHistogram(data, 'fake_d0', 'LC_d0', (50), r'Track $d_0$ [mm]', 'Normalized Count', x_range=(-6,6), y_scale='log')

    elif plot_type == 'nhits':
        PlotHistogram(data, 'fake_nhits', 'LC_nhits', np.arange(-0.5, 26, 1), r'Track $n_{hits}$', 'Normalized Count', y_scale='linear', custom_data_func=lambda d: (ak.flatten(d['fake_nhits']), ak.flatten(d['LC_nhits'])))


def plot_efficiencies(results, min_value, max_value, xlabel=None, ylim=None, bib=False, labels="", misctext='', savepdf=False, savename=''):
    plt.figure(figsize=(8, 6))

    # Assign plotting parameters
    alpha = 1.0
    size = 10

    # When plotting multiple datasets, loop through label names, decrease the size and alpha of the markers
    for i, (bin_centers, efficiencies, errors, bin_widths) in enumerate(results):
        label = labels[i] if labels else f"Dataset {i+1}"
        plt.errorbar(bin_centers, efficiencies, yerr=errors, xerr=bin_widths, fmt='o', label=label, markersize=size, alpha=alpha)
        size -= 2
        alpha -= 0.2

    # Set x and y limits (with some leeway)

    #if "bib_eff_vs_pt" in savename: max_value=1000
    #if "nobib_eff_vs_pt_barrel" in savename: max_value=5000
    #if "nobib_eff_vs_pt_endcap" in savename: max_value=3000
    if "vs_theta" in savename:
        minvalue = 0
        maxvalue = 180
    if "vs_pt" in savename :
        plt.xscale('log')
        x, y = [min_value, max_value], [1, 1]
        plt.xlim(min_value, max_value)
    else :
        plt.xlim(min_value-10, max_value+10)
        x, y = [min_value-10, max_value+10], [1, 1]
    plt.ylim(0, 1.35)


    plt.plot(x, y, linestyle="dashed")

    if xlabel is not None: plt.xlabel(xlabel, loc='right', fontsize=20)
    plt.ylabel('Reconstruction Efficiency', loc='top', fontsize=20)
    plt.title(r'$\it{MAIA}$ Detector Concept', fontsize=20, loc = 'right')

    # Custom Muon Collider text
    plt.text(0.04, 0.92, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
    label = "Simulation"
    if bib: label += ', with BIB'
    else  : label += ', no BIB'
    plt.text(0.04, 0.85, label, transform=plt.gca().transAxes)
    plt.text(0.04, 0.77, r'Lattice v04, $\sqrt{s}$ = 10 TeV', transform=plt.gca().transAxes)
    if len(misctext): plt.text(0.8-len(misctext)*0.004, 0.92, misctext, transform=plt.gca().transAxes)



    plt.xticks(fontsize=20)
    if "nobib_eff_vs_pt_barrel" in savename:
        plt.xticks(ticks=[1,10,100,1000,5000],fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(loc='lower left', fontsize=20)

    if len(savename)>0:
        plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')


    plt.show()
