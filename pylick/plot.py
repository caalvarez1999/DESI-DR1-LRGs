import matplotlib.pyplot as plt 
import matplotlib as mpl 
import numpy as np

from pylick._config import __style__, __dirPlot__




def _update_rcParams():
    mpl.rcParams["lines.linewidth"] = 1.
    mpl.rcParams["axes.linewidth"] = 1.
    mpl.rcParams["axes.labelsize"] = 16.
    mpl.rcParams["xtick.top"] = True
    mpl.rcParams["xtick.labelsize"] = 14
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['xtick.major.width'] = 1.
    mpl.rcParams["ytick.right"] = True
    mpl.rcParams["ytick.labelsize"] = 14
    mpl.rcParams["ytick.direction"] = "in"
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['ytick.major.width'] = 1.

    try:
        plt.style.use(__style__)
    except:
        pass
    




def plot_spectrum(ID, x, y, yerr, mask, regions, latex_names, units, done, settings={}):

    if not 'figsize' in settings:
        settings['figsize'] = (8,4.8)

    if not 'xlab' in settings:
        settings['xlab'] = r'Restframe wavelength [$\AA$]'
    
    if not 'ylab' in settings:
        settings['ylab'] = r'Flux [unitless]'

    if not 'title' in settings:
        settings['title'] = ID

    if not 'spec_colors' in settings:
        settings['spec_cols'] = 'navy', 'deepskyblue', (.4,.4,.4,1), (.4,.4,.4,1)

    if not 'inspect' in settings:
        settings['inspect'] = False
    else:
        settings['spec_cols'] = 'navy', 'deepskyblue', 'darkgoldenrod'


    if not 'outpath' in settings:
        settings['outpath'] = __dirPlot__ 

    if not 'format' in settings:
        settings['format'] = '.pdf'
        

    _update_rcParams()

    fig, ax = plt.subplots(figsize=settings['figsize'])

    ax.step(x, y, where='mid', color=settings['spec_cols'][0])
    ax.fill_between(x, y-yerr, y+yerr, step='mid',color=settings['spec_cols'][1])

    for i in range(np.size(done)):
        if (units[i] == 'A' or units[i] == 'mag'):
            l, r = regions[i][1][0],regions[i][1][1]

            if done[i]:
                ax.axvspan(l, r , alpha=0.3, facecolor=settings['spec_cols'][2], edgecolor=settings['spec_cols'][3], lw=1, zorder=0)
                ax.annotate(r"${:s}$".format(latex_names[i]), (0.5*(l+r), 0.05), 
                            xycoords=("data", "axes fraction"), rotation=90, ha='center',va='bottom', fontsize=14)
            
            else:
                if settings['inspect']:
                    ax.axvspan(l, r , alpha=0.3, color='darkred', zorder=0)
                    ax.annotate(r"${:s}$".format(latex_names[i]), (0.5*(l+r), 0.05), 
                            xycoords=("data", "axes fraction"), rotation=90, ha='center',va='bottom', fontsize=14)

            

    if settings['inspect'] is True:
        markbad = x[mask.astype(bool)]
        for llambda in markbad:
            ax.axvline(llambda, color='lightgray', lw=1, zorder=0)

    plt.title(settings['title'])
    plt.xlabel(settings['xlab'])
    plt.ylabel(settings['ylab'])
    # plt.xlim(x.min()+3, x.max()-3)


    plt.subplots_adjust(left=0.125,bottom=0.135,right=0.95,top=0.92)

    if settings['inspect']:
        plt.show()
    else:
        plt.savefig(settings['outpath']+ID+settings['format'])
        plt.close()
