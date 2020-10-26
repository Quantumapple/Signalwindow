import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import pandas as pd

def func(x, a):
    return a/x

def func2(x, a, b):
    return a/(x+b)

names = [
    'dphiL123', 'dphiL124', 'dphiL134', 'dphiL234', 'dphiL12D1', 'dphiL13D1',
    'dphiL23D1', 'dphiL12D2', 'dphiL1D12', 'dphiL2D12', 'dphiL1D13', 'dphiL1D23',
    'dphiD123', 'dphiD124', 'dphiD134', 'dphiD234', 'dphiD235', 'dphiD245', 'dphiD345'
]

pp = PdfPages('/home/jongho/fitting-scatter.pdf')
for inum, iname in enumerate(names):
    #if not inum == 0: continue
    txt = np.loadtxt('points/'+iname+'.txt')

    xdata = np.asarray(txt[0])
    ydata = np.asarray(txt[1])
    print(xdata.size)

    popt, pocv = curve_fit(func, xdata, ydata)
    popt2, pocv2 = curve_fit(func2, xdata, ydata)
    print(iname, '{:.4f}'.format(popt[0]))
    print(iname, '{:.4f}, {:.4f}'.format(*popt2))

    ### Chi-square test in python
    res = stats.chisquare(f_exp=ydata, f_obs=func(xdata, *popt))
    res2 = stats.chisquare(f_exp=ydata, f_obs=func2(xdata, *popt2))
    print('chi-square value: {:.5f}, p-value {:.5f}'.format(res[0], res[1]))
    print('chi-square value: {:.5f}, p-value {:.5f}'.format(res2[0], res2[1]))
    print()

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.grid(which='major')
    ax.grid(which='minor')
    ax.scatter(xdata, ydata)
    ax.plot(xdata, func(xdata, *popt), 'r-', label='fit: a=%.4f' % tuple(popt))
    ax.plot(xdata, func2(xdata, *popt2), 'g-', label='fit: a=%.4f, b=%.4f' % tuple(popt2))
    ax.set_xlabel('$\Delta\phi$', fontsize=15)
    ax.set_ylabel('Gen-level electron $p_T$', fontsize=15)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_xlim(0.0005, 0.2)
    #ax.set_xlim(None, 0.2)
    #ax.set_ylim(None, 10)
    ax.legend(fontsize=15)

    #ax.text(0.001, 2.5, r'$p_T = \frac{%.3f}{\Delta\phi}$' % tuple(popt), style='italic', fontsize=30, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})
    ax.text(0.006, 0.5, 'combination: %s' %(iname.split('dphi')[1]), style='italic', fontsize=15, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})
    #ax.text(0.006, 1, r'$\chi^2/ndf  (100 GeV) : {:1.4f}$'.format(res[0]/(xdata.size-1)), style='italic', fontsize=12, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})
    #ax.text(0.006, 0.75, r'$\chi^2/ndf  (10 GeV) : {:1.4f}$'.format(res2[0]/(xdata.size-1)), style='italic', fontsize=12, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})

    pp.savefig(fig)
    plt.close(fig)

pp.close()
