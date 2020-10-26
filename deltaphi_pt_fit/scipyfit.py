import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle

def func(x, a, b):
    return a * np.power(x, b)

def func2(x, a):
    return a/x

def func3(x, a, b):
    return a/(x+b)

def func4(x, a, b):
    return a/x + b

txt = np.loadtxt('points.txt')

names = [
    'dphiL123', 'dphiL124', 'dphiL134', 'dphiL234', 'dphiL12D1', 'dphiL13D1',
    'dphiL23D1', 'dphiL12D2', 'dphiL1D12', 'dphiL2D12', 'dphiL1D13', 'dphiL1D23',
    'dphiD123', 'dphiD124', 'dphiD134', 'dphiD234', 'dphiD235', 'dphiD245', 'dphiD345'
]

pp = PdfPages('/home/jongho/fitting.pdf')
for inum, iname in enumerate(names):
    #if not inum == 0: continue
    #print(iname)
    #print('x:', txt[4*inum])
    #print('x_err:', txt[4*inum+1])
    #print('y:', txt[4*inum+2])
    #print('y_err:', txt[4*inum+3])

    xdata = np.asarray(txt[4*inum])
    ydata = np.asarray(txt[4*inum+2])

    xdata10GeV = np.asarray(txt[4*inum][:19])
    ydata10GeV = np.asarray(txt[4*inum+2][:19])

    #popt, pocv = curve_fit(func, xdata, ydata, bounds=([0.0, -1.0001], [0.1, -0.99999]) )
    popt, pocv = curve_fit(func2, xdata, ydata)
    popt2, pocv2 = curve_fit(func2, xdata10GeV, ydata10GeV)
    #popt, pocv = curve_fit(func3, xdata, ydata)
    #popt, pocv = curve_fit(func4, xdata, ydata)
    print(iname, '{:.3f}'.format(popt[0]))

    ### Chi-square test in python
    #res = stats.chisquare(f_obs=ydata, f_exp=func2(xdata, *popt))
    res = stats.chisquare(f_exp=ydata, f_obs=func2(xdata, *popt))
    #res2 = stats.chisquare(f_obs=ydata10GeV, f_exp=func2(xdata10GeV, *popt2))
    res2 = stats.chisquare(f_exp=ydata10GeV, f_obs=func2(xdata10GeV, *popt2))
    #res = stats.chisquare(f_obs=ydata, f_exp=func3(xdata, *popt))
    #res = stats.chisquare(f_obs=ydata, f_exp=func4(xdata, *popt))
    print('chi-square value: {:.5f}, p-value {:.5f}'.format(res[0], res[1]))
    print('chi-square value: {:.5f}, p-value {:.5f}'.format(res2[0], res2[1]))
    print()

    #fig, ax = plt.subplots(1,1, figsize=(18,16))
    fig, ax = plt.subplots(figsize=(10, 8))
    #ax.locator_params(axis='x', tight=True, nbins=20)
    ax.grid(which='major')
    ax.grid(which='minor')
    #ax.xaxis.set_major_locator(plt.MaxNLocator(10))
    #ax.xaxis.set_minor_locator(plt.MaxNLocator(10))
    ax.errorbar(xdata, ydata, xerr=txt[4*inum+1], yerr=txt[4*inum+3], linestyle='None', marker='^', ms=2,  mfc='b', ecolor='k', elinewidth=0.3, capsize=2)
    #ax.plot(xdata, func(xdata, *popt), 'r-', label='fit: a=%2.3f, b=%2.2f' % tuple(popt))
    ax.plot(xdata, func2(xdata, *popt), 'r-', label='100 GeV fit: a=%.4f' % tuple(popt), color='r', linewidth=5, linestyle=':')
    ax.plot(xdata, func2(xdata, *popt), 'r-', label='10 GeV fit: a=%.4f' % tuple(popt2), color='g', linewidth=2)
    #ax.plot(xdata, func3(xdata, *popt), 'r-', label='fit: a=%.5f, b=%.5f' % tuple(popt))
    #ax.plot(xdata, func4(xdata, *popt), 'r-', label='fit: a=%.5f, b=%.5f' % tuple(popt))
    ax.set_xlabel('$\Delta\phi$', fontsize=15)
    ax.set_ylabel('Gen-level electron $p_T$', fontsize=15)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.0005, 0.2)
    #ax.set_xlim(None, 0.2)
    ax.set_ylim(None, 10)
    ax.legend(fontsize=15)

    ax.text(0.001, 2.5, r'$p_T = \frac{%.3f}{\Delta\phi}$' % tuple(popt), style='italic', fontsize=30, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})
    ax.text(0.001, 1.5, 'combination: %s' %(iname.split('dphi')[1]), style='italic', fontsize=15, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})
    ax.text(0.001, 0.75, r'$\chi^2/ndf  (100 GeV) : {:1.4f}$'.format(res[0]/(xdata.size-1)), style='italic', fontsize=12, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})
    ax.text(0.001, 0.5, r'$\chi^2/ndf  (10 GeV) : {:1.4f}$'.format(res2[0]/(xdata10GeV.size-1)), style='italic', fontsize=12, bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10})

    pp.savefig(fig)
    plt.close(fig)

pp.close()
