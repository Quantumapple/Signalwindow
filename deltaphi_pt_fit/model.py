import ROOT as rt
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm
import pickle

infile = rt.TFile.Open('results.root')
tree = infile.Get('l1PiXTRKTree/L1PiXTRKTree')

def deltaPhi(list1, list2, list3):
    v1 = rt.TVector3(list2[0][0] - list1[0][0], list2[0][1] - list1[0][1], list2[0][2] - list1[0][2])
    v2 = rt.TVector3(list3[0][0] - list2[0][0], list3[0][1] - list2[0][1], list3[0][2] - list2[0][2])

    dphi = v2.Phi() - v1.Phi()
    if dphi >= math.pi:
        dphi -= 2.*math.pi
    elif dphi < -math.pi:
        dphi += 2.*math.pi

    #return abs(dphi)
    return dphi

def errorPoints(pt, combi):
    ### Define outputs ###
    x = []
    x_err = []
    y = []
    y_err = []

    for k in range(len(pt)-1):
        #if not k == 0: continue
        dphi = []
        #print(pt[k], pt[k+1])

        for j in range(len(combi)):
            if combi[j][0] > pt[k] and combi[j][0] < pt[k+1]:
                dphi.append(combi[j][1])

        mu, std = norm.fit(np.asarray(dphi))
        x.append(abs(mu))
        x_err.append(std)
        y.append((pt[k]+pt[k+1])/2.)
        y_err.append(pt[k+1] - ((pt[k]+pt[k+1])/2.))
        #if pt[k] <= 1.:
        #    #### Median of array and its standard error ####
        #    #print('value:', np.median(dphi), 'error:', 1.253*np.std(dphi)/np.sqrt(len(dphi)))
        #    x.append(np.median(dphi))
        #    x_err.append(1.253*np.std(dphi)/np.sqrt(len(dphi)))
        #    y.append((pt[k]+pt[k+1])/2.)

        #if pt[k] > 1.:
        #    mu, std = norm.fit(np.asarray(dphi))
        #    x.append(mu)
        #    x_err.append(std)
        #    y.append((pt[k]+pt[k+1])/2.)
        #    #print('value from fit', mu, 'error (std)', std)
        #    #print()

    return x, x_err, y, y_err

def drawHist(ax, ptmin, ptmax, combi, dofit=True):

    arr = []
    for j in range(len(combi)):
        if combi[j][0] > ptmin and combi[j][0] < ptmax:
            arr.append(combi[j][1])

    ax.hist(arr, bins=50, density=True, alpha=0.6)
    ax.set(xlabel='$\Delta\phi$ between gen-level $p_T$ {:1.1f} (GeV) and {:1.1f} (GeV)'.format(ptmin, ptmax))
    xmin, xmax = ax.get_xlim()
    x = np.linspace(xmin, xmax, 100)

    if dofit:
        try:
            mu, std = norm.fit(np.asarray(arr))
            pdf = norm.pdf(x, mu, std)
            ax.plot(x, pdf, 'r', linewidth=2, label='Fit results mu = %.5f, std = %.5f' % (mu, std))
            ax.legend()
        except:
            pass

    return ax

#### delta-phi lists ####
dphiL123 = []
dphiL124 = []
dphiL134 = []
dphiL234 = []

dphiL12D1 = []
dphiL13D1 = []
dphiL23D1 = []
dphiL12D2 = []
dphiL1D12 = []
dphiL2D12 = []
dphiL1D13 = []
dphiL1D23 = []

dphiD123 = []
dphiD124 = []
dphiD134 = []
dphiD234 = []
dphiD235 = []
dphiD245 = []
dphiD345 = []

event = 0
##### Divide pixels into each lists #####
for entry in tree:

    #print(str(event)+'th event')
    if event % 100000 == 0: print(str(event)+'th event')
    #if event > 10000: break

    bp1st = []
    bp2nd = []
    bp3rd = []
    bp4th = []

    fp1st = []
    fp2nd = []
    fp3rd = []
    fp4th = []
    fp5th = []

    #### Declare branches ####
    bx = entry.bRecHitGx
    by = entry.bRecHitGy
    bz = entry.bRecHitGz
    fx = entry.fRecHitGx
    fy = entry.fRecHitGy
    fz = entry.fRecHitGz

    layer = entry.bRecHitLayer
    disk = entry.fRecHitDisk

    genpt = entry.propgenElPartPt

    for i in range(bx.size()):
        if layer.at(i) == 1:
            bp1st.append([ bx.at(i), by.at(i), bz.at(i) ])
        if layer.at(i) == 2:
            bp2nd.append([ bx.at(i), by.at(i), bz.at(i) ])
        if layer.at(i) == 3:
            bp3rd.append([ bx.at(i), by.at(i), bz.at(i) ])
        if layer.at(i) == 4:
            bp4th.append([ bx.at(i), by.at(i), bz.at(i) ])

    for j in range(fx.size()):
        if disk.at(j) == 1:
            fp1st.append([ fx.at(j), fy.at(j), fz.at(j) ])
        if disk.at(j) == 2:
            fp2nd.append([ fx.at(j), fy.at(j), fz.at(j) ])
        if disk.at(j) == 3:
            fp3rd.append([ fx.at(j), fy.at(j), fz.at(j) ])
        if disk.at(j) == 4:
            fp4th.append([ fx.at(j), fy.at(j), fz.at(j) ])
        if disk.at(j) == 5:
            fp5th.append([ fx.at(j), fy.at(j), fz.at(j) ])

    #### Playing with pixel lists ####

    #print(len(bp1st), len(bp2nd), len(bp3rd), len(bp4th))
    #print(len(fp1st), len(fp2nd), len(fp3rd), len(fp4th), len(fp5th))
    #print()

    ### region1 
    if len(bp1st)+len(bp2nd)+len(bp3rd) >= 3:
        if genpt.at(0) < 100.:
            dphiL123.append([ genpt.at(0), deltaPhi(bp1st, bp2nd, bp3rd) ])
    if len(bp1st)+len(bp2nd)+len(bp4th) >= 3:
        if genpt.at(0) < 100.:
            dphiL124.append([ genpt.at(0), deltaPhi(bp1st, bp2nd, bp4th) ])
    if len(bp1st)+len(bp3rd)+len(bp4th) >= 3:
        if genpt.at(0) < 100.:
            dphiL134.append([ genpt.at(0), deltaPhi(bp1st, bp3rd, bp4th) ])
    if len(bp2nd)+len(bp3rd)+len(bp4th) >= 3:
        if genpt.at(0) < 100.:
            dphiL234.append([ genpt.at(0), deltaPhi(bp2nd, bp3rd, bp4th) ])


    ### region2, region3, region4
    if len(bp1st)+len(bp2nd)+len(fp1st) >= 3:
        if genpt.at(0) < 100.:
            dphiL12D1.append([ genpt.at(0), deltaPhi(bp1st, bp2nd, fp1st) ])
    if len(bp1st)+len(bp3rd)+len(fp1st) >= 3:
        if genpt.at(0) < 100.:
            dphiL13D1.append([ genpt.at(0), deltaPhi(bp1st, bp3rd, fp1st) ])
    if len(bp2nd)+len(bp3rd)+len(fp1st) >= 3:
        if genpt.at(0) < 100.:
            dphiL23D1.append([ genpt.at(0), deltaPhi(bp2nd, bp3rd, fp1st) ])
    if len(bp1st)+len(bp2nd)+len(fp2nd) >= 3:
        if genpt.at(0) < 100.:
            dphiL12D2.append([ genpt.at(0), deltaPhi(bp1st, bp2nd, fp2nd) ])
    if len(bp1st)+len(fp1st)+len(fp2nd) >= 3:
        if genpt.at(0) < 100.:
            dphiL1D12.append([ genpt.at(0), deltaPhi(bp1st, fp1st, fp2nd) ])
    if len(bp2nd)+len(fp1st)+len(fp2nd) >= 3:
        if genpt.at(0) < 100.:
            dphiL2D12.append([ genpt.at(0), deltaPhi(bp2nd, fp1st, fp2nd) ])
    if len(bp1st)+len(fp1st)+len(fp3rd) >= 3:
        if genpt.at(0) < 100.:
            dphiL1D13.append([ genpt.at(0), deltaPhi(bp1st, fp1st, fp3rd) ])
    if len(bp1st)+len(fp2nd)+len(fp3rd) >= 3:
        if genpt.at(0) < 100.:
            dphiL1D23.append([ genpt.at(0), deltaPhi(bp1st, fp2nd, fp3rd) ])

    ### region5, region6
    if len(fp1st)+len(fp2nd)+len(fp3rd) >= 3:
        if genpt.at(0) < 100.:
            dphiD123.append([ genpt.at(0), deltaPhi(fp1st, fp2nd, fp3rd) ])
    if len(fp1st)+len(fp2nd)+len(fp4th) >= 3:
        if genpt.at(0) < 100.:
            dphiD124.append([ genpt.at(0), deltaPhi(fp1st, fp2nd, fp4th) ])
    if len(fp1st)+len(fp3rd)+len(fp4th) >= 3:
        if genpt.at(0) < 100.:
            dphiD134.append([ genpt.at(0), deltaPhi(fp1st, fp3rd, fp4th) ])
    if len(fp2nd)+len(fp3rd)+len(fp4th) >= 3:
        if genpt.at(0) < 100.:
            dphiD234.append([ genpt.at(0), deltaPhi(fp2nd, fp3rd, fp4th) ])
    if len(fp2nd)+len(fp3rd)+len(fp5th) >= 3:
        if genpt.at(0) < 100.:
            dphiD235.append([ genpt.at(0), deltaPhi(fp2nd, fp3rd, fp5th) ])
    if len(fp2nd)+len(fp4th)+len(fp5th) >= 3:
        if genpt.at(0) < 100.:
            dphiD245.append([ genpt.at(0), deltaPhi(fp2nd, fp4th, fp5th) ])
    if len(fp3rd)+len(fp4th)+len(fp5th) >= 3:
        if genpt.at(0) < 100.:
            dphiD345.append([ genpt.at(0), deltaPhi(fp3rd, fp4th, fp5th) ])

    event += 1

#### Fitting ####

ptrange = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
for pt in range(1, 101):
    if pt <= 10:
        ptrange.append(pt)
    if pt > 10 and pt <= 30 and pt % 2 == 0:
        ptrange.append(pt)
    if pt > 30 and pt <= 50 and pt % 5 == 0:
        ptrange.append(pt)
    if pt > 50 and pt <= 100 and pt % 10 == 0:
        ptrange.append(pt)

print(ptrange)

dphilists = [
    dphiL123, dphiL124, dphiL134, dphiL234, dphiL12D1, dphiL13D1, dphiL23D1,
    dphiL12D2, dphiL1D12, dphiL2D12, dphiL1D13, dphiL1D23, dphiD123, dphiD124,
    dphiD134, dphiD234, dphiD235, dphiD245, dphiD345
]

names = [
    'dphiL123', 'dphiL124', 'dphiL134', 'dphiL234', 'dphiL12D1', 'dphiL13D1',
    'dphiL23D1', 'dphiL12D2', 'dphiL1D12', 'dphiL2D12', 'dphiL1D13', 'dphiL1D23',
    'dphiD123', 'dphiD124', 'dphiD134', 'dphiD234', 'dphiD235', 'dphiD245', 'dphiD345'
]

with open('points.txt', 'wb') as filehandle:

    for inum, idphi in enumerate(dphilists):
        #if not inum < 2: continue
        x = []
        x_err = []
        y = []
        y_err = []
        if len(idphi) > 0:
            x, x_err, y, y_err = errorPoints(ptrange, idphi)
        else:
            continue

        np.savetxt(filehandle, [np.asarray(x)])
        np.savetxt(filehandle, [np.asarray(x_err)])
        np.savetxt(filehandle, [np.asarray(y)], fmt='%3.3f')
        np.savetxt(filehandle, [np.asarray(y_err)], fmt='%3.3f')

pp = PdfPages('/home/jongho/example_l12d1.pdf')
for ival in range(len(ptrange)-1):
    fig, ax = plt.subplots()
    plt.figure(figsize=(12,8))
    minval = ptrange[ival]
    maxval = ptrange[ival+1]
    drawHist(ax, minval, maxval, dphiL12D1)
    pp.savefig(fig)
    plt.close('all')
pp.close()

pp = PdfPages('/home/jongho/example_l13d1.pdf')
for ival in range(len(ptrange)-1):
    fig, ax = plt.subplots()
    plt.figure(figsize=(12,8))
    minval = ptrange[ival]
    maxval = ptrange[ival+1]
    drawHist(ax, minval, maxval, dphiL13D1)
    pp.savefig(fig)
    plt.close('all')
pp.close()

pp = PdfPages('/home/jongho/example_l1d12.pdf')
for ival in range(len(ptrange)-1):
    fig, ax = plt.subplots()
    plt.figure(figsize=(12,8))
    minval = ptrange[ival]
    maxval = ptrange[ival+1]
    drawHist(ax, minval, maxval, dphiL1D12)
    pp.savefig(fig)
    plt.close('all')
pp.close()

pp = PdfPages('/home/jongho/example_l2d12.pdf')
for ival in range(len(ptrange)-1):
    fig, ax = plt.subplots()
    plt.figure(figsize=(12,8))
    minval = ptrange[ival]
    maxval = ptrange[ival+1]
    drawHist(ax, minval, maxval, dphiL2D12)
    pp.savefig(fig)
    plt.close('all')
pp.close()

#plt.figure(figsize=(12,8))
#plt.grid()
#plt.errorbar(x, y, yerr=None, xerr=x_err, linestyle='None', marker='^', ms=5, mec='r', mfc='r', ecolor='k', elinewidth=2)
#plt.ylim(0.05, 100)
#plt.ylim(0.05, 10)
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
