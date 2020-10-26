import os
import ROOT as rt
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
import scipy.stats as stats

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

    return abs(dphi)
    #return dphi

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
    if event > 1000000: break

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

dphilists = [
    dphiL123, dphiL124, dphiL134, dphiL234, dphiL12D1, dphiL13D1, dphiL23D1,
    dphiL12D2, dphiL1D12, dphiL2D12, dphiL1D13, dphiL1D23, dphiD123, dphiD124,
    dphiD134, dphiD234, dphiD235, dphiD245, dphiD345
]

path = '/home/jongho/Physics/Delphes/pixel_delphes/deltaPhi_pt_fit/points'

names = [
    'dphiL123', 'dphiL124', 'dphiL134', 'dphiL234', 'dphiL12D1', 'dphiL13D1',
    'dphiL23D1', 'dphiL12D2', 'dphiL1D12', 'dphiL2D12', 'dphiL1D13', 'dphiL1D23',
    'dphiD123', 'dphiD124', 'dphiD134', 'dphiD234', 'dphiD235', 'dphiD245', 'dphiD345'
]

try:
    os.mkdir(path)
except:
    pass

for inum, idphi in enumerate(dphilists):
    filename = names[inum]+'.txt'
    with open(path+'/'+filename, 'wb') as filehandle:
        dphi = []
        pt = []

        for ival in range(len(idphi)):
            if ival > 3000: break
            if idphi[ival][0] > 10. or idphi[ival][0] < 0.5: continue
            dphi.append(idphi[ival][1])
            pt.append(idphi[ival][0])

        np.savetxt(filehandle, [np.asarray(dphi)])
        np.savetxt(filehandle, [np.asarray(pt)], fmt='%3.2f')
