import ROOT as rt
import numpy as np

infile = rt.TFile.Open('../deltaPhi_pt_fit/results.root')
tree = infile.Get('l1PiXTRKTree/L1PiXTRKTree')

#### parameters ####
B = 3.8 ### (T)
speedoflight = 299792458 ### (m/s)
#e_charge = 1.6e-19
e_charge = -1
e_mass = 0.511e-3
####################

def momentum(genphi, list1, list2, bfield, sol):

    if genphi < 0.:
        genphi += np.pi

    rot_x1 = list1[0][0]*np.cos(-genphi) - list1[0][1]*np.sin(-genphi)
    rot_x2 = list2[0][0]*np.cos(-genphi) - list2[0][1]*np.sin(-genphi)
    rot_y1 = list1[0][0]*np.sin(-genphi) + list1[0][1]*np.cos(-genphi)
    rot_y2 = list2[0][0]*np.sin(-genphi) + list2[0][1]*np.cos(-genphi)

    r1 = np.sqrt(rot_x1**2+rot_y1**2)
    r2 = np.sqrt(rot_x2**2+rot_y2**2)

    radius1 = (r1**2 - r2**2)/(2.*(rot_y1 - rot_y2))
    radius2 = (r1**2 + r2**2)/(2.*(rot_y1 + rot_y2))

    pt1 = bfield*abs(radius1)*sol*1e-11
    pt2 = bfield*abs(radius2)*sol*1e-11

    return pt1, pt2

event = 0
for entry in tree:

    if event > 30: break

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

    genpt = entry.propgenElPartPt.at(0)
    genphi = entry.propgenElPartPhi.at(0)
    geneta = entry.propgenElPartEta.at(0)

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

    if genpt > 100 or genpt < 0.5:
        continue

    if abs(geneta) > 3.0:
        continue

    print('Gen-level eta:', geneta)
    ### region1 
    if len(bp1st)+len(bp2nd)+len(bp3rd) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, bp3rd, B, speedoflight)
        print('L123', pt1, pt2, 'genpt:', genpt)
    if len(bp1st)+len(bp2nd)+len(bp4th) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, bp4th, B, speedoflight)
        print('L124', pt1, pt2, 'genpt:', genpt)
    if len(bp1st)+len(bp3rd)+len(bp4th) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, bp4th, B, speedoflight)
        print('L134', pt1, pt2, 'genpt:', genpt)
    if len(bp2nd)+len(bp3rd)+len(bp4th) >= 3:
        pt1, pt2 = momentum(genphi, bp2nd, bp4th, B, speedoflight)
        print('L234', pt1, pt2, 'genpt:', genpt)

    ### region2, region3, region4
    if len(bp1st)+len(bp2nd)+len(fp1st) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, fp1st, B, speedoflight)
        print('L12D1', pt1, pt2, 'genpt:', genpt)
    if len(bp1st)+len(bp3rd)+len(fp1st) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, fp1st, B, speedoflight)
        print('L13D1', pt1, pt2, 'genpt:', genpt)
    if len(bp2nd)+len(bp3rd)+len(fp1st) >= 3:
        pt1, pt2 = momentum(genphi, bp2nd, fp1st, B, speedoflight)
        print('L23D1', pt1, pt2, 'genpt:', genpt)
    if len(bp1st)+len(bp2nd)+len(fp2nd) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, fp2nd, B, speedoflight)
        print('L12D2', pt1, pt2, 'genpt:', genpt)
    if len(bp1st)+len(fp1st)+len(fp2nd) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, fp2nd, B, speedoflight)
        print('L1D12', pt1, pt2, 'genpt:', genpt)
    if len(bp2nd)+len(fp1st)+len(fp2nd) >= 3:
        pt1, pt2 = momentum(genphi, bp2nd, fp2nd, B, speedoflight)
        print('L2D12', pt1, pt2, 'genpt:', genpt)
    if len(bp1st)+len(fp1st)+len(fp3rd) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, fp3rd, B, speedoflight)
        print('L1D13', pt1, pt2, 'genpt:', genpt)
    if len(bp1st)+len(fp2nd)+len(fp3rd) >= 3:
        pt1, pt2 = momentum(genphi, bp1st, fp3rd, B, speedoflight)
        print('L1D23', pt1, pt2, 'genpt:', genpt)

    ### region5, region6
    if len(fp1st)+len(fp2nd)+len(fp3rd) >= 3:
        pt1, pt2 = momentum(genphi, fp1st, fp3rd, B, speedoflight)
        print('D123', pt1, pt2, 'genpt:', genpt)
    if len(fp1st)+len(fp2nd)+len(fp4th) >= 3:
        pt1, pt2 = momentum(genphi, fp1st, fp4th, B, speedoflight)
        print('D124', pt1, pt2, 'genpt:', genpt)
    if len(fp1st)+len(fp3rd)+len(fp4th) >= 3:
        pt1, pt2 = momentum(genphi, fp1st, fp4th, B, speedoflight)
        print('D134', pt1, pt2, 'genpt:', genpt)
    if len(fp2nd)+len(fp3rd)+len(fp4th) >= 3:
        pt1, pt2 = momentum(genphi, fp2nd, fp4th, B, speedoflight)
        print('D234', pt1, pt2, 'genpt:', genpt)
    if len(fp2nd)+len(fp3rd)+len(fp5th) >= 3:
        pt1, pt2 = momentum(genphi, fp2nd, fp5th, B, speedoflight)
        print('D235', pt1, pt2, 'genpt:', genpt)
    if len(fp2nd)+len(fp4th)+len(fp5th) >= 3:
        pt1, pt2 = momentum(genphi, fp2nd, fp5th, B, speedoflight)
        print('D245', pt1, pt2, 'genpt:', genpt)
    if len(fp3rd)+len(fp4th)+len(fp5th) >= 3:
        pt1, pt2 = momentum(genphi, fp3rd, fp5th, B, speedoflight)
        print('D345', pt1, pt2, 'genpt:', genpt)

    print()
    event += 1
