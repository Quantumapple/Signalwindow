import ROOT as rt
import numpy as np

#rt.gStyle.SetOptStat(1110)
#rt.gStyle.SetOptFit(1111)

txt = np.loadtxt('points.txt')

names = [
    'dphiL123', 'dphiL124', 'dphiL134', 'dphiL234', 'dphiL12D1', 'dphiL13D1',
    'dphiL23D1', 'dphiL12D2', 'dphiL1D12', 'dphiL2D12', 'dphiL1D13', 'dphiL1D23',
    'dphiD123', 'dphiD124', 'dphiD134', 'dphiD234', 'dphiD235', 'dphiD245', 'dphiD345'
]

c1 = rt.TCanvas('c1', '', 1200, 1000)
c1.SetLeftMargin(0.14)
c1.SetBottomMargin(0.13)
c1.SetLogx()
c1.SetLogy()
c1.SetGrid()

for inum, iname in enumerate(names):
    #if not inum == 0: continue
    #print(iname)
    #print('x:', txt[4*inum])
    #print('x_err:', txt[4*inum+1])
    #print('y:', txt[4*inum+2])
    #print('y_err:', txt[4*inum+3])
    #print('')

    x = txt[4*inum]
    xerr = txt[4*inum+1]
    y = txt[4*inum+2]
    yerr = txt[4*inum+3]

    #print(y[:19])

    g1 = rt.TGraphErrors()
    for i in range(20):
        g1.SetPoint(i, x[i], y[i])
        g1.SetPointError(i, xerr[i], yerr[i])

    g1.SetMarkerStyle(22)
    g1.SetMarkerSize(2.4)
    g1.SetMarkerColor(rt.kAzure)
    g1.GetXaxis().SetTitle('#Delta#phi')
    g1.GetXaxis().CenterTitle()
    g1.GetXaxis().SetMoreLogLabels()
    g1.GetXaxis().SetTitleOffset(1.3)
    g1.GetYaxis().SetTitle('Gen-level p_{T} (GeV)')
    g1.GetYaxis().SetMoreLogLabels()
    g1.SetMinimum(0.05)
    g1.Draw("APE")
    c1.Update()

    f1 = rt.TF1('f1', '[0]/x', 0.001, 0.4)
    g1.Fit(f1)
    print('1st parameter {:3.4f}'.format(f1.GetParameter(0)))

    r_ = c1.GetLeftMargin()
    t_ = c1.GetBottomMargin()
    combi = rt.TLatex(r_+0.3, t_+0.45*t_, iname.split('dphi')[1])
    combi.SetNDC(True)
    combi.SetTextFont(42)
    combi.SetTextAlign(31)
    combi.SetTextSize(0.5*t_)
    combi.Draw()

    c1.Update()

    path = '/home/jongho/temp/roofit/'
    filename = iname.split('dphi')[1]+'_roofit'

    c1.Print(path+filename+'.png')

c1.Close()
