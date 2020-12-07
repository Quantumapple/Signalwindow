import ROOT as rt

rt.gStyle.SetOptStat(0)

c1 = rt.TCanvas('c1', '', 1000, 900)
c1.SetGrid()

sig = rt.TFile.Open('../option1/sig/sig_final.root')
bkg = rt.TFile.Open('../option1/bkg/bkg_final.root')

r1sig = sig.Get('h1')
r6sig = sig.Get('h6')
r1bkg = bkg.Get('h1')
r6bkg = bkg.Get('h6')

r1sig.SetDirectory(0)
r1bkg.SetDirectory(0)
r6sig.SetDirectory(0)
r6bkg.SetDirectory(0)

sig.Close()
bkg.Close()

r1sig.Scale(1/r1sig.Integral())
r1bkg.Scale(1/r1bkg.Integral())
r6sig.Scale(1/r6sig.Integral())
r6bkg.Scale(1/r6bkg.Integral())

g1 = rt.TGraph()
g2 = rt.TGraph()

for i in range(1, 52):
    g1.SetPoint(i-1, r1sig.Integral(1,i), 1. - r1bkg.Integral(1,i))
    g2.SetPoint(i-1, r6sig.Integral(1,i), 1. - r6bkg.Integral(1,i))

g1.SetMarkerStyle(22)
g2.SetMarkerStyle(22)
g1.SetMarkerSize(2.)
g2.SetMarkerSize(2.)
g1.SetMarkerColor(rt.kOrange)
g2.SetMarkerColor(rt.kOrange)

g1.SetLineStyle(9)
g2.SetLineStyle(9)
g1.SetLineWidth(3)
g2.SetLineWidth(3)
g1.SetLineColor(rt.kOrange)
g2.SetLineColor(rt.kOrange)

#### input ####
#finput = rt.TFile.Open('r1_roc.root')
#mva1 = finput.Get('r1_roc_BDTG')
#mva2 = finput.Get('r1_roc_MLPS')
#mva3 = finput.Get('r1_roc_DNN')
#mva4 = finput.Get('r1_roc_Ld')

finput = rt.TFile.Open('r6_roc.root')
mva1 = finput.Get('r6_roc_BDTG')
mva2 = finput.Get('r6_roc_MLPS')
mva3 = finput.Get('r6_roc_DNN')
mva4 = finput.Get('r6_roc_Ld')
###############

mva1.SetLineWidth(3)
mva2.SetLineWidth(3)
mva3.SetLineWidth(3)
mva4.SetLineWidth(3)

mva1.SetLineColor(rt.kBlack)
mva2.SetLineColor(rt.kRed)
mva3.SetLineColor(rt.kBlue)
mva4.SetLineColor(rt.kGreen+1)

mva1.SetLineStyle(1)
mva2.SetLineStyle(4)
mva3.SetLineStyle(9)
mva4.SetLineStyle(10)

mva1.SetTitle('')
mva1.SetMinimum(0.)
mva1.GetXaxis().SetTitle('Signal eff')
mva1.GetXaxis().SetRangeUser(0.90, 1.0)
mva1.GetYaxis().SetTitle('Backgr rejection (1-eff)')

mva1.Draw('ac')
mva2.Draw('csame')
mva3.Draw('csame')
mva4.Draw('csame')
g2.Draw('csame')

# legend
x0L = 0.15
y0H = 0.899
dxL = 0.557-x0L
dyH = 0.22
y0H = 1 - y0H + dyH + 0.07

legend = rt.TLegend( x0L, y0H-dyH, x0L+dxL, y0H )
legend.SetHeader( "MVA Method:" )
legend.SetMargin( 0.4 )
legend.AddEntry(g1, "Cut based", "l")
legend.AddEntry(mva1, "BDT", "l")
legend.AddEntry(mva2, "MLP", "l")
legend.AddEntry(mva3, "Linear Discriminator", "l")
legend.AddEntry(mva4, "DNN", "l")
legend.Draw("same")

c1.Update()
c1.Print('roc_r1.png')
