import pandas as pd
import numpy as np
import uproot

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure,  kOrange, kFullCircle, kFullSquare, gROOT, TF1, TGaxis, TTree, TH1D, TH2D, gStyle
from DfUtils import GetGraphErrorsFromCSV
from StyleFormatter import SetObjectStyle, SetGlobalStyle
import sys

if __name__=='__main__':
    tree=uproot.open("Beta/data/output/BetaOutput.root")["BetaTree"]
    outfilename = 'Beta/data/output/Time_resolution.root'
    df=tree.arrays(library='pd')
    hTime=TH1D("hTime","hTime",19,-350,550.)
    theshold2=80
    theshold3=80
    for a2,a3,t2,t3 in zip(df["Amp2"],df["Amp3"],df["ToA2"],df["ToA3"]): 
        if (a2 >= theshold2 and a3>= theshold3):
            hTime.Fill(1000*(t2-t3))
            print(t2-t3)
    gaus=TF1("gaus","[0]*exp(-((x-[1])^2)/(2*[2]^2))",-100,300)
    gaus.SetLineColor(kRed)
    #gaus.SetParLimits(0,3100,4000)
    gaus.SetParameter(2,50)
    gaus.SetParameter(1,80)
    gaus.SetParameter(0,3500)
    hTime.Fit(gaus,"rml+")
    gStyle.SetOptFit(11111111)
    timeres=gaus.GetParameter(2)
    print("risolution UFSD: ",timeres/np.sqrt(2))
    canvas=TCanvas("canvas","canvas",1900,1500)
    canvas.cd()
    hTime.SetTitle("#Delta_{t} distribution, theshold:80 mV ")
    hTime.GetXaxis().SetTitle("#Delta_{t} (ps)")
    hTime.GetYaxis().SetTitle("Counts")
    hTime.Draw("hist")
    gaus.Draw("same")
    legend1 = TLegend(0.6, 0.3, 0.85, 0.5)
    legend1.AddEntry(hTime,'Data','l')
    legend1.AddEntry(gaus,'p0*exp(#frac{-(x-p1)^{2}}{2*p2^{2}})','l')
    legend1.SetTextSize(0.045)
    legend1.Draw("same")
    print("chi/dof: ",gaus.GetChisquare(),"/",gaus.GetNDF())
    print("input")
    input()
    outfile = TFile(outfilename, 'recreate')
    hTime.Write()
    gaus.Write()
    canvas.Write()
    canvas.SaveAs('Beta/data/output/Time_distribution.pdf')
    outfile.Close()