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
    hTimeFit=TH1D("hTimeFit","hTimeFit",46,-350,550.)
    theshold2=80
    theshold3=80
    for a2,a3,t2,t3,t2fit,t3fit in zip(df["Amp2"],df["Amp3"],df["ToA2"],df["ToA3"], df["ToA2f"], df["ToA3f"]): 
        if (a2 >= theshold2 and a3>= theshold3):
            hTime.Fill(1000*(t2-t3))
            hTimeFit.Fill(1000*(t2fit-t3fit))
    gaus=TF1("gaus","[0]*exp(-((x-[1])^2)/(2*[2]^2))",-100,300)
    gaus.SetLineColor(kRed)
    gaus.SetParameter(2,50)
    gaus.SetParameter(1,80)
    gaus.SetParameter(0,3500)
    gStyle.SetOptFit(11111111)
    
    hTime.Fit(gaus,"rml+")
    timeres=gaus.GetParameter(2)
    print("par2=",gaus.GetParameter(2))
    gaus2=TF1("gaus2","[0]*exp(-((x-[1])^2)/(2*[2]^2))",0,240)
    gaus2.SetLineColor(kAzure)
    gaus2.SetParameter(2,50)
    gaus2.SetParameter(1,80)
    gaus2.SetParameter(0,1800)
    gStyle.SetOptFit(11111111)
    hTimeFit.Fit(gaus2,"rml+")
    timeresfit=gaus2.GetParameter(2)
   
    print("risolution UFSD: ",timeres/np.sqrt(2))
    print("risolution UFSD with fit: ",timeresfit/np.sqrt(2))
    canvas=TCanvas("canvas","canvas",1900,1500)
    canvas.Divide(2,1)
    canvas.cd(1)
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
    print("chi/dof : ",gaus.GetChisquare(),"/",gaus.GetNDF())
    canvas.cd(2)
    hTimeFit.SetTitle("#Delta_{t} distribution, theshold:80 mV ")
    hTimeFit.GetXaxis().SetTitle("#Delta_{t} (ps)")
    hTimeFit.GetYaxis().SetTitle("Counts")
    hTimeFit.Draw("hist")
    gaus2.Draw("same")
    legend2 = TLegend(0.6, 0.3, 0.85, 0.5)
    legend2.AddEntry(hTimeFit,'Data','l')
    legend2.AddEntry(gaus2,'p0*exp(#frac{-(x-p1)^{2}}{2*p2^{2}})','l')
    legend2.SetTextSize(0.045)
    legend2.Draw("same")
    print("chi/dof fit: ",gaus2.GetChisquare(),"/",gaus2.GetNDF())
    print("input")
    input()
    outfile = TFile(outfilename, 'recreate')
    hTimeFit.Write()
    hTime.Write()
    gaus.Write()
    gaus2.Write()
    canvas.Write()
    canvas.SaveAs('Beta/data/output/Time_distribution.pdf')
    outfile.Close()