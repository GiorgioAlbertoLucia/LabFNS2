import pandas as pd
import numpy as np
import uproot

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure,  kOrange, kFullCircle, kFullSquare, gROOT, TF1, TGaxis, TTree, TH1D, TH2D
from DfUtils import GetGraphErrorsFromCSV
from StyleFormatter import SetObjectStyle, SetGlobalStyle
import sys

if __name__=='__main__':
    tree=uproot.open("Beta/data/input/BetaOutput.root")["BetaTree"]
    outfilename = 'Beta/data/output/Time_resolution.root'
    df=tree.arrays(library='pd')
    hTime=TH1D("hTime","hTime",150,25.6,25.9)
    theshold=50
    for a2,a3,t2,t3 in zip(df["Amp2"],df["Amp3"],df["ToA2"],df["ToA3"]): 
        #if (a2 >= theshold and a3>= theshold):
        hTime.Fill(abs(t2-t3))
    timeres=hTime.GetRMS()
    print("risolution UFSD: ",timeres/np.sqrt(2))
    canvas=TCanvas("canvas","canvas",1500,1500)
    canvas.cd()
    hTime.GetXaxis().SetTitle("Time (ns)")
    hTime.GetYaxis().SetTitle("Counts")
    hTime.Draw()

    outfile = TFile(outfilename, 'recreate')
    hTime.Write()
    canvas.Write()
    canvas.SaveAs('Beta/data/output/Time_distribution.pdf')
    outfile.Close()