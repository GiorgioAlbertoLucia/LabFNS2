import pandas as pd
import numpy as np

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure, kFullCircle, kFullSquare, gROOT
from DfUtils import GetGraphErrorsFromCSV
from StyleFormatter import SetObjectStyle, SetGlobalStyle

if __name__=='__main__':
    
    gROOT.SetBatch()
    SetGlobalStyle(padleftmargin=0.12, padbottommargin=0.12, padrightmargin=0.05, padtopmargin=0.1, titleoffsety=1.2, titleoffsetx=0.9, titleoffset= 0.7, opttitle=1)

    infileLGAD = 'Probe-station/data/input/C_vs_f_LGAD.csv'
    infilePin = 'Probe-station/data/input/C_vs_f_pin.csv'
    #infileLGAD = 'Probe-station/data/input/C_vs_f_LGAD_ex.csv'
    #infilePin = 'Probe-station/data/input/C_vs_f_LGAD_ex.csv'
    outfilename = 'Probe-station/data/output/C_vs_f.root'
    
    gLGAD = GetGraphErrorsFromCSV(infileLGAD)
    gLGAD.SetName("gLGAD")
    gLGAD.SetTitle("LGAD; Frequency [kHz]; Capacitance (pF)")
    gLGAD.SetMarkerStyle(kFullCircle)
    gLGAD.SetMarkerSize(1)
    gLGAD.SetMarkerColor(kRed+1)
    gLGAD.SetLineColor(kRed+1)
    
    gPin = GetGraphErrorsFromCSV(infilePin)
    gPin.SetName("gPin")
    gPin.SetTitle("Pin; Frequency [kHz]; Capacitance (pF)")
    gPin.SetMarkerStyle(kFullSquare)
    gPin.SetMarkerSize(1)
    gPin.SetMarkerColor(kAzure + 3)
    gPin.SetLineColor(kAzure + 3)
    
    canvas = TCanvas("canvas","canvas",1000,1000)
    hFrame = canvas.cd().DrawFrame(0,1,50,100,"C-f curves of FBK-UFSD2 LGAD and PiN pads; Frequency [kHz]; Capacitance (pF)")
    canvas.SetLogx()
    gLGAD.Draw("pa")
    gPin.Draw("pa")

    legend = TLegend(0.7, 0.7, 0.86, 0.8)
    legend.AddEntry(gLGAD,'LGAD','p')
    legend.AddEntry(gPin,'PiN','p')
    legend.Draw("same")

    canvas.Modified()
    canvas.Update()
    
    canvas.SaveAs('Probe-station/data/output/C_vs_f.pdf')

    outfile = TFile(outfilename, 'recreate')
    gLGAD.Write()
    gPin.Write()
    canvas.Write()
    outfile.Close()
