import pandas as pd
import numpy as np

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure, kOrange, kFullCircle,kOpenDiamond, kFullSquare, gROOT
from DfUtils import GetGraphErrorsFromCSV
from StyleFormatter import SetObjectStyle, SetGlobalStyle

def ProduceGraphJittervsGain(infile, DetectorName, GraphColor, GraphMarker):
    df=pd.read_csv(infile,comment='#',sep='\t')
    df['Jitter']=df['N']/df['dV/dt']
    df['err_Jitter']=df['Jitter']*np.sqrt((df['err_N']/df['N'])**2+(df['err_dV/dt']/df['dV/dt'])**2)
    df['Gain']=df['AtotLGAD']/df['AtotPin']
    df['err_Gain']=df['Gain']*np.sqrt((df['err_AtotLGAD']/df['AtotLGAD'])**2+(df['err_AtotPin']/df['AtotPin'])**2)
    gJittervsGain=TGraphErrors(len(df['Gain']), np.asarray(df['Gain'], dtype=float), np.asarray(df['Jitter'], dtype=float), np.asarray(df['err_Gain'], dtype=float), np.asarray(df['err_Jitter'], dtype=float))
    gJittervsGain.SetName("gJittervsGain"+DetectorName)
    gJittervsGain.SetTitle(";Gain; Jitter (ns)")
    gJittervsGain.SetMarkerStyle(GraphMarker)
    gJittervsGain.SetMarkerSize(1)
    gJittervsGain.SetMarkerColor(GraphColor)
    gJittervsGain.SetLineColor(GraphColor)
    return gJittervsGain


if __name__=='__main__':
    
    gROOT.SetBatch()
    SetGlobalStyle(padleftmargin=0.12, padbottommargin=0.12, padrightmargin=0.05, padtopmargin=0.1, titleoffsety=1.2, titleoffsetx=0.9, titleoffset= 0.7, opttitle=1)

    infiles = ['TCT/data/input/Jitter_vs_gain.csv']*2
    DetectorNames = ['LGAD']*2
    GraphColors = [kRed+1,kAzure+3]
    GraphMarkers = [kFullCircle, kOpenDiamond]
    outfilename = 'TCT/data/output/Jitter_vs_gain.root'
    legendmax = 0.85
    legendstep = 0.05
    
    canvas = TCanvas("canvas","canvas",1000,1000)
    hFrame = canvas.cd().DrawFrame(0,0,20,1,";Gain; Jitter (ns)")
    canvas.SetLogx()

    Graphs = []
    legend = TLegend(0.6, legendmax-len(infiles)*legendstep, 0.8, legendmax)
    legend.SetTextSize(0.03)
    for file, name, color, marker in zip(infiles, DetectorNames, GraphColors, GraphMarkers):
        Graphs.append(ProduceGraphJittervsGain(file, name, color, marker))
        Graphs[-1].Draw("p,same")
        legend.AddEntry(Graphs[-1],name,'p')

    legend.Draw("same")

    canvas.Modified()
    canvas.Update()
    
    canvas.SaveAs('TCT/data/output/Jitter_vs_gain.pdf')

    outfile = TFile(outfilename, 'recreate')
    for graph in Graphs:
        graph.Write()
    canvas.Write()
    outfile.Close()
