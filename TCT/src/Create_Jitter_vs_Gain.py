import pandas as pd
import numpy as np

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure, kOrange, kFullCircle,kOpenDiamond, kFullSquare, gROOT
from DfUtils import GetGraphErrorsFromCSV, GetGraphErrorsFromPandas
from StyleFormatter import SetObjectStyle, SetGlobalStyle

def ProduceGraphJittervsGain(infile, DetectorName, dfGain, GraphsColors, GraphMarker):
    df=pd.read_csv(infile,comment='#',sep='\t')
    print(df.columns)
    df['Jitter_analitical']=df['N']/df['dV/dt']*1000 # in ps
    df['err_Jitter_analitical']=df['Jitter_analitical']*np.sqrt((df['err_N']/df['N'])**2+(df['err_dV/dt']/df['dV/dt'])**2) # in ps
    df['Gain']=dfGain['Gain']
    df['err_Gain']=[0]*len(dfGain['Gain_err'])
    # FIXME: Gain error is not propagated to the jitter
    gJittervsGainAnalitical=TGraphErrors(len(df['Gain']), np.asarray(df['Gain'], dtype=float), np.asarray(df['Jitter_analitical'], dtype=float), np.asarray([0]*len(df), dtype=float), np.asarray(df['err_Jitter_analitical'], dtype=float))
    gJittervsGainAnalitical.SetName("gJittervsGainAnalitical"+DetectorName[0])
    gJittervsGainAnalitical.SetTitle(";Gain; Jitter (ps)")
    gJittervsGainAnalitical.SetMarkerStyle(GraphMarker[0])
    gJittervsGainAnalitical.SetMarkerSize(1)
    gJittervsGainAnalitical.SetMarkerColor(GraphsColors[0])
    gJittervsGainAnalitical.SetLineColor(GraphsColors[0])
    gJittervsGainMeas = GetGraphErrorsFromPandas(df[['Gain','err_Gain','std_dev','err_dev']])
    gJittervsGainMeas.SetName("gJittervsGainMeas"+DetectorName[1])
    gJittervsGainMeas.SetTitle(";Gain; Jitter (ps)")
    gJittervsGainMeas.SetMarkerStyle(GraphMarker[1])
    gJittervsGainMeas.SetMarkerSize(1)
    gJittervsGainMeas.SetMarkerColor(GraphsColors[1])
    gJittervsGainMeas.SetLineColor(GraphsColors[1])
    return gJittervsGainAnalitical, gJittervsGainMeas


if __name__=='__main__':
    
    gROOT.SetBatch()
    SetGlobalStyle(padleftmargin=0.12, padbottommargin=0.12, padrightmargin=0.05, padtopmargin=0.1, titleoffsety=1.2, titleoffsetx=0.9, titleoffset= 0.7, opttitle=1)

    infiles = ['TCT/data/input/Jitter_vs_gain.csv']
    gainfile = 'TCT/data/output/Gain_vs_Bias_with_gain.csv'
    DetectorNames = [('Measured','Analitical')]
    GraphColors = [(kRed+1,kAzure+3)]
    GraphMarkers = [(kFullCircle, kOpenDiamond)]
    outfilename = 'TCT/data/output/Jitter_vs_gain.root'
    legendmax = 0.85
    legendstep = 0.05
    
    canvas = TCanvas("canvas","canvas",1000,1000)
    hFrame = canvas.cd().DrawFrame(10,0,500,150,";Gain; Jitter (ps)")

    Graphs = []
    legend = TLegend(0.6, legendmax-len(infiles)*legendstep, 0.8, legendmax)
    legend.SetTextSize(0.03)
    for file, name, color, marker in zip(infiles, DetectorNames, GraphColors, GraphMarkers):
        gJittervsGainAnalitical, gJittervsGainMeas = ProduceGraphJittervsGain(file, name, pd.read_csv(gainfile, sep=',', comment='#'), color, marker)
        Graphs.append(gJittervsGainAnalitical)
        Graphs.append(gJittervsGainMeas)
        Graphs[-1].Draw("p,same")
        Graphs[-2].Draw("p,same")
        legend.AddEntry(Graphs[-1],name[0],'p')
        legend.AddEntry(Graphs[-2],name[1],'p')

    legend.Draw("same")

    canvas.Modified()
    canvas.Update()
    
    canvas.SaveAs('TCT/data/output/Jitter_vs_gain.pdf')

    outfile = TFile(outfilename, 'recreate')
    for graph in Graphs:
        graph.Write()
    canvas.Write()
    outfile.Close()
