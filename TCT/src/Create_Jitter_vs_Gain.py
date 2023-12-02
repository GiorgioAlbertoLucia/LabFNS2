import pandas as pd
import numpy as np

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure, kBlack, kFullCircle,kOpenDiamond, kFullSquare, gROOT
from DfUtils import GetGraphErrorsFromCSV, GetGraphErrorsFromPandas
from StyleFormatter import SetObjectStyle, SetGlobalStyle

def ProduceGraphJittervsGain(infile, DetectorName, dfGain, GraphsColors, GraphMarker):
    df=pd.read_csv(infile,comment='#',sep='\t')
    print(df.columns)
    df['Jitter_analitical']=df['N']/df['dV/dt']*1000 # in ps
    df['err_Jitter_analitical']=df['Jitter_analitical']*np.sqrt((df['err_N']/df['N'])**2+(df['err_dV/dt']/df['dV/dt'])**2) # in ps
    df['Gain']=dfGain['Gain']
    df['err_Gain']=dfGain['Gain_err']
    df["err_dev"]=df["std_dev"]/np.sqrt(2*(df["N"]-1))
    gJittervsGainAnalitical=TGraphErrors(len(df['Gain']), np.asarray(df['Gain'], dtype=float), np.asarray(df['Jitter_analitical'], dtype=float), np.asarray(df['err_Gain'], dtype=float), np.asarray(df['err_Jitter_analitical'], dtype=float))
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

    gNoise=TGraphErrors(len(df['N']), np.asarray(df['V'], dtype=float),np.asarray(df['N']), np.ones(len(df['V'])), np.asarray(df['err_N'], dtype=float))
    gNoise.SetName("gNoise")
    gNoise.SetTitle("Noise; Reverse Bias (V); Noise")
    gNoise.SetMarkerStyle(kFullCircle)
    gNoise.SetMarkerSize(3)
    gNoise.SetMarkerColor(kBlack)
    gNoise.SetLineColor(kBlack)
    gNoise.GetYaxis().SetRangeUser(0, 30)
    gNoise.GetXaxis().SetRangeUser(0, 300)

    return gJittervsGainAnalitical, gJittervsGainMeas, gNoise


if __name__=='__main__':
    
    gROOT.SetBatch()
    SetGlobalStyle(padleftmargin=0.12, padbottommargin=0.12, padrightmargin=0.05, padtopmargin=0.1, titleoffsety=1.2, titleoffsetx=0.9, titleoffset= 0.7, opttitle=1)

    infiles = ['TCT/data/input/Jitter_vs_gain.csv']
    gainfile = 'TCT/data/output/Gain_vs_Bias_with_gain.csv'
    DetectorNames = [('Measured','Analitical')]
    GraphColors = [(kRed+1,kAzure+3)]
    GraphMarkers = [(kFullCircle, kOpenDiamond)]
    outfilename = 'TCT/data/output/Jitter_vs_gain.root'
    legendmax = 0.75
    legendstep = 0.07
    outputfile2 = 'TCT/data/output/Noise_vs_Bias.root'
    
    canvas = TCanvas("canvas","canvas",1000,1000)
    hFrame = canvas.cd().DrawFrame(10,0,500,300,"Jitter vs Gain for LGAD sensor;Gain; Jitter (ps)")

    Graphs = []
    legend = TLegend(0.5, legendmax-len(infiles)*legendstep, 0.7, legendmax)
    legend.SetTextSize(0.04)
    for file, name, color, marker in zip(infiles, DetectorNames, GraphColors, GraphMarkers):
        gJittervsGainAnalitical, gJittervsGainMeas, gNoise = ProduceGraphJittervsGain(file, name, pd.read_csv(gainfile, sep=',', comment='#'), color, marker)
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

    canvas1 = TCanvas("canvas1","canvas1",1000,1000)
    hFrame = canvas1.cd().DrawFrame(120,10,280,25,"Noise for LGAD sensor;Bias (V); Noise (mV)")
    gNoise.Draw("p,same")
    legend1 = TLegend(0.6, 0.6, 0.8, 0.8)
    legend1.AddEntry(gNoise,'Noise','p')
    legend1.SetTextSize(0.045)
    legend1.Draw("same")
    canvas1.SaveAs('TCT/data/output/Noise_vs_gain.pdf')
    outfile = TFile(outfilename, 'recreate')
    for graph in Graphs:
        graph.Write()
    canvas.Write()
    outfile.Close()
