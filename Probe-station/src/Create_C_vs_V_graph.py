import pandas as pd
import numpy as np

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure,  kOrange, kFullCircle, kFullSquare, gROOT, TF1, TGaxis
from DfUtils import GetGraphErrorsFromCSV
from StyleFormatter import SetObjectStyle, SetGlobalStyle


if __name__=='__main__':
    
    gROOT.SetBatch()
    SetGlobalStyle(padleftmargin=0.12, padbottommargin=0.12, padrightmargin=0.05, padtopmargin=0.1, titleoffsety=1.2, titleoffsetx=0.9, titleoffset= 0.7, opttitle=1)

    infileLGAD = 'Probe-station/data/input/C_vs_V_LGAD.csv'
    infilePin = 'Probe-station/data/input/C_vs_V_pin.csv'
    infileStrip = 'Probe-station/data/input/C_vs_V_strip.csv'
    outfilename = 'Probe-station/data/output/C_vs_V.root'

    dfS=pd.read_csv(infileStrip)
    dfP=pd.read_csv(infilePin)
    dfL=pd.read_csv(infileLGAD)
    '''gLGAD=TGraphErrors(len(dfL['V']), np.asarray(-4.2*dfL['V'], dtype=float), np.asarray(dfL['C'], dtype=float), np.asarray(dfL['V_err'], dtype=float), np.asarray(dfL['C_err'], dtype=float))
    gPin=TGraphErrors(len(dfP['V']), np.asarray(-21*dfP['V'], dtype=float), np.asarray(9.2857*dfP['C'], dtype=float), np.asarray(dfP['V_err'], dtype=float), np.asarray(dfP['C_err'], dtype=float))

    #gLGAD = GetGraphErrorsFromCSV(infileLGAD) #change voltage sign in csv (from - to +), but the correct values are negative
    gLGAD.SetName("gLGAD")
    gLGAD.SetTitle("LGAD; Reverse Bias Voltage (V); Capacitance Strip and LGDA (pF)")
    gLGAD.SetMarkerStyle(kFullCircle)
    gLGAD.SetMarkerSize(1)
    gLGAD.SetMarkerColor(kRed+1)
    gLGAD.SetLineColor(kRed+1)
    gLGAD.GetXaxis().SetRangeUser(-55., 10)
    
    #gPin = GetGraphErrorsFromCSV(infilePin)
    gPin.SetName("gPin")
    gPin.SetTitle("Pin; Reverse Bias Voltage (V); Capacitance Strip and LGDA (pF)")
    gPin.SetMarkerStyle(kFullSquare)
    gPin.SetMarkerSize(1)
    gPin.SetMarkerColor(kAzure + 3)
    gPin.SetLineColor(kAzure + 3)
    gPin.GetXaxis().SetRangeUser(-55., 10)

    gStrip = GetGraphErrorsFromCSV(infileStrip)
    gStrip.SetName("gStrip")
    gStrip.SetTitle("Strip; Reverse Bias Voltage Strip (V); Capacitance (pF)")
    gStrip.SetMarkerStyle(kFullSquare)
    gStrip.SetMarkerSize(1)
    gStrip.SetMarkerColor(kOrange - 3)
    gStrip.SetLineColor(kOrange - 3)
    
    canvas = TCanvas("canvas","canvas",1000,1000)
    canvas.cd()
    hFrame = canvas.cd().DrawFrame(-11.5,0,210,650,"C-V curves of FBK-UFSD2 LGAD and PiN pads, strip Hamamatsu baby 2x0-4; Reverse Bias Voltage Strip (V); Strip and LGDA (pF)")
    #canvas.SetLogx()
    gLGAD.Draw("p,same")
    gPin.Draw("p,same")
    gStrip.Draw("p,same")
    #from here to line 83 usefull to set different scale (is necessary modify the value in the graph)
    fL = TF1("fL","3*x",0,50);
    A1 = TGaxis(0,650,210,650,"fL",510,"-");
    A1.SetLineColor(kRed+1)
    A1.SetLabelColor(kRed+1)
    A1.SetTitleColor(kRed+1)
    A1.SetTitle("LGDA Reverse Bias Voltage (V)");
    A1.Draw();
    fP = TF1("fP","3*x",-0.5,10);
    A2 = TGaxis(-11.5,-50,210,-50,"fP",510,"-");
    A2.SetTitle("Pin Reverse Bias Voltage (V)")
    A2.SetLineColor(kAzure + 3)
    A2.SetLineWidth(2)
    A2.SetLabelColor(kAzure + 3)
    A2.SetTitleColor(kAzure + 3)
    A2.Draw()
    fC = TF1("fC","3*x",0,70);
    A3 = TGaxis(0,0,0,650,"fC",510,"-");
    A3.SetTitle("Pin  Capacitance (pF) ")
    A3.SetLineColor(kAzure + 3)
    A3.SetLineWidth(2)
    A3.SetLabelColor(kAzure + 3)
    A3.SetTitleColor(kAzure + 3)
    A3.Draw()'''

    gLGAD=TGraphErrors(len(dfL['V']), np.asarray(-1*dfL['V'], dtype=float), np.asarray(dfL['C'], dtype=float), np.asarray(dfL['V_err'], dtype=float), np.asarray(dfL['C_err'], dtype=float))
    gPin=TGraphErrors(len(dfP['V']), np.asarray(-1*dfP['V'], dtype=float), np.asarray(dfP['C'], dtype=float), np.asarray(dfP['V_err'], dtype=float), np.asarray(dfP['C_err'], dtype=float))

    gLGAD.SetName("gLGAD")
    gLGAD.SetTitle("LGAD; Reverse Bias Voltage (V); Capacitance Strip and LGDA (pF)")
    gLGAD.SetMarkerStyle(kFullSquare)
    gLGAD.SetMarkerSize(1)
    gLGAD.SetMarkerColor(kRed+1)
    gLGAD.SetLineColor(kRed+1)
    
    #gPin = GetGraphErrorsFromCSV(infilePin)
    gPin.SetName("gPin")
    gPin.SetTitle("Pin; Reverse Bias Voltage (V); Capacitance Strip and LGDA (pF)")
    gPin.SetMarkerStyle(kFullSquare)
    gPin.SetMarkerSize(1)
    gPin.SetMarkerColor(kAzure + 3)
    gPin.SetLineColor(kAzure + 3)

    gStrip = GetGraphErrorsFromCSV(infileStrip)
    gStrip.SetName("gStrip")
    gStrip.SetTitle("Strip; Reverse Bias Voltage Strip (V); Capacitance (pF)")
    gStrip.SetMarkerStyle(kFullCircle)
    gStrip.SetMarkerSize(1.5)
    gStrip.SetMarkerColor(kOrange - 3)
    gStrip.SetLineColor(kOrange - 3)
    
    canvas = TCanvas("canvas","canvas",2100,1000)
    canvas.Divide(3,1)
    canvas.cd(3)
    hFrame = canvas.cd(3).DrawFrame(-4,0,50.5,350,"C-V curves of FBK-UFSD2 LGAD pad; Reverse Bias Voltage (V); Capacitance (pF)")
    gLGAD.Draw("p,same")
    legend3 = TLegend(0.25, 0.5, 0.55, 0.83)
    legend3.AddEntry(gLGAD,'LGAD, -1 * V, f = 5kHz, rev. pol.','p')
    legend3.SetTextSize(0.045)
    legend3.Draw("same")
    
    canvas.cd(1)
    hFrame = canvas.cd(1).DrawFrame(-2,0,230,650,"C-V curves of Hamamatsu Stipt ; Reverse Bias Voltage  (V); Capacitance (pF)")
    gStrip.Draw("same,p")
    legend1 = TLegend(0.25, 0.5, 0.55, 0.83)
    legend1.AddEntry(gStrip,'Strip, f = 1 kHz, dir. pol.','p')
    legend1.SetTextSize(0.045)
    legend1.Draw("same")

    canvas.cd(2)
    hFrame = canvas.cd(2).DrawFrame(-1,0,10.5,70,"C-V curves of FBK-UFSD2 PiN pad ; Reverse Bias Voltage  (V); Capacitance (pF)")
    gPin.Draw("same,p")
    legend2 = TLegend(0.25, 0.5, 0.55, 0.83)
    legend2.AddEntry(gPin,'PiN, -1 * V, f = 10 kHz,  rev. pol.','p')
    legend2.SetTextSize(0.045)
    legend2.Draw("same")
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('Probe-station/data/output/C_vs_V.pdf')

    outfile = TFile(outfilename, 'recreate')
    gLGAD.Write()
    gPin.Write()
    canvas.Write()
    

    
    outfile.Close()
