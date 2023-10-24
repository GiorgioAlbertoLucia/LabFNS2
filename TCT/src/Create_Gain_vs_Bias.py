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

    infile = 'TCT/data/input/Gain_vs_Bias.csv'
    infilePin = 'TCT/data/input/Gain_vs_Bias_Pin.csv'
    infileLGAD = 'TCT/data/input/Gain_vs_Bias_LGAD.csv'
    outfilename = 'TCT/data/output/Gain_vs_Bias.root'

    #if is necessary only one data file
    df=pd.read_csv(infile,comment='#')
    df['QPin']=(df['AtotPin']-df['AbasPin'])
    df['QPin_err']=(df['Atot_errPin']+df['Abas_errPin'])
    df['QLGAD']=(df['AtotLGAD']-df['AbasLGAD'])
    df['QLGAD_err']=(df['Atot_errLGAD']+df['Abas_errLGAD'])
    df['Gain']=df['QLGAD']/df['QPin']
    df['Gain_err']=np.sqrt((df['QLGAD']/(df['QPin']*df['QPin'])*df['QPin_err'])**2+ (df['QLGAD_err']/df['QPin'])**2)

    gGain=TGraphErrors(len(df['Gain']), np.asarray(*df['Gain'], dtype=float), np.asarray(df['Gain'], dtype=float), np.asarray(df['RevB'], dtype=float), np.asarray(df['RevB_err'], dtype=float))
    gGain.SetName("gGain")
    gGain.SetTitle("Gain; Gain; Reverse Bias (V)")
    gGain.SetMarkerStyle(kFullSquare)
    gGain.SetMarkerSize(1)
    gGain.SetMarkerColor(kRed+1)
    gGain.SetLineColor(kRed+1)

    #if are necessary two different data files
    dfPin=pd.read_csv(infilePin,comment='#')
    dfLGAD=pd.read_csv(infileLGAD,comment='#')
    dfd=pd.DataFrame()
    dfd['QPin']=(dfPin['AtotPin']-dfPin['AbasPin'])
    dfd['QPin_err']=(dfPin['Atot_errPin']+dfPin['Abas_errPin'])
    dfd['QLGAD']=(dfLGAD['AtotLGAD']-dfLGAD['AbasLGAD'])
    dfd['QLGAD_err']=(dfLGAD['Atot_errLGAD']+dfLGAD['Abas_errLGAD'])
    dfd['Gain']=dfd['QLGAD']/dfd['QPin']
    dfd['Gain_err']=np.sqrt((dfd['QLGAD']/(dfd['QPin']*dfd['QPin'])*dfd['QPin_err'])**2+ (dfd['QLGAD_err']/dfd['QPin'])**2)

    gGain2=TGraphErrors(len(dfd['Gain']), np.asarray(*dfd['Gain'], dtype=float), np.asarray(dfPin['Bias'], dtype=float), np.asarray(dfPin['Bias_err'], dtype=float))
    gGain2.SetName("gGain2")
    gGain2.SetTitle("Gain; Gain; Reverse Bias (V)")
    gGain2.SetMarkerStyle(kFullSquare)
    gGain2.SetMarkerSize(1)
    gGain2.SetMarkerColor(kRed+1)
    gGain2.SetLineColor(kRed+1)

    canvas= TCanvas()
    hFrame = canvas.DrawFrame(0,0,100,500,"C-V curves of FBK-UFSD2 ; Gain; Reverse Bias Voltage (V)")
    gGain.Draw("same,p")
    #gGain2.Draw("same,p")
    legend2 = TLegend(0.25, 0.5, 0.55, 0.83)
    legend2.AddEntry(gGain,'Gain Pin/LGAD','p')
    #legend2.AddEntry(gGain2,'Gain, -1 * V, f = 10 kHz','p')
    legend2.SetTextSize(0.045)
    legend2.Draw("same")
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('TCT/data/output/Gain_vs_Bias.pdf')

    outfile = TFile(outfilename, 'recreate')
    gGain.Write()
    #gGain2.Write()
    canvas.Write()
    

    
    outfile.Close()
