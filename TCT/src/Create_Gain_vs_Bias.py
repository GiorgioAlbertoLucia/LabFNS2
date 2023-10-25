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

    infile = 'TCT/data/input/Gain_vs_Bias_5MIP.csv'
    infilePin = 'TCT/data/input/Gain_vs_Bias_Pin_5MIP.csv'
    infileLGAD = 'TCT/data/input/Gain_vs_Bias_LGAD_5MIP.csv'
    outfilename = 'TCT/data/output/Gain_vs_Bias_5MIP.root'

    #if is necessary only one data file
    dfPiN = pd.read_csv(infilePin,comment='#')
    dfLGAD = pd.read_csv(infileLGAD,comment='#')

    dfPiN['QPin']=(dfPiN['AtotPin']-dfPiN['AbasPin'])
    dfPiN['QPin_err']=(dfPiN['Atot_errPin']+dfPiN['Abas_errPin'])
    dfLGAD['QLGAD']=(dfLGAD['AtotLGAD']-dfLGAD['AbasLGAD'])
    dfLGAD['QLGAD_err']=(dfLGAD['Atot_errLGAD']+dfLGAD['Abas_errLGAD'])

    df = pd.DataFrame()
    df['Gain']=dfLGAD['QLGAD']/dfPiN['QPin']
    df['Gain_err']=np.sqrt((dfLGAD['QLGAD']/(dfPiN['QPin']*dfPiN['QPin'])*dfPiN['QPin_err'])**2+ (dfLGAD['QLGAD_err']/dfPiN['QPin'])**2)
    df['Bias'] = dfPiN['Bias']
    df['Bias_err'] = dfPiN['Bias_err']

    df.to_csv('TCT/data/output/Gain_vs_Bias_with_gain_5MIP.csv', index=False)

    gGain=TGraphErrors(len(df['Gain']), np.asarray(df['Bias'], dtype=float), np.asarray(df['Gain'], dtype=float), np.asarray(df['Bias_err'], dtype=float), np.asarray(df['Gain_err'], dtype=float),)
    gGain.SetName("gGain")
    gGain.SetTitle("Gain; Reverse Bias (V); Gain")
    gGain.SetMarkerStyle(kFullSquare)
    gGain.SetMarkerSize(1)
    gGain.SetMarkerColor(kRed+1)
    gGain.SetLineColor(kRed+1)

    #if are necessary two different data files
    #dfPin=pd.read_csv(infilePin,comment='#')
    #dfLGAD=pd.read_csv(infileLGAD,comment='#')
    #dfd=pd.DataFrame()
    #dfd['QPin']=(dfPin['AtotPin']-dfPin['AbasPin'])
    #dfd['QPin_err']=(dfPin['Atot_errPin']+dfPin['Abas_errPin'])
    #dfd['QLGAD']=(dfLGAD['AtotLGAD']-dfLGAD['AbasLGAD'])
    #dfd['QLGAD_err']=(dfLGAD['Atot_errLGAD']+dfLGAD['Abas_errLGAD'])
    #dfd['Gain']=dfd['QLGAD']/dfd['QPin']
    #dfd['Gain_err']=np.sqrt((dfd['QLGAD']/(dfd['QPin']*dfd['QPin'])*dfd['QPin_err'])**2+ (dfd['QLGAD_err']/dfd['QPin'])**2)
#
    #gGain2=TGraphErrors(len(dfd['Gain']), np.asarray(*dfd['Gain'], dtype=float), np.asarray(dfPin['Bias'], dtype=float), np.asarray(dfPin['Bias_err'], dtype=float))
    #gGain2.SetName("gGain2")
    #gGain2.SetTitle("Gain; Gain; Reverse Bias (V)")
    #gGain2.SetMarkerStyle(kFullSquare)
    #gGain2.SetMarkerSize(1)
    #gGain2.SetMarkerColor(kRed+1)
    #gGain2.SetLineColor(kRed+1)
#
    #canvas= TCanvas()
    #hFrame = canvas.DrawFrame(0,0,100,500,"C-V curves of FBK-UFSD2 ; Gain; Reverse Bias Voltage (V)")
    #gGain.Draw("same,p")
    ##gGain2.Draw("same,p")
    #legend2 = TLegend(0.25, 0.5, 0.55, 0.83)
    #legend2.AddEntry(gGain,'Gain Pin/LGAD','p')
    ##legend2.AddEntry(gGain2,'Gain, -1 * V, f = 10 kHz','p')
    #legend2.SetTextSize(0.045)
    #legend2.Draw("same")
    #canvas.Modified()
    #canvas.Update()
    #canvas.SaveAs('TCT/data/output/Gain_vs_Bias.pdf')

    outfile = TFile(outfilename, 'recreate')
    gGain.Write()
    #gGain2.Write()
    #canvas.Write()
    

    
    outfile.Close()
