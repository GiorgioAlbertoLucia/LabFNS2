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

    infilePin5 = 'TCT/data/input/Gain_vs_Bias_Pin_5MIP.csv'
    infileLGAD5 = 'TCT/data/input/Gain_vs_Bias_LGAD_5MIP.csv'
    infilePin = 'TCT/data/input/Gain_vs_Bias_Pin.csv'
    infileLGAD = 'TCT/data/input/Gain_vs_Bias_LGAD.csv'
    outfilename = 'TCT/data/output/Gain_vs_Bias.root'

    dfPiN = pd.read_csv(infilePin,comment='#')
    dfLGAD = pd.read_csv(infileLGAD,comment='#')
    dfPiN5 = pd.read_csv(infilePin5,comment='#')
    dfLGAD5 = pd.read_csv(infileLGAD5,comment='#')

    dfPiN['QPin']=(dfPiN['AtotPin']-dfPiN['AbasPin'])
    dfPiN['QPin_err']=(dfPiN['Atot_errPin']+dfPiN['Abas_errPin'])
    dfLGAD['QLGAD']=(dfLGAD['AtotLGAD']-dfLGAD['AbasLGAD'])
    dfLGAD['QLGAD_err']=(dfLGAD['Atot_errLGAD']+dfLGAD['Abas_errLGAD'])

    dfPiN5['QPin']=(dfPiN5['AtotPin']-dfPiN5['AbasPin'])
    dfPiN5['QPin_err']=(dfPiN5['Atot_errPin']+dfPiN5['Abas_errPin'])
    dfLGAD5['QLGAD']=(dfLGAD5['AtotLGAD']-dfLGAD5['AbasLGAD'])
    dfLGAD5['QLGAD_err']=(dfLGAD5['Atot_errLGAD']+dfLGAD5['Abas_errLGAD'])

    df = pd.DataFrame()
    df['Gain']=dfLGAD['QLGAD']/dfPiN['QPin']
    df['Gain_err']=np.sqrt((dfLGAD['QLGAD']/(dfPiN['QPin']*dfPiN['QPin'])*dfPiN['QPin_err'])**2+ (dfLGAD['QLGAD_err']/dfPiN['QPin'])**2)
    df['Bias'] = dfPiN['Bias']
    df['Bias_err'] = dfPiN['Bias_err']
    #commento a caso p
    df.to_csv('TCT/data/output/Gain_vs_Bias_with_gain.csv', index=False)

    df5 = pd.DataFrame()
    df5['Gain']=dfLGAD5['QLGAD']/dfPiN5['QPin']
    df5['Gain_err']=np.sqrt((dfLGAD5['QLGAD']/(dfPiN5['QPin']*dfPiN5['QPin'])*dfPiN5['QPin_err'])**2+ (dfLGAD5['QLGAD_err']/dfPiN5['QPin'])**2)
    df5['Bias'] = dfPiN5['Bias']
    df5['Bias_err'] = dfPiN5['Bias_err']

    df5.to_csv('TCT/data/output/Gain_vs_Bias_with_gain_5MIP.csv', index=False)

    gGain=TGraphErrors(len(df['Gain']), np.asarray(df['Bias'], dtype=float), np.asarray(df['Gain'], dtype=float), np.asarray(df['Bias_err'], dtype=float), np.asarray(df['Gain_err'], dtype=float),)
    gGain.SetName("gGain")
    gGain.SetTitle("Gain; Reverse Bias (V); Gain")
    gGain.SetMarkerStyle(kFullSquare)
    gGain.SetMarkerSize(1)
    gGain.SetMarkerColor(kRed+1)
    gGain.SetLineColor(kRed+1)

    gGain5=TGraphErrors(len(df5['Gain']), np.asarray(df5['Bias'], dtype=float), np.asarray(df5['Gain'], dtype=float), np.asarray(df5['Bias_err'], dtype=float), np.asarray(df5['Gain_err'], dtype=float),)
    gGain5.SetName("gGain5")
    gGain5.SetTitle("Gain; Reverse Bias (V); Gain")
    gGain5.SetMarkerStyle(kFullSquare)
    gGain5.SetMarkerSize(1)
    gGain5.SetMarkerColor(kAzure)
    gGain5.SetLineColor(kAzure)
    canvas = TCanvas("canvas","canvas",2100,1000)
    canvas.cd()
    hFrame = canvas.cd().DrawFrame(0,0,500,4000,"Gain vs Bias")
    gGain.Draw("same,p")
    legend1 = TLegend(0.25, 0.5, 0.55, 0.83)
    legend1.AddEntry(gGain,'Gain:LGAD/PiN 1 MIP','p')
    legend1.AddEntry(gGain5,'Gain:LGAD/PiN 5 MIP','p')
    legend1.SetTextSize(0.045)
    legend1.Draw("same")
    

    outfile = TFile(outfilename, 'recreate')
    gGain.Write()
    gGain5.Write()
    canvas.Write()
    canvas.SaveAs('TCT/data/output/Gain_vs_Bias.pdf')

    
    outfile.Close()
