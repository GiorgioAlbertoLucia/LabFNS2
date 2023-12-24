#script to study linearity versus VH for pixel 5 for baseline, falltime and amplitude
import pandas as pd
import numpy as np
import glob
import re

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, kRed, kAzure, kOrange, kFullCircle, gROOT, TLatex, gStyle
from StyleFormatter import SetObjectStyle, SetGlobalStyle

def SetGraph(graph, name, title, color, marker):
    graph.SetName(name)
    graph.SetTitle(title)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.5)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)

if __name__=='__main__':
    dir_path='ITS3/Data/VH'
    txt_files=glob.glob(dir_path+'/*.txt')
    VHvalue=np.array([])
    baseline=np.array([])
    baseline_err=np.array([])
    amplitude=np.array([])
    amplitude_err=np.array([])
    falltime=np.array([])
    falltime_err=np.array([])

    for file_name in txt_files:
        parts=file_name.split('_')
        VH=parts[2]
        match= re.search(r'(\d+\.\d+)(?=V)',VH)
        VHvalue=np.append(VHvalue, float(match.group()))

        with open(file_name,'r') as file:
            for line in file:
                data=line.split()
                baseline=np.append(baseline,float(data[0])*1000)
                baseline_err=np.append(baseline_err,float(data[1])*1000)
                amplitude=np.append(amplitude,float(data[2])*1000)
                amplitude_err=np.append(amplitude_err,float(data[3])*1000)
                falltime=np.append(falltime,float(data[4])*1.e9)
                falltime_err=np.append(falltime_err,float(data[5])*1.e9)
  
    gBaseline=TGraphErrors(len(VHvalue), np.asarray(VHvalue, dtype=float), np.asarray(baseline, dtype=float), np.zeros(len(VHvalue)), np.asarray(baseline_err, dtype=float))
    SetGraph(gBaseline, "gBaseline", ";VH (V); Baseline (mV)", kRed+1, kFullCircle)
    gAmplitude=TGraphErrors(len(VHvalue), np.asarray(VHvalue, dtype=float), np.asarray(amplitude, dtype=float), np.zeros(len(VHvalue)), np.asarray(amplitude_err, dtype=float))
    SetGraph(gAmplitude, "gAmplitude", ";VH (V); Amplitude (mV)", kAzure+3, kFullCircle)
    gFalltime=TGraphErrors(len(VHvalue), np.asarray(VHvalue, dtype=float), np.asarray(falltime, dtype=float), np.zeros(len(VHvalue)), np.asarray(falltime_err, dtype=float))
    SetGraph(gFalltime, "gFalltime", ";VH (V); Falltime (ns)", kOrange+1, kFullCircle)

    canvas = TCanvas("canvas","canvas",2200,600)
    canvas.Divide(3,1)
    canvas.cd(1)
    hFrame = canvas.cd(1).DrawFrame(0,210,1.3,240,"Baseline vs VH; VH (V); Baseline (mV)")
    gBaseline.Draw("p,same")
    text =TLatex(0.15, 0.8,"APTS_AO10P_B6, Pixel 5")
    text.SetNDC()
    text.SetTextSize(gStyle.GetTextSize())
    text.SetTextFont(42)
    text.Draw()

    canvas.cd(2)
    hFrame = canvas.cd(2).DrawFrame(0,0,1.3,50,"Amplitude vs VH; VH (V); Amplitude (mV)")
    gAmplitude.Draw("p,same")
    text2 =TLatex(0.15, 0.80,"APTS_AO10P_B6, Pixel 5")
    text2.SetNDC()
    text2.SetTextSize(gStyle.GetTextSize())
    text2.SetTextFont(42)
    text2.Draw()
    #canvas.SaveAs('ITS3/data/output/Base_Ampli_VH.pdf')

    canvas.cd(3)
    hFrame = canvas.cd(3).DrawFrame(0,0,1.3,7,"Falltime vs VH; VH (V); Falltime (ns)")
    gFalltime.Draw("p,same")
    text3 =TLatex(0.4, 0.80,"APTS_AO10P_B6, Pixel 5")
    text3.SetNDC()
    text3.SetTextSize(gStyle.GetTextSize())
    text3.SetTextFont(42)
    text3.Draw()
    canvas.SaveAs('ITS3/data/output/VH.pdf')