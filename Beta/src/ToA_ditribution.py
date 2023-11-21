import pandas as pd
import numpy as np
import uproot

import sys
sys.path.append('utils')

from ROOT import TGraphErrors, TFile, TCanvas, TLegend, kRed, kAzure,  kOrange, kFullCircle, kFullSquare, gROOT, TF1, TGaxis, TTree, TH1D, TH2D, gStyle
from DfUtils import GetGraphErrorsFromCSV
from StyleFormatter import SetObjectStyle, SetGlobalStyle
import sys
import itertools

def ComputeSystOnThresholds(df, th2, th3, Up=False):
    cuts = itertools.product(th2, th3)
    resolutions = []
    resolutionsfit = []

    hResvsCut = TH1D("hResvsCut", "hResvsCut", len(th2)*len(th3), 0, len(th2)*len(th3))
    hResFitvsCut = TH1D("hResFitvsCut", "hResFitvsCut", len(th2)*len(th3), 0, len(th2)*len(th3))
    hResvsTh = TH2D("hRes", "hRes", len(th2), min(th2), max(th2), len(th3), min(th3), max(th3))
    hResFitvsTh = TH2D("hResFit", "hResFit", len(th2), min(th2), max(th2), len(th3), min(th3), max(th3))
    hResDistr = TH1D("hResDistr", "hResDistr", 100, 0, 100)
    hResFitDistr = TH1D("hResFitDistr", "hResFitDistr", 100, 0, 100)
    for idx, cut in enumerate(cuts):
        th2 = cut[0]
        th3 = cut[1]
        if(Up==True):
            _, _, _, _, timeres, timeresfit = ComputeResolutionUp(df, th2, th3)
        else:   
            _, _, _, _, timeres, timeresfit = ComputeResolution(df, th2, th3)
        timeres = timeres/np.sqrt(2)
        timeresfit = timeresfit/np.sqrt(2)
        resolutions.append(timeres)
        resolutionsfit.append(timeresfit)
        hResvsTh.SetBinContent(hResvsTh.FindBin(th2, th3), timeres)
        hResFitvsTh.SetBinContent(hResFitvsTh.FindBin(th2, th3), timeresfit)
        hResvsCut.SetBinContent(idx, timeres)
        hResFitvsCut.SetBinContent(idx, timeresfit)
        hResDistr.Fill(timeres)
        hResFitDistr.Fill(timeresfit)
        
    return resolutions, resolutionsfit, hResvsCut, hResFitvsCut, hResvsTh, hResFitvsTh, hResDistr, hResFitDistr

def ComputeResolution(df, th2, th3):

    hTime=TH1D("hTime","hTime",19,-350,550.)
    hTimeFit=TH1D("hTimeFit","hTimeFit",46,-350,550.)
    for a2,a3,t2,t3,t2fit,t3fit in zip(df["Amp2"],df["Amp3"],df["ToA2"],df["ToA3"], df["ToA2f"], df["ToA3f"]): 
        if (a2 >= th2 and a3>= th3):
            hTime.Fill(1000*(t2-t3))
            hTimeFit.Fill(1000*(t2fit-t3fit))
    gaus = TF1("gaus","[0]*exp(-((x-[1])^2)/(2*[2]^2))",-100,300)
    gaus.SetLineColor(kRed)
    gaus.SetParameter(2,50)
    gaus.SetParameter(1,80)
    gaus.SetParameter(0,3500)
    gStyle.SetOptFit(11111111)
    hTime.Fit(gaus,"qrml+")
    timeres=gaus.GetParameter(2)

    
    gausFit=TF1("gausFit","[0]*exp(-((x-[1])^2)/(2*[2]^2))",0,240)
    gausFit.SetLineColor(kAzure)
    gausFit.SetParameter(2,50)
    gausFit.SetParameter(1,80)
    gausFit.SetParameter(0,1800)
    gStyle.SetOptFit(11111111)
    hTimeFit.Fit(gausFit,"qrml+")
    timeresfit=gausFit.GetParameter(2)

    return hTime, hTimeFit, gaus, gausFit, timeres, timeresfit

def ComputeResolutionUp(df, th2Up, th3Up):
    
        hTime1=TH1D("hTime1","hTime1",19,-350,550.)
        hTimeFit1=TH1D("hTimeFit1","hTimeFit1",46,-350,550.)
        for a2,a3,t2,t3,t2fit,t3fit in zip(df["Amp2"],df["Amp3"],df["ToA2"],df["ToA3"], df["ToA2f"], df["ToA3f"]): 
            if (a2 <= th2Up and a3<= th3Up and a2>80 and a3>80):
                hTime1.Fill(1000*(t2-t3))
                hTimeFit1.Fill(1000*(t2fit-t3fit))
        gaus1 = TF1("gaus","[0]*exp(-((x-[1])^2)/(2*[2]^2))",-100,300)
        gaus1.SetLineColor(kRed)
        gaus1.SetParameter(2,50)
        gaus1.SetParameter(1,80)
        gaus1.SetParameter(0,3500)
        gStyle.SetOptFit(11111111)
        hTime1.Fit(gaus1,"qrml+")
        timeres1=gaus1.GetParameter(2)

        gausFit1=TF1("gausFit","[0]*exp(-((x-[1])^2)/(2*[2]^2))",0,240)
        gausFit1.SetLineColor(kAzure)
        gausFit1.SetParameter(2,50)
        gausFit1.SetParameter(1,80)
        gausFit1.SetParameter(0,1800)
        gStyle.SetOptFit(11111111)
        hTimeFit1.Fit(gausFit1,"qrml+")
        timeresfit1=gausFit1.GetParameter(2)

        return hTime1, hTimeFit1, gaus1, gausFit1, timeres1, timeresfit1
   

if __name__=='__main__':
    gROOT.SetBatch(True)
    tree=uproot.open("Beta/data/output/BetaOutput.root")["BetaTree"]
    outfilename = 'Beta/data/output/Time_resolution.root'
    df=tree.arrays(library='pd')
    theshold2=80
    theshold3=80
    theshold2Syst = np.arange(0, 150, 10)
    theshold3Syst = np.arange(30, 150,10)
    theshold2SystUp = np.arange(240, 360, 10)
    theshold3SystUp = np.arange(250, 360, 10)
    #theshold2Syst = np.flip(theshold2Syst)
    #theshold3Syst = np.flip(theshold3Syst)
    _, _, hResvsCut, hResFitvsCut, hResvsTh, hResFitvsTh, hResDistr, hResFitDistr  = ComputeSystOnThresholds(df, theshold2Syst, theshold3Syst)
    _, _, hResvsCutUp, hResFitvsCutUp, hResvsThUp, hResFitvsThUp, hResDistrUp, hResFitDistrUp  = ComputeSystOnThresholds(df, theshold2SystUp, theshold3SystUp, Up=True)
    hTime, hTimeFit, gaus, gaus2, timeres, timeresfit = ComputeResolution(df, theshold2, theshold3)
    print("risolution UFSD: ",timeres/np.sqrt(2))
    print("risolution UFSD with fit: ",timeresfit/np.sqrt(2))
    canvas=TCanvas("canvas","canvas",1900,1500)
    canvas.Divide(2,1)
    canvas.cd(1)
    hTime.SetTitle("#Delta_{t} distribution, theshold:80 mV ")
    hTime.GetXaxis().SetTitle("#Delta_{t} (ps)")
    hTime.GetYaxis().SetTitle("Counts")
    hTime.Draw("hist")
    gaus.Draw("same")
    legend1 = TLegend(0.6, 0.3, 0.85, 0.5)
    legend1.AddEntry(hTime,'Data','l')
    legend1.AddEntry(gaus,'p0*exp(#frac{-(x-p1)^{2}}{2*p2^{2}})','l')
    legend1.SetTextSize(0.045)
    legend1.Draw("same")
    print("chi/dof : ",gaus.GetChisquare(),"/",gaus.GetNDF())
    canvas.cd(2)
    hTimeFit.SetTitle("#Delta_{t} distribution, theshold:80 mV ")
    hTimeFit.GetXaxis().SetTitle("#Delta_{t} (ps)")
    hTimeFit.GetYaxis().SetTitle("Counts")
    hTimeFit.Draw("hist")
    gaus2.Draw("same")
    legend2 = TLegend(0.6, 0.3, 0.85, 0.5)
    legend2.AddEntry(hTimeFit,'Data','l')
    legend2.AddEntry(gaus2,'p0*exp(#frac{-(x-p1)^{2}}{2*p2^{2}})','l')
    legend2.SetTextSize(0.045)
    legend2.Draw("same")
    print("chi/dof fit: ",gaus2.GetChisquare(),"/",gaus2.GetNDF())
    outfile = TFile(outfilename, 'recreate')
    hTimeFit.Write()
    hTime.Write()
    gaus.Write()
    gaus2.Write()
    canvas.Write()
    canvas.SaveAs('Beta/data/output/Time_distribution.pdf')

    canvasSys=TCanvas("canvasSys","canvasSys",1900,1500)
    canvasSys.Divide(2,2)
    canvasSys.cd(1).DrawFrame(0, 0, len(theshold2Syst)*len(theshold3Syst), 130, "Resolution vs cut;Cut;Resolution")
    hResvsCut.SetLineColor(kRed)
    hResvsCut.Draw("hist,same")
    hResFitvsCut.SetLineColor(kAzure)
    hResFitvsCut.Draw("hist,same")
    legendResvsCut = TLegend(0.6, 0.7, 0.85, 0.9)
    legendResvsCut.AddEntry(hResvsCut,'w/o fit','l')
    legendResvsCut.AddEntry(hResFitvsCut,'w fit','l')
    legendResvsCut.SetTextSize(0.045)
    legendResvsCut.Draw("same")
    canvasSys.cd(2)
    hResvsTh.SetTitle("Resolution vs threshold")
    hResvsTh.GetXaxis().SetTitle("Threshold 2")
    hResvsTh.GetYaxis().SetTitle("Threshold 3")
    hResvsTh.SetStats(0)
    hResvsTh.SetAxisRange(50.,60.,"z")
    hResvsTh.Draw("colz")
    canvasSys.cd(3)
    hResFitvsTh.SetTitle("Resolution (fit) vs threshold")
    hResFitvsTh.GetXaxis().SetTitle("Threshold 2")
    hResFitvsTh.GetYaxis().SetTitle("Threshold 3")
    hResFitvsTh.SetStats(0)
    hResFitvsTh.SetAxisRange(35., 45.,"z")
    hResFitvsTh.Draw("colz")
    canvasSys.cd(4).DrawFrame(0, 0, 100, 150, "Resolution distribution;Resolution;Counts")
    hResDistr.SetLineColor(kRed)
    hResDistr.Draw("hist,same")
    hResFitDistr.SetLineColor(kAzure)
    hResFitDistr.Draw("hist,same")
    legendResDistr = TLegend(0.6, 0.7, 0.85, 0.9)
    legendResDistr.AddEntry(hResDistr,'w/o fit','l')
    legendResDistr.AddEntry(hResFitDistr,'w fit','l')
    legendResDistr.SetTextSize(0.045)
    legendResDistr.Draw("same")
    canvasSys.Write()
    canvasSys.SaveAs('Beta/data/output/Time_resolution.pdf')

    canvasSysUp=TCanvas("canvasSysUp","canvasSysUp",1900,1500)
    canvasSysUp.Divide(2,2)
    canvasSysUp.cd(1).DrawFrame(0, 0,len(theshold2SystUp)*len(theshold3SystUp), 130, "Resolution vs cut;Cut;Resolution")
    hResvsCutUp.SetLineColor(kRed)
    hResvsCutUp.Draw("hist,same")
    hResFitvsCutUp.SetLineColor(kAzure)
    hResFitvsCutUp.Draw("hist,same")
    legendResvsCutUp = TLegend(0.6, 0.7, 0.85, 0.9)
    legendResvsCutUp.AddEntry(hResvsCutUp,'w/o fit','l')
    legendResvsCutUp.AddEntry(hResFitvsCutUp,'w fit','l')
    legendResvsCutUp.SetTextSize(0.045)
    legendResvsCutUp.Draw("same")
    canvasSysUp.cd(2)
    hResvsThUp.SetTitle("Resolution vs threshold")
    hResvsThUp.GetXaxis().SetTitle("High Threshold 2")
    hResvsThUp.GetYaxis().SetTitle("High Threshold 3")
    hResvsThUp.SetStats(0)
    hResvsThUp.SetAxisRange(46.,58.,"z")
    hResvsThUp.Draw("colz")
    canvasSysUp.cd(3)
    hResFitvsThUp.SetTitle("Resolution (fit) vs threshold")
    hResFitvsThUp.GetXaxis().SetTitle("High Threshold 2")
    hResFitvsThUp.GetYaxis().SetTitle("High Threshold 3")
    hResFitvsThUp.SetStats(0)
    hResFitvsThUp.SetAxisRange(35., 42.,"z")
    hResFitvsThUp.Draw("colz")
    canvasSysUp.cd(4).DrawFrame(20, 0, 80, 150, "Resolution distribution;Resolution;Counts")
    hResDistrUp.SetLineColor(kRed)
    hResDistrUp.Draw("hist,same")
    hResFitDistrUp.SetLineColor(kAzure)
    hResFitDistrUp.Draw("hist,same")
    legendResDistrUp = TLegend(0.6, 0.7, 0.85, 0.9)
    legendResDistrUp.AddEntry(hResDistrUp,'w/o fit','l')
    legendResDistrUp.AddEntry(hResFitDistrUp,'w fit','l')
    legendResDistrUp.SetTextSize(0.045)
    legendResDistrUp.Draw("same")
    canvasSysUp.Write()
    canvasSysUp.SaveAs('Beta/data/output/Time_resolutionUp.pdf')




    outfile.Close()