from ROOT import TFile, TH1D, TCanvas, TLegend, TLatex, TF1, gStyle, TGraphErrors, TMultiGraph, gPad
import uproot
import pandas as pd
import yaml

def drawXProjection(inFile):

    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)

    hist = inFile.Get('chip113_projX')
    hist.SetTitle('Source acquisition - x projection; x (pixel); Counts')
    hist.SetFillColorAlpha(797, 0.5)
    hist.SetLineColor(797)
    hist.SetFillStyle(3356)

    canvas = TCanvas('canvas', '^{90}Sr acquisition - x projection ', 800, 600)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.15)
    canvas.SetBottomMargin(0.15)

    hframe = canvas.DrawFrame(-200, 0, 1024, 5000, r'^{90}Sr acquisition - x projection; x (pixel); Counts')
    #hframe = canvas.DrawFrame(-200, 0, 1024, 3500, r'^{90}Sr acquisition - x projection; x (pixel); Counts')

    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetLabelSize(0.04)
    hframe.GetYaxis().SetLabelSize(0.04)


    fit = hist.GetFunction('projX')
    fit.SetLineColor(863)

    legend = TLegend(0.2, 0.32, 0.5, 0.57)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.AddEntry(hist, 'Data', 'f')
    legend.AddEntry(fit, 'Fit', 'l')

    canvas.cd()
    hist.Draw('hist same')
    fit.Draw('same')
    legend.Draw('same')

    latex = TLatex()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatexNDC(0.2, 0.85, '^{90}Sr acquisition')
    latex.DrawLatexNDC(0.2, 0.8, 'Chip 113')
    latex.DrawLatexNDC(0.2, 0.75, f'Threshold:')
    #latex.DrawLatexNDC(0.2, 0.7, f'Configuration 2')
    latex.DrawLatexNDC(0.2, 0.7, f'Configuration 3')
    latex.DrawLatexNDC(0.2, 0.65, f'x position = ({fit.GetParameter(1):.0f} #pm {fit.GetParameter(2):.0f})')
    

    #canvas.SaveAs('../Data/output/sourceXProjection.pdf')
    canvas.SaveAs('../Data/output/sourceXProjection_lowTH.pdf')
    del canvas

def drawYProjection(inFile):

    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)

    hist = inFile.Get('chip113_projY')
    hist.SetTitle('Source acquisition - y projection; y (pixel); Counts')
    hist.SetFillColorAlpha(797, 0.5)
    hist.SetLineColor(797)
    hist.SetFillStyle(3356)

    canvas = TCanvas('canvas', '^{90}Sr acquisition - y projection ', 900, 600)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.15)
    canvas.SetBottomMargin(0.15)
    hframe = canvas.DrawFrame(-300, 0, 512, 1700, r'^{90}Sr acquisition - y projection; y (pixel); Counts')
    #hframe = canvas.DrawFrame(-300, 0, 512, 1200, r'^{90}Sr acquisition - y projection; y (pixel); Counts')

    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetLabelSize(0.04)
    hframe.GetYaxis().SetLabelSize(0.04)

    fit = hist.GetFunction('projY')
    fit.SetLineColor(863)

    legend = TLegend(0.2, 0.27, 0.5, 0.52)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.AddEntry(hist, 'Data', 'f')
    legend.AddEntry(fit, 'Fit', 'l')

    canvas.cd()
    hist.Draw('same hist')
    fit.Draw('same')
    legend.Draw('same')

    latex = TLatex()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatexNDC(0.2, 0.8, '^{90}Sr acquisition')
    latex.DrawLatexNDC(0.2, 0.75, 'Chip 113')
    latex.DrawLatexNDC(0.2, 0.7, f'Threshold:')
    #latex.DrawLatexNDC(0.2, 0.65, f'Configuration 2')
    latex.DrawLatexNDC(0.2, 0.65, f'Configuration 3')
    latex.DrawLatexNDC(0.2, 0.6, f'x range = [505, 575]')
    latex.DrawLatexNDC(0.2, 0.55, f'y position = ({fit.GetParameter(1):.0f} #pm {fit.GetParameter(2):.0f})')
    

    #canvas.SaveAs('../Data/output/sourceYProjection.pdf')    
    canvas.SaveAs('../Data/output/sourceYProjection_lowTH.pdf')
    del canvas

def drawAcquisition(inFile):

    hist = inFile.Get('chip113')
    hist.SetTitle('^{90}Sr acquisition - pixel map; x (pixel); y (pixel)')

    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)

    # eliminate -1
    for xbin in range(1, hist.GetNbinsX()+1):
        for ybin in range(1, hist.GetNbinsY()+1):
            if hist.GetBinContent(xbin, ybin) < 0:    hist.SetBinContent(xbin, ybin, 0)

    gStyle.SetPalette(53)

    canvas = TCanvas('canvas', '^{90}Sr acquisition', 900, 600)
    
    canvas.cd()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.15)
    canvas.SetBottomMargin(0.15)
    #canvas.SetLogz()

    hframe = canvas.DrawFrame(-200, -50, 1024, 562, r'^{90}Sr acquisition - pixel map; x (pixel);y (pixel)')

    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetLabelSize(0.04)
    hframe.GetYaxis().SetLabelSize(0.04)

    hist.Draw('colz same')

    latex = TLatex()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatexNDC(0.2, 0.85, '^{90}Sr acquisition')
    latex.DrawLatexNDC(0.2, 0.75, 'Background subtracted')
    latex.DrawLatexNDC(0.2, 0.8, 'Chip 113')
    #latex.DrawLatexNDC(0.2, 0.17, f'Threshold: Configuration 2')
    latex.DrawLatexNDC(0.2, 0.17, f'Threshold: Configuration 3')
    
    #canvas.SaveAs('../Data/output/source.pdf')
    canvas.SaveAs('../Data/output/source_lowTH.pdf')
    del canvas

def drawHitsAcquisition(inFile):

    hist = inFile.Get('chip113_hits')
    hist.SetTitle('^{90}Sr acquisition; x (pixel); y (pixel)')

    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)

    # eliminate -1
    for xbin in range(1, hist.GetNbinsX()+1):
        for ybin in range(1, hist.GetNbinsY()+1):
            if hist.GetBinContent(xbin, ybin) < 0:    hist.SetBinContent(xbin, ybin, 0)

    gStyle.SetPalette(53)

    canvas = TCanvas('canvas', '^{90}Sr acquisition', 900, 600)
    
    canvas.cd()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.15)
    canvas.SetBottomMargin(0.15)
    #canvas.SetLogz()

    hframe = canvas.DrawFrame(-200, -50, 1024, 562, r'^{90}Sr acquisition; x (pixel);y (pixel)')

    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetLabelSize(0.04)
    hframe.GetYaxis().SetLabelSize(0.04)

    hist.Draw('colz same')

    latex = TLatex()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatexNDC(0.2, 0.85, '^{90}Sr acquisition')
    latex.DrawLatexNDC(0.2, 0.75, 'Background subtracted')
    latex.DrawLatexNDC(0.2, 0.8, 'Chip 113')
    #latex.DrawLatexNDC(0.2, 0.17, f'Threshold: Configuration 2')
    latex.DrawLatexNDC(0.2, 0.17, f'Threshold: Configuration 3')
    
    #canvas.SaveAs('../Data/output/source_hits.pdf')
    canvas.SaveAs('../Data/output/source_hits_lowTH.pdf')
    del canvas

def drawClusterSize(inFile):

    hist = inFile.Get('clusterSize')
    hist.SetTitle('^{90}Sr acquisition; Cluster size; Counts')

    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)

    hist.SetFillColorAlpha(797, 0.5)
    hist.SetLineColor(797)
    hist.SetFillStyle(3356)

    canvas = TCanvas('canvas', '^{90}Sr acquisition', 900, 600)
    
    canvas.cd()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.15)
    canvas.SetBottomMargin(0.15)
    #canvas.SetLogz()
    #canvas.SetLogy()

    hframe = canvas.DrawFrame(-0.5, 100, 10.5, 1e8, r'^{90}Sr acquisition; Cluster size; Counts')
    canvas.SetLogy()

    hframe.GetXaxis().SetTitleSize(0.05)
    hframe.GetYaxis().SetTitleSize(0.05)
    hframe.GetXaxis().SetLabelSize(0.04)
    hframe.GetYaxis().SetLabelSize(0.04)

    hist.Draw('hist same')

    latex = TLatex()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatexNDC(0.55, 0.75, '^{90}Sr acquisition')
    #latex.DrawLatexNDC(0.55, 0.7, f'Threshold: Configuration 2')
    latex.DrawLatexNDC(0.55, 0.7, f'Threshold: Configuration 3')
    
    #canvas.SaveAs('../Data/output/source_clusterSize.pdf')
    canvas.SaveAs('../Data/output/source_clusterSize_lowTH.pdf')
    del canvas

def main():

    #inFile = TFile.Open('../Data/output/source_123218_analysis_chip113.root')
    inFile = TFile.Open('../Data/output/source_130551_analysis_chip113.root')
    #inFileCluster = TFile.Open('../Data/output/source_123218.root')
    inFileCluster = TFile.Open('../Data/output/source_130551.root')

    drawXProjection(inFile)
    drawYProjection(inFile)
    drawAcquisition(inFile)
    drawHitsAcquisition(inFile)
    drawClusterSize(inFileCluster)

if __name__ == '__main__':

    main()