#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TString.h>
#include <TLatex.h>
#include <TF1.h>
#include <TLegend.h>

#include <Riostream.h>
#include <vector>
#include <algorithm>

#include "rise_time.hh"
#include "rise_time_Analysis.hh"

void doubleCanvas(const char * outputPath)
{
    TFile * preprocessedFile = TFile::Open(outputPath);

    TH1D * RThist2 = (TH1D*)gDirectory->Get("RThist2");
    TH1D * RThist3 = (TH1D*)gDirectory->Get("RThist3");

    TCanvas * c1 = new TCanvas("c1", "Rise Time", 2100, 1800);
    c1->Divide(2,1);

    // ------------------------------------------------------------

    TPad * pad2 = (TPad*)c1->cd(1);

    TF1 * fit2 = new TF1("gaus", "gaus", 440, 545);
    fit2->SetLineColor(kRed);
    fit2->SetParameter(1, RThist2->GetMaximumBin());
    fit2->SetParameter(2, RThist2->GetRMS());
    RThist2->Fit(fit2, "rm+");
    std::cout << "chi2 / NDF  = " << fit2->GetChisquare() << " / " << fit2->GetNDF() << "\n";
    RThist2->SetFillStyle(3356);
    
    TH1F * frame2 = pad2->DrawFrame(350, 0, 900, 440, "Rise Time - Channel 2; Rise Time (ps); Counts");
    //frame2->GetYaxis()->SetTitleOffset(1.5);
    frame2->GetYaxis()->SetTitleSize(0.03);
    frame2->GetXaxis()->SetTitleSize(0.03);
    frame2->GetXaxis()->SetLabelSize(0.03);
    frame2->GetYaxis()->SetLabelSize(0.03);

    RThist2->Draw("same");

    TLatex * latex2 = new TLatex();
    latex2->SetTextSize(0.035);
    latex2->DrawLatexNDC(0.45, 0.75, "^{90}Sr acquisition");
    latex2->DrawLatexNDC(0.45, 0.7, Form("Rise Time = %.f #pm %.f ps", fit2->GetParameter(1), fit2->GetParameter(2)));
    latex2->DrawLatexNDC(0.45, 0.5, "Amplitude threshold = 80 mV");

    TLegend * leg2 = new TLegend(0.45, 0.55, 0.8, 0.65);
    leg2->AddEntry(RThist2, "Channel 2", "f");
    leg2->AddEntry(fit2, "Gaussian fit", "l");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->SetFillStyle(0);
    leg2->Draw();

    // ---------------------------------------

    TPad * pad3 = (TPad*)c1->cd(2);

    TF1 * fit3 = new TF1("gaus", "gaus", 440, 550);
    fit3->SetLineColor(kRed);
    fit3->SetParameter(1, RThist3->GetMaximumBin());
    //fit3->SetParLimits(1, 350, 390);
    fit3->SetParameter(2, 2*RThist3->GetRMS());
    RThist3->Fit(fit3, "rm+");
    std::cout << "chi2 / NDF  = " << fit3->GetChisquare() << " / " << fit3->GetNDF() << "\n";
    RThist3->SetFillStyle(3356);

    TH1F * frame3 = pad3->DrawFrame(300, 0, 1000, 4400, "Rise Time - Channel 3; Rise Time (ps); Counts");
    frame3->GetYaxis()->SetTitleOffset(1.7);
    frame3->GetYaxis()->SetTitleSize(0.03);
    frame3->GetXaxis()->SetTitleSize(0.03);
    frame3->GetXaxis()->SetLabelSize(0.03);
    frame3->GetYaxis()->SetLabelSize(0.03);

    RThist3->Draw("same");

    TLatex * latex3 = new TLatex();
    latex3->SetTextSize(0.035);
    latex3->DrawLatexNDC(0.45, 0.75, "^{90}Sr acquisition");
    latex3->DrawLatexNDC(0.45, 0.7, Form("Rise Time = %.f #pm %.f ps", fit3->GetParameter(1), fit3->GetParameter(2)));
    latex3->DrawLatexNDC(0.45, 0.5, "Amplitude threshold = 80 mV");

    TLegend * leg3 = new TLegend(0.45, 0.55, 0.8, 0.65);
    leg3->AddEntry(RThist3, "Channel 3", "f");
    leg3->AddEntry(fit3, "Gaussian fit", "l");
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.04);
    leg3->SetFillStyle(0);
    leg3->Draw();

    c1->SaveAs("Beta/data/output/RiseTime80mV.pdf");

    //delete c1;
    //delete fit2;
    //delete fit3;
    //delete RThist2;
    //delete RThist3;
    //delete latex2;
    //delete latex3;
    //delete leg2;
    //delete leg3;
    //delete dir;
    //delete preprocessedFile;
}

/**
 * @brief 
 * 
 */
void ComputeRiseTime(const bool drawOnly = false) 
{

    const char * wfmPath = "Beta/data/input/BetaMeas_Lab2.root";
    const char * wfmTreeName = "wfm";
    const char * preprocessedPath = "Beta/data/output/BetaOutput.root";
    const char * preprocessedTreeName = "BetaTree";
    const char * outputPath = "Beta/data/output/RiseTime.root";

    if(!drawOnly)
    {
        RiseTimeAnalysis RTanalysis(wfmPath, wfmTreeName, preprocessedPath, preprocessedTreeName, outputPath);

        RTanalysis.buildRiseTime("w2", "t2", "baseline2", 0.2, 10, 0.08, "Beta/data/output/RiseTimeCh2.pdf");
        RTanalysis.saveRiseTime("rt2");
        RTanalysis.analyseRiseTime(2);

        RTanalysis.buildRiseTime("w3", "t3", "baseline3", 0.2, 10, 0.08, "Beta/data/output/RiseTimeCh3.pdf");
        RTanalysis.saveRiseTime("rt3");
        RTanalysis.analyseRiseTime(3);
    }

    doubleCanvas(outputPath);
    
}