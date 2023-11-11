#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TString.h>
#include <TLatex.h>
#include <TF1.h>
#include <TLegend.h>

#include <vector>
#include <algorithm>

#include "rise_time.hh"
#include "rise_time_Analysis.hh"

void doubleCanvas(const char * preprocessedPath)
{
    TFile * preprocessedFile = TFile::Open(preprocessedPath);
    TDirectory * dir = preprocessedFile->GetDirectory("RiseTimeAnalysis");
    dir->cd();

    TH1D * RThist2 = (TH1D*)gDirectory->Get("RThist2");
    TH1D * RThist3 = (TH1D*)gDirectory->Get("RThist3");

    TCanvas * c1 = new TCanvas("c1", "Rise Time", 800, 600);
    c1->Divide(2,1);
    
    c1->cd(1);
    RThist2->Draw();
    TF1 * fit2 = RThist2->GetFunction("gaus");

    TLatex * latex2 = new TLatex();
    latex2->SetTextSize(0.04);
    latex2->SetTextColor(kBlue);
    latex2->DrawLatexNDC(0.6, 0.85, "^{90}Sr acquisition");
    latex2->DrawLatexNDC(0.6, 0.8, Form("Rise Time = %.2f #pm %.2f ps", fit2->GetParameter(1), fit2->GetParameter(2)));

    TLegend * leg2 = new TLegend(0.6, 0.6, 0.8, 0.7);
    leg2->AddEntry(RThist2, "Channel 2", "f");
    leg2->AddEntry(fit2, "Gaussian fit", "l");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->Draw();


    c1->cd(2);
    RThist3->Draw();
    TF1 * fit3 = RThist3->GetFunction("gaus");

    TLatex * latex3 = new TLatex();
    latex3->SetTextSize(0.04);
    latex3->SetTextColor(kBlue);
    latex3->DrawLatexNDC(0.6, 0.85, "^{90}Sr acquisition");
    latex3->DrawLatexNDC(0.6, 0.8, Form("Rise Time = %.2f #pm %.2f ps", fit3->GetParameter(1), fit3->GetParameter(2)));

    TLegend * leg3 = new TLegend(0.6, 0.6, 0.8, 0.7);
    leg3->AddEntry(RThist3, "Channel 3", "f");
    leg3->AddEntry(fit3, "Gaussian fit", "l");
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->Draw();

    c1->SaveAs("Beta/data/output/RiseTime.pdf");

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
void RunRiseTime() 
{
    const char * wfmPath = "Beta/data/input/BetaMeas_Lab2.root";
    const char * wfmTreeName = "wfm";
    const char * preprocessedPath = "Beta/data/output/BetaOutput.root";
    const char * preprocessedTreeName = "BetaTree";

    RiseTimeAnalysis RTanalysis(wfmPath, wfmTreeName, preprocessedPath, preprocessedTreeName);

    RTanalysis.buildRiseTime("w2", "t2", "baseline2", 0.2, 10, "Beta/data/output/RiseTimeCh2.pdf");
    RTanalysis.saveRiseTime("rt2");
    RTanalysis.analyseRiseTime(2);

    RTanalysis.buildRiseTime("w3", "t3", "baseline3", 0.2, 10, "Beta/data/output/RiseTimeCh3.pdf");
    RTanalysis.saveRiseTime("rt3");
    RTanalysis.analyseRiseTime(3);

    doubleCanvas(preprocessedPath);
    
}