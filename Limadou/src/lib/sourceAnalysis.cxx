/*

*/

#include "sourceAnalysis.hh"

#include <Riostream.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>

#include "../../../yaml/Yaml.hpp"

/*  PUBLIC  */

SourceAnalysis::SourceAnalysis(const int chipID, const char * inFilePath, const char * outFilePath):
    fHits(NULL),
    fChipID(chipID),
    fInFilePath(inFilePath)
{
    if (strcmp(outFilePath, "default") == 0)
    {
        TString inFileName(inFilePath);
        inFileName.ReplaceAll(".root", "");
        outFilePath = Form("%s_analysis_chip%d.root", inFileName.Data(), fChipID);
    }
    fOutFile = new TFile(outFilePath, "recreate");

    auto inFile = TFile::Open(inFilePath);
    fHits = (TH2D*)inFile->Get(Form("chip%d", fChipID));
}

SourceAnalysis::~SourceAnalysis()
{
    fOutFile->Close();
    delete fHits;
}   

void SourceAnalysis::SubtractBackground(const char * bkgFilePath)
{
    auto bkgFile = TFile::Open(bkgFilePath);
    auto bkgHist = (TH2D*)bkgFile->Get(Form("chip%d", fChipID));

    // normalize the background histogram to the same number of entries as the hits histogram
    bkgHist->Scale(fHits->GetEntries() / bkgHist->GetEntries());
    
    fHits->Add(bkgHist, -1);
    for (int xbin = 1; xbin <= fHits->GetNbinsX(); xbin++)
    {
        for (int ybin = 1; ybin <= fHits->GetNbinsY(); ybin++)
        {
            if (fHits->GetBinContent(xbin, ybin) < 0)    fHits->SetBinContent(xbin, ybin, 0);
        }
    }

    delete bkgHist;
}

void SourceAnalysis::FitHits(const char * cfgFitFile, const char * outLogPath)
{
    if (strcmp(outLogPath, "default") == 0)
    {
        TString logFileName(fInFilePath);
        logFileName.ReplaceAll(".root", "");
        outLogPath = Form("%s_analysis_chip%d.log", logFileName.Data(), fChipID);
    }

    Yaml::Node cfgFits;
    Yaml::Parse(cfgFits, cfgFitFile);

    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::ofstream outputFile(outLogPath, std::ios_base::out);
    std::cout.rdbuf(outputFile.rdbuf());

    std::cout << "fits:" << std::endl;

    auto cfgFitX = cfgFits["fits"][0];
    TH1D * projX = fHits->ProjectionX(Form("chip%d_projX", fChipID), cfgFitX["projRange"][0].As<int>(), cfgFitX["projRange"][1].As<int>());
    projX->Rebin(cfgFitX["rebin"].As<int>());
    for (int ibin = 0; ibin < projX->GetNbinsX(); ibin++)   projX->SetBinError(ibin, sqrt(projX->GetBinContent(ibin)));
    SourceAnalysis::FitFromConfig(cfgFitX, projX);

    auto cfgFitY = cfgFits["fits"][1];
    TH1D * projY = fHits->ProjectionY(Form("chip%d_projY", fChipID), cfgFitY["projRange"][0].As<int>(), cfgFitY["projRange"][1].As<int>());
    projY->Rebin(cfgFitY["rebin"].As<int>());
    for (int ibin = 0; ibin < projY->GetNbinsX(); ibin++)   projY->SetBinError(ibin, sqrt(projY->GetBinContent(ibin)));
    SourceAnalysis::FitFromConfig(cfgFitY, projY);

    fOutFile->cd();
    projX->Write();
    projY->Write();

    delete projX;
    delete projY;

    std::cout.rdbuf(originalCoutBuffer);
}

void SourceAnalysis::Save()
{
    fOutFile->cd();
    fHits->Write();
}

/*  PROTECTED   */
void SourceAnalysis::FitFromConfig(Yaml::Node & cfgFit, TH1D * hist)
{
    TF1 fit(cfgFit["name"].As<std::string>().c_str(), cfgFit["func"].As<std::string>().c_str(), cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());
    for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    fit.SetParameter(iParam, cfgFit["param"][iParam].As<double>());
    hist->Fit(&fit, "LRM+", "", cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());

    // write fit results to log file
    std::cout << "  - name: " << cfgFit["name"].As<std::string>().c_str() << std::endl;
    std::cout << "    params:" << std::endl;
    for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    std::cout << "      - " << fit.GetParameter(iParam) << std::endl;
    std::cout << "    parerrors:" << std::endl;
    for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    std::cout << "      - " << fit.GetParError(iParam) << std::endl;
    std::cout << "    chi2: " << fit.GetChisquare() << std::endl;
    std::cout << "    NDF: " << fit.GetNDF() << std::endl;
    std::cout << std::endl;

    fit.SetLineColor(cfgFit["color"].As<int>());
    fit.SetRange(cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());
}