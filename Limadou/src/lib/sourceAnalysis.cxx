/*

*/

#include "sourceAnalysis.hh"

#include <Riostream.h>
#include <TString.h>
#include <TH1D.h>

#include "../../../yaml/Yaml.h"

/*  PUBLIC  */

SourceAnalysis::SourceAnalysis(const int chipID, const char * inFilePath, const char * outFilePath):
    fChipID(chipID),
    fInFilePath(inFilePath)
{
    if (outFilePath == "default")
    {
        TString inFileName(inFilePath);
        inFileName.ReplaceAll(".root", "");
        outFilePath = Form("%s_analysis_chip%d.root", inFileName.Data(), fChipID);
    }
    fOutFile = new TFile(outFilePath, "recreate");

    auto inFile = TFile::Open(inFilePath);
    fHits = (TH2I*)inFile->Get(Form("cihip%d", fChipID));
}

SourceAnalysis::~SourceAnalysis()
{
    fOutFile->Close();
    delete fHits;
}   

void SourceAnalysis::FitHits(const char * cfgFitFile, const char * outLogPath)
{
    if (outLogPath == "default")
    {
        TString logFileName(fInFilePath);
        logFileName.ReplaceAll(".root", "");
        outLogPath = Form("%s_analysis_chip%d.log", logFileName.Data(), fChipID);
    }

    YAML::Node cfgFits = YAML::LoadFile(cfgFitFile);

    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::ofstream outputFile(outLogPath, std::ios_base::app);
    std::cout.rdbuf(outputFile.rdbuf());

    std::cout << "  fits:" << std::endl;

    auto cfgFit = cfgFits["fits"][0];
    TH1D * projX = fHits->ProjectionX();
    SourceAnalysis::FitFromConfig(cfgFit, *projX);

    auto cfgFit = cfgFits["fits"][1];
    TH1D * projY = fHits->ProjectionY();
    SourceAnalysis::FitFromConfig(cfgFit, *projY);

    fOutFile->cd();
    projX->Write();
    projY->Write();

    delete projX;
    delete projY;

    std::cout.rdbuf(originalCoutBuffer);
}

/*  PROTECTED   */
void SourceAnalysis::FitFromConfig(YAML::Node & cfgFit, TH1D & hist)
{
    TF1 fit(cfgFit["name"].As<std::string>().c_str(), ["func"].As<std::string>().c_str(), cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());
    for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    fit.SetParameter(iParam, cfgFit["param"][iParam].As<double>());
    hist->Fit(&fit, "RM+", "", cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());

    // write fit results to log file
    std::cout << "    - name: " << cfgFit["name"].As<std::string>.c_str() << std::endl;
    std::cout << "      params:" << std::endl;
    for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    std::cout << "        - " << fit.GetParameter(iParam) << std::endl;
    std::cout << "      parerrors:" << std::endl;
    for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    std::cout << "        - " << fit.GetParError(iParam) << std::endl;
    std::cout << "      chi2: " << fit.GetChisquare() << std::endl;
    std::cout << "      NDF: " << fit.GetNDF() << std::endl;
    std::cout << std::endl;

    fit.SetLineColor(cfgFit["color"].As<int>());
    fit.SetRange(cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());
}