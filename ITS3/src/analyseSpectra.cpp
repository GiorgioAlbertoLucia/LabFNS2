<<<<<<< HEAD:ITS3/src/testSpectrum.cpp
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include "../yaml/Yaml.hpp"

#include "spectrum.hh"

void testSpectrum()
{
    const char * inFilePath = "ITS3/Data/run175174828_230428174901_preprocessed.root";
    const char * inTreeName = "PreprocessedData";

    const char * outFilePath = "ITS3/Data/output/spectrumOutput.root";
    const char * cfgPath = "ITS3/src/config/simpleFitCfg.yml";
    const char * outLogPath = "ITS3/Data/output/spectrumOutput.log";

    auto inFile = TFile::Open(inFilePath);
    auto inTree = (TTree*)inFile->Get(inTreeName);
    auto outFile = new TFile(outFilePath, "recreate");

    int pixels[] = {5, 6, 9, 10};
    for (int iPixel: pixels)
    {
        Yaml::Node cfg;
        Yaml::Parse(cfg, cfgPath);
        
        auto spectrum = new TH1D(Form("spectrumPx%d", iPixel), "", 160, 1, 81);
        inTree->Draw(Form("pixel%d.amplitude>>spectrumPx%d", iPixel, iPixel));
        for (int iBin = 1; iBin < spectrum->GetEntries(); iBin++)   spectrum->SetBinError(iBin, sqrt(spectrum->GetBinContent(iBin)));

        auto s = new Spectrum(spectrum, outFile);

        s->EstimateBackground();
        //s->DrawSpectrumAndBackground();

        s->SimpleFitSpectrum(outLogPath, cfg[Form("pixel%d", iPixel)]);

        delete s;
    }

    outFile->Close();
=======
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include "../yaml/Yaml.hpp"

#include "spectrum.hh"

void analyseSpectra()
{
    const char * inFilePath = "ITS3/Data/run175174828_230428174901_preprocessed.root";
    const char * inTreeName = "PreprocessedData";

    const char * outFilePath = "ITS3/Data/output/spectrumOutput.root";
    const char * cfgPath = "ITS3/src/config/simpleFitCfg.yml";
    const char * outLogPath = "ITS3/Data/output/spectrumOutput.log";

    auto inFile = TFile::Open(inFilePath);
    auto inTree = (TTree*)inFile->Get(inTreeName);
    auto outFile = new TFile(outFilePath, "recreate");

    std::ofstream createOutput(outLogPath);
    createOutput.close();

    int pixels[] = {5, 6, 9, 10};
    for (int iPixel: pixels)
    {
        Yaml::Node cfg;
        Yaml::Parse(cfg, cfgPath);
        
        auto spectrum = new TH1D(Form("spectrumPx%d", iPixel), "", cfg[Form("pixel%d", iPixel)]["spectrumSpec"][0].As<int>(), 
            cfg[Form("pixel%d", iPixel)]["spectrumSpec"][1].As<double>(), cfg[Form("pixel%d", iPixel)]["spectrumSpec"][2].As<double>());
        
        Yaml::Node cfgHist = cfg[Form("pixel%d", iPixel)]["spectrumSpec"];
        inTree->Draw(Form("pixel%d.amplitude>>spectrumPx%d", iPixel, iPixel));
        for (int iBin = 1; iBin < spectrum->GetEntries(); iBin++)   spectrum->SetBinError(iBin, sqrt(spectrum->GetBinContent(iBin)));

        auto s = new Spectrum(spectrum, outFile);

        s->EstimateBackground();
        //s->DrawSpectrumAndBackground();

        std::streambuf* originalCoutBuffer = std::cout.rdbuf();
        std::ofstream outputFile(outLogPath, std::ios::app);
        std::cout.rdbuf(outputFile.rdbuf());
        std::cout << "pixel" << iPixel << ":" << std::endl;
        std::cout.rdbuf(originalCoutBuffer);

        s->SimpleFitSpectrum(outLogPath, cfg[Form("pixel%d", iPixel)]);

        delete s;
    }

    outFile->Close();
    inFile->Close();
>>>>>>> 426078767d4e7cb224e63f70c50666d91a8e9623:ITS3/src/analyseSpectra.cpp
}