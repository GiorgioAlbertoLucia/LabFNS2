#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include "../../yaml/Yaml.hpp"
#include "lib/spectrum.hh"

void spectralAnalysis()
{
    const char * inFilePath = "ITS3/Data/run175174828_230428174901_seed_and_cluster.root";
    const char * inTreeName = "SeedAndCluster";

    const char * outFilePath = "ITS3/Data/output/spectralAnalysisOutput.root";
    const char * cfgPath = "ITS3/src/config/spectrumFitCfg.yml";
    const char * outLogPath = "ITS3/Data/output/spectralAnalysisOutput.log";

    auto inFile = TFile::Open(inFilePath);
    auto inTree = (TTree*)inFile->Get(inTreeName);
    auto outFile = new TFile(outFilePath, "recreate");

    std::ofstream createOutput(outLogPath);
    createOutput.close();

    Yaml::Node cfg;
    Yaml::Parse(cfg, cfgPath);
    
    // -------------------------------------
    // seed 
    Yaml::Node cfgSeed = cfg["seed"];

    // -- seed drawing

    auto spectrumSeedFull = new TH1D("spectrumSeedFull",  "^{55}Fe spectrum - seed pixel; Seed signal (e); Counts (a.u.)", 
        cfgSeed["spectrumSpec"][0].As<int>(), cfgSeed["spectrumSpec"][1].As<double>(), cfgSeed["spectrumSpec"][2].As<double>());
    
    inTree->Draw("Seed.electrons>>spectrumSeedFull");
    for (int iBin = 1; iBin < spectrumSeedFull->GetEntries(); iBin++)   spectrumSeedFull->SetBinError(iBin, sqrt(spectrumSeedFull->GetBinContent(iBin)));
    outFile->cd();
    spectrumSeedFull->Write();

    // -- seed analysis

    auto spectrumSeed = new TH1D("spectrumSeed",  "^{55}Fe spectrum - seed pixel; Seed signal (e); Counts (a.u.)", 
        cfgSeed["spectrumSpec"][0].As<int>(), cfgSeed["spectrumSpec"][1].As<double>(), cfgSeed["spectrumSpec"][2].As<double>());
    
    inTree->Draw("Seed.electrons>>spectrumSeed", cfgSeed["selection"].As<std::string>().c_str());
    for (int iBin = 1; iBin < spectrumSeed->GetEntries(); iBin++)   spectrumSeed->SetBinError(iBin, sqrt(spectrumSeed->GetBinContent(iBin)));
    
    auto sSeed = new Spectrum(spectrumSeed, outFile);
    sSeed->EstimateBackground();
    //sSeed->DrawSpectrumAndBackground();
    
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::ofstream outputFile(outLogPath, std::ios::app);
    std::cout.rdbuf(outputFile.rdbuf());
    std::cout << "Seed:" << std::endl;
    std::cout.rdbuf(originalCoutBuffer);

    sSeed->SimpleFitSpectrum(outLogPath, cfgSeed);

    delete sSeed;

    // -------------------------------------
    // cluster
    Yaml::Node cfgCluster = cfg["cluster"];

    auto spectrumCluster = new TH1D("spectrumCluster",  "^{55}Fe spectrum - cluster pixel; Cluster signal (e); Counts (a.u.)", 
        cfgCluster["spectrumSpec"][0].As<int>(), cfgCluster["spectrumSpec"][1].As<double>(), cfgCluster["spectrumSpec"][2].As<double>());
    
    inTree->Draw("Cluster.electrons>>spectrumCluster", cfgCluster["selection"].As<std::string>().c_str());
    for (int iBin = 1; iBin < spectrumCluster->GetEntries(); iBin++)   spectrumCluster->SetBinError(iBin, sqrt(spectrumCluster->GetBinContent(iBin)));
    
    auto sCluster = new Spectrum(spectrumCluster, outFile);
    sCluster->EstimateBackground();
    //s->DrawSpectrumAndBackground();
    
    std::cout.rdbuf(outputFile.rdbuf());
    std::cout << "Cluster:" << std::endl;
    std::cout.rdbuf(originalCoutBuffer);

    sCluster->SimpleFitSpectrum(outLogPath, cfgCluster);

    delete sCluster;
    
    outFile->Close();
    inFile->Close();
}