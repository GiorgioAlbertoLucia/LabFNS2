#include "rise_time_Analysis.hh"
#include "rise_time.hh"

#include <Riostream.h>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>

ClassImp(RiseTimeAnalysis)

RiseTimeAnalysis::RiseTimeAnalysis(const char * wfmPath, const char * wfmTreeName, const char * preprocessedPath, const char * preprocessedTreeName, const char * outputPath):
fWfmPath(wfmPath),
fWfmTreeName(wfmTreeName),
fPreprocessedPath(preprocessedPath),
fPreprocessedTreeName(preprocessedTreeName),
fOutputPath(outputPath),
fRiseTimeBranchName("")
{
    fRiseTime = std::vector<double>();
    TFile outputFile(fOutputPath.c_str(), "recreate");
    outputFile.Close();
}

RiseTimeAnalysis::~RiseTimeAnalysis()
{

}

/*
    PUBLIC METHODS
*/

/**
 * @brief Build rise time
 * 
 * @param w_branchname  Name of the waveform branch in the waveform file
 * @param t_branchname  Name of the time frame branch in the waveform file
 * @param b_branchname  Name of the baseline branch in the preprocessed file
 * @param threshold     Threshold for rise time calculation in percent of max amplitude
 * @param window        Window for rise time calculation in samples (in timeframes)
 * @param wfmDrawPath   Path to the waveform plot
 */
void RiseTimeAnalysis::buildRiseTime(const char* w_branchname, const char* t_branchname, const char* b_branchname, const double threshold, const int window, const double minAmplitude, const char* wfmDrawPath)
{
    // input waveform
    TFile* wfmFile = TFile::Open(fWfmPath.c_str());
    TTree* wfmTree = (TTree*)wfmFile->Get(fWfmTreeName.c_str());

    TBranch *b_waveform, *b_timeFrame;
    std::vector<double> *waveform = 0;
    std::vector<double> *timeFrame = 0;
    wfmTree->SetBranchAddress(w_branchname, &waveform, &b_waveform);
    wfmTree->SetBranchAddress(t_branchname, &timeFrame, &b_timeFrame);

    // input baseline
    TFile* preprocessedFile = TFile::Open(fPreprocessedPath.c_str());
    TTree* preprocessedTree = (TTree*)preprocessedFile->Get(fPreprocessedTreeName.c_str());

    TBranch *b_baseline;
    double baseline = 0;
    preprocessedTree->SetBranchAddress(b_branchname, &baseline, &b_baseline);

    fRiseTime.clear();
    fRiseTime.reserve(wfmTree->GetEntries());

    // Loop over the entries in the TTree
    for (Long64_t ientry = 0; ientry < wfmTree->GetEntries(); ++ientry) 
    {
        if(ientry%1000==0 || ientry==wfmTree->GetEntries()-1)  std::cout << "Event " << ientry << std::endl;
        wfmTree->GetEntry(ientry);
        preprocessedTree->GetEntry(ientry);
        
        // convert timeframes in ns
        const double factor = 1e9;
        std::transform(timeFrame->begin(), timeFrame->end(), timeFrame->begin(), [factor](double element) {return element*factor;});

        // subtract baseline
        std::transform(waveform->begin(), waveform->end(), waveform->begin(), [baseline](double element) {return element-baseline;});
        
        RiseTime rt(threshold, window, *waveform, *timeFrame);
        double rise_time = rt.findRiseTime(minAmplitude);
        
        fRiseTime.push_back(rise_time);
        if(ientry == 100)   
        {
            rt.drawWaveform(wfmDrawPath, ientry, minAmplitude);
            TString wfmDrawPathFit = wfmDrawPath;
            wfmDrawPathFit.ReplaceAll(".pdf", "_fit.pdf"); 
            rt.drawLinearFit(wfmDrawPathFit.Data(), ientry, minAmplitude);
            TString wfmDrawPathDerivative = wfmDrawPath;
            wfmDrawPathDerivative.ReplaceAll(".pdf", "_derivative.pdf");
            rt.drawDerivative(wfmDrawPathDerivative.Data(), ientry, minAmplitude);
        }
    }

    preprocessedFile->Close();
    wfmFile->Close();

    //delete b_waveform;
    //delete b_timeFrame;
    //delete waveform;
    //delete timeFrame;
    //delete wfmTree;
    //delete wfmFile;
}

/**
 * @brief Add double branch in an existing TTree
 * 
 * @param branchname
 */
void RiseTimeAnalysis::saveRiseTime(const char * branchname)
{
    TFile preprocessedFile(fPreprocessedPath.c_str(), "update");
    TTree* preprocessedTree = (TTree*)preprocessedFile.Get(fPreprocessedTreeName.c_str());

    TBranch* existingBranch = preprocessedTree->GetBranch(branchname);
    if (existingBranch) {
        // If the branch exists, delete it
        preprocessedTree->GetListOfBranches()->Remove(existingBranch);
        delete existingBranch;
    }

    double var = 0;
    TBranch *b_var = preprocessedTree->Branch(branchname, &var, Form("%s/D", branchname));
    preprocessedTree->SetBranchAddress(branchname, &var, &b_var);

    for(auto ientry=0; ientry<preprocessedTree->GetEntries(); ++ientry)
    {
        var = fRiseTime.at(ientry);
        b_var->Fill();
    }

    preprocessedTree->Print();
    preprocessedTree->Write("", TObject::kOverwrite);

    //delete b_var;
    //delete preprocessedTree;
    preprocessedFile.Close();
}


/**
 * @brief Analysis on rise time
 * 
 */
void RiseTimeAnalysis::analyseRiseTime(const int channel)
{
    auto riseTimeHist = new TH1D(Form("RThist%d", channel), Form("Rise Time - Channel %d; Rise Time (ps); Counts", channel), 320, -2000, 2000);
    for(unsigned long ientry=0; ientry<fRiseTime.size(); ++ientry)  riseTimeHist->Fill(fRiseTime.at(ientry)*1000);

    riseTimeHist->SetStats(0);
    riseTimeHist->SetFillStyle(2);
    riseTimeHist->SetFillColor(kBlue-10);
    riseTimeHist->SetLineColor(kBlue);

    // save results
    TFile* outputFile = TFile::Open(fOutputPath.c_str(), "update");
    riseTimeHist->Write();
    outputFile->Close();

    //delete gaus;
    //delete riseTimeHist;
    //delete dir;
    //delete outputFile;
}