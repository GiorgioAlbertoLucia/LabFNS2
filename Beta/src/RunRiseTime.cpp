#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TBranch.h>

#include <vector>
#include <algorithm>

#include "rise_time.hh"

TH1D* BuildRiseTime(const char* filename, const char* treename, const char* wfm_branchname, const char* tf_branchname) 
{
    
    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get(treename);

    TBranch *b_waveform, *b_timeFrame;
    std::vector<double> *waveform = 0;
    std::vector<double> *timeFrame = 0;
    tree->SetBranchAddress(wfm_branchname, &waveform, &b_waveform);
    tree->SetBranchAddress(tf_branchname, &timeFrame, &b_timeFrame);

    auto riseTime = new TH1D("riseTime", "Rise Time; Rise Time (ns); Counts", 20, 0, 1);
    riseTime->SetDirectory(0);

    // Loop over the entries in the TTree
    for (Long64_t ientry = 0; ientry < tree->GetEntries(); ++ientry) 
    {
        tree->GetEntry(ientry);
        
        // convert timeframes in ns
        const double factor = 1e9;
        std::transform(timeFrame->begin(), timeFrame->end(), timeFrame->begin(), [factor](double element) {return element*factor;});
        
        RiseTime rt(0.1, 10, *waveform, *timeFrame);
        double rise_time = rt.findRiseTime();
        
        riseTime->Fill(rise_time);
        if(ientry == 105)   rt.drawWaveform("Beta/data/output/RiseTime.pdf", ientry);
    }

    file->Close();

    return riseTime;
}

void RunRiseTime() 
{
    auto riseTime = BuildRiseTime("Beta/data/input/BetaMeas_Lab2.root", "wfm", "w3", "t3");
    riseTime->Draw();
}