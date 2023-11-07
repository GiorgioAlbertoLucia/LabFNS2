#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TBranch.h>

#include <vector>
#include <algorithm>

#include "rise_time.hh"

TH1D* BuildRiseTime(const char* filename, const char* treename, const char* wfm_branchname, const char* tf_branchname,
                    const char* baseline_filename, const char* baseline_treename, const char* b_branchname) 
{
    // input waveform
    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get(treename);

    TBranch *b_waveform, *b_timeFrame;
    std::vector<double> *waveform = 0;
    std::vector<double> *timeFrame = 0;
    tree->SetBranchAddress(wfm_branchname, &waveform, &b_waveform);
    tree->SetBranchAddress(tf_branchname, &timeFrame, &b_timeFrame);

    // input baseline
    TFile* file_baseline = TFile::Open(baseline_filename);
    TTree* tree_baseline = (TTree*)file_baseline->Get(baseline_treename);

    TBranch *b_baseline;
    double baseline = 0;
    tree_baseline->SetBranchAddress(b_branchname, &baseline, &b_baseline);

    auto riseTime = new TH1D("riseTime", "Rise Time; Rise Time (ns); Counts", 20, 0, 1);
    riseTime->SetDirectory(0);

    // Loop over the entries in the TTree
    for (Long64_t ientry = 0; ientry < tree->GetEntries(); ++ientry) 
    {
        if(ientry%1000==0 || ientry==tree->GetEntries()-1)  std::cout << "Event " << ientry << std::endl;
        tree->GetEntry(ientry);
        tree_baseline->GetEntry(ientry);
        
        // convert timeframes in ns
        const double factor = 1e9;
        std::transform(timeFrame->begin(), timeFrame->end(), timeFrame->begin(), [factor](double element) {return element*factor;});

        // subtract baseline
        std::transform(waveform->begin(), waveform->end(), waveform->begin(), [baseline](double element) {return element-baseline;});
        
        RiseTime rt(0.2, 10, *waveform, *timeFrame);
        double rise_time = rt.findRiseTime();
        
        riseTime->Fill(rise_time);
        if(ientry == 105)   rt.drawWaveform("Beta/data/output/RiseTime.pdf", ientry);
    }

    file->Close();

    return riseTime;
}

void RunRiseTime() 
{
    const char * filename = "Beta/data/input/BetaMeas_Lab2.root";
    const char * treename = "wfm";
    const char * wfm_branchname = "w3";
    const char * tf_branchname = "t3";
    const char * baseline_filename = "Beta/data/input/BetaMeas_Lab2.root";
    const char * baseline_treename = "wfm";
    const char * b_branchname = "baseline3";
    
    auto riseTime = BuildRiseTime(filename, treename, wfm_branchname, tf_branchname, baseline_filename, baseline_treename, b_branchname);
    riseTime->Draw();
}