#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TBranch.h>

#include <vector>
#include <algorithm>

#include "rise_time.hh"

/**
 * @brief Extract rise time from the waveform and return it as a vector
 * 
 * @param filename 
 * @param treename 
 * @param wfm_branchname 
 * @param tf_branchname 
 * @param baseline_filename 
 * @param baseline_treename 
 * @param b_branchname 
 * @return TH1D* 
 */
std::vector<double> BuildRiseTime(const char* filename, const char* treename, const char* wfm_branchname, const char* tf_branchname,
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

    std::vector<double> riseTime;
    riseTime.reserve(tree->GetEntries());

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
        
        riseTime.push_back(rise_time);
        if(ientry == 105)   rt.drawWaveform("Beta/data/output/RiseTime.pdf", ientry);
    }

    file->Close();
    return riseTime;
}

/**
 * @brief Add double branch in an existing TTree
 * 
 * @param variable
 * @param filename
 * @param tree_name
 * @param branchname
 */
void saveRiseTime(std::vector<double> const & variable, const char * filename, const char* tree_name, const char * branchname)
{
    TFile* file = new TFile(filename, "update");
    TTree* tree = (TTree*)file->Get(tree_name);

    double var = 0;
    TBranch *b_var = tree->Branch(branchname, &var, Form("%s/D", branchname));
    tree->SetBranchAddress(branchname, &var, &b_var);

    for(auto ientry=0; ientry<tree->GetEntries(); ++ientry)
    {
        var = variable.at(ientry);
        b_var->Fill();
    }

    tree->Print();
    tree->Write();

    delete tree;
    file->Close();
    delete file;
}

/**
 * @brief Analysis on rise time
 * 
 */
void analysisRiseTime(std::vector<double> const & riseTime)
{
    auto riseTimeHist = new TH1D("riseTime", "Rise Time; Rise Time (ps); Counts", 72, -200, 1000);
    for(auto& rt : riseTime)   riseTimeHist->Fill(rt*1000);

    // analysis
    riseTimeHist->Draw();
}

/**
 * @brief 
 * 
 */
void RunRiseTime() 
{
    const char * filename = "Beta/data/input/BetaMeas_Lab2.root";
    const char * treename = "wfm";
    const char * wfm_branchname = "w3";
    const char * tf_branchname = "t3";
    const char * baseline_filename = "Beta/data/output/BetaOutput.root";
    const char * baseline_treename = "BetaTree";
    const char * b_branchname = "baseline3";
    
    auto riseTime = BuildRiseTime(filename, treename, wfm_branchname, tf_branchname, baseline_filename, baseline_treename, b_branchname);
    analysisRiseTime(riseTime);
    saveRiseTime(riseTime, baseline_filename, baseline_treename, "rt3");
}