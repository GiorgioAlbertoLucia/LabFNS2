#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TString.h>
#include <TLatex.h>
#include <TF1.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TROOT.h>

/**
 * @brief 
 * 
 */
void RunRiseTime(TString myopt = "", const bool drawOnly = false) 
{

    gSystem->AddIncludePath((string("-I ")+gSystem->GetWorkingDirectory()+"Beta/src/build").c_str());
    
    TString opt;
    if (myopt.Contains("force"))    opt = "kfg-";
    else                            opt = "kg-";
    if (myopt.Contains("clean"))    gSystem->Exec("./Beta/src/clean.sh");

    gSystem->CompileMacro("Beta/src/rise_time.cxx",opt.Data(),"","Beta/src/build");
    gSystem->CompileMacro("Beta/src/rise_time_Analysis.cxx",opt.Data(),"","Beta/src/build");
    gSystem->CompileMacro("Beta/src/ComputeRiseTime.cpp",opt.Data(),"","Beta/src/build");

    if(!drawOnly)
    {
        gSystem->CompileMacro("Beta/src/analysis.C",opt.Data(),"","Beta/src/build");
        gROOT->ProcessLine("analysis a; a.Loop();");
    }

    gROOT->ProcessLine("ComputeRiseTime();");
}