#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

void Run(TString myopt="")
{
    TString opt;
    if (myopt.Contains("force"))
        opt="kfg-";          //k mantiene il .so, f forza la compilazione (come mettere il ++), g serve se vogliamo usare un debugger
    else
        opt="kg-";
    if (myopt.Contains("clean"))
        gSystem->Exec("./clean.sh");
    //gSystem->CompileMacro("Beta/src/analysis.C",opt.Data(),"","Build");

    gROOT->ProcessLine(".L Beta/src/analysis.C");
    gROOT->ProcessLine("analysis a; a.Loop();");
}