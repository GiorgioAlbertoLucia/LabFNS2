#include <TString.h>
#include <TSystem.h>

#include <string>

void loadMacros(TString myopt="fast")
{
    gSystem->AddIncludePath((std::string("-I ")+gSystem->GetWorkingDirectory()+"ITS3/src/build").c_str());
    TString opt;
    if(myopt.Contains("force"))   opt = "kfg";
    else                          opt = "kg";
  
    // YAML
    //gSystem->CompileMacro("yaml/Yaml.cpp",opt.Data(), "", "ITS3/src/build");

    // FILEMANAGER
    gSystem->CompileMacro("Limadou/src/lib/fileManager.cxx",opt.Data(), "", "Limadou/src/build");

    // CLUSTERER
    gSystem->CompileMacro("Limadou/src/lib/clusterer.cxx",opt.Data(), "", "Limadou/src/build");
}