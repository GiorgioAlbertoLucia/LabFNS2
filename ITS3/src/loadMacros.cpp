void loadMacros(TString myopt="fast")
{
    gSystem->AddIncludePath((string("-I ")+gSystem->GetWorkingDirectory()+"ITS3/src/build").c_str());
    TString opt;
    if(myopt.Contains("force"))   opt = "kfg";
    else                          opt = "kg";
  
    // YAML
    gSystem->CompileMacro("ITS3/yaml/Yaml.cpp",opt.Data(), "", "ITS3/src/build");

    // PREPROCESSING
    gSystem->CompileMacro("ITS3/src/preprocessor.cxx",opt.Data(),"","ITS3/src/build");

    // SPECTRUM
    gSystem->CompileMacro("ITS3/src/spectrum.cxx",opt.Data(),"","ITS3/src/build");

}