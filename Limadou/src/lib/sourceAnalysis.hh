/*

*/

#ifndef SOURCEANALYSIS_HH
#define SOURCEANALYSIS_HH

#include <string>

#include <TH2D.h>
#include <TFile.h>

#include "../../../yaml/Yaml.hpp"

class SourceAnalysis 
{
    public:
        SourceAnalysis(const int chipID, const char * inFilePath, const char * outFilePath = "default");
        ~SourceAnalysis();

        void RemoveNoisyPixels(const double threshold = 99999.);
        void SubtractBackground(const char * bkgFilePath);

        void FitHits(const char * cfgFitFile, const char * outLogPath = "default");
        void Save();

    protected: 
        void FitFromConfig(Yaml::Node & cfgFit, TH1D * hist);

    private:    

        int fChipID;
        TH2D * fHits;
        TH2D * fPixelMap;

        std::string fInFilePath;
        TFile * fOutFile;

};


#endif