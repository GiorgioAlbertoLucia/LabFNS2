/*

*/

#ifndef SOURCEANALYSIS_HH
#define SOURCEANALYSIS_HH

#include <string>

#include <TH2I.h>
#include <TFile.h>

class SourceAnalysis 
{
    public:
        SourceAnalysis(const int chipID, const char * inFilePath, const char * outFilePath = "default");
        ~SourceAnalysis();

        void FitHits(const char * cfgFitFile, const char * outLogPath = "default");
        void SubtractBackground(const char * bkgFilePath);

    private:    

        int fChipID;
        TH2I * fHits;

        std::string fInFilePath;
        TFile * fOutFile;

};


#endif