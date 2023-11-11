#ifndef RISE_TIME_ANALYSIS_HH
#define RISE_TIME_ANALYSIS_HH

#include <string>
#include <vector>

#include <TDirectory.h>

#include "rise_time.hh"


class RiseTimeAnalysis
{
    public:
        RiseTimeAnalysis(const char * wfmPath, const char * wfmTreeName, const char * preprocessedPath, const char * preprocessedTreeName);
        ~RiseTimeAnalysis();

        void buildRiseTime(const char* w_branchname, const char* t_branchname, const char* b_branchname, 
                            const double threshold = 0.2, const int window = 10, const char* wfmDrawPath = "Beta/data/output/RiseTime.pdf"); 
        void saveRiseTime(const char * rt_branchname);
        void analyseRiseTime(const int channel);


    private:
        std::string fWfmPath;               // Path to the waveform file
        std::string fWfmTreeName;           // Name of the tree in the waveform file
        std::string fPreprocessedPath;      // Path to the preprocessed file
        std::string fPreprocessedTreeName;  // Name of the tree in the preprocessed file
        std::string fRiseTimeBranchName;    // Name of the rise time branch in the preprocessed file

        std::vector<double> fRiseTime;      // Vector to store the rise time

        std::string fDirectoryName;         // Name of the directory with results

};


#endif