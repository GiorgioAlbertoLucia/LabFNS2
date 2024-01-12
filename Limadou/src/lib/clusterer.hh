/*


*/

#ifndef _CLUSTERER_HH_
#define _CLUSTERER_HH_

#include <TFile.h>

#include "hit.hh"
#include "fileManager.hh"

#include <vector>

class Clusterer
{

    public:
        Clusterer(const char * inFilePath, const char * outFilePath);
        ~Clusterer();

        void FindHits();
        void VisualizeHits();
        void VisualizeClusterSize();

    private:

        FileManager * fInputData;
        TFile * fOutFile;

        std::vector<Hit> fHits;
        int fNEvents;
        
        int fNChips;
        int * fChipIDs;


};

#endif