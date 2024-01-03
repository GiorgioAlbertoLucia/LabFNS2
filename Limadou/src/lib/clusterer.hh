/*


*/

#ifndef _CLUSTERER_HH_
#define _CLUSTERER_HH_

#include "fileManager.hh"

class Clusterer
{

    public:
        Clusterer(const char * inFilePath);
        ~Clusterer();

        void VisualizeHits(const char * outFilePath);

    private:

        FileManager * fInputData;
        int fNEvents;
        
        int fNChips;
        int * fChipIDs;


};

#endif