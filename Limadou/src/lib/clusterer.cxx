/*

*/

#include <TFile.h>
#include <TH2I.h>

#include "clusterer.hh"

Clusterer::Clusterer(const char * inFilePath):
    fInputData(new FileManager(inFilePath)),
    fNEvents(0),
    fNChips(0),
    fChipIDs(NULL)
{
    fNEvents = fInputData->max("event_count");
    double * chipIDs = fInputData->unique("chip_id(decimal)", fNChips);
    
    fChipIDs = new int[fNChips];
    for (int i = 0; i < fNChips; i++)   fChipIDs[i] = (int) chipIDs[i];

    delete [] chipIDs;
}

Clusterer::~Clusterer() 
{   
    delete fInputData;
  
}

void Clusterer::VisualizeHits(const char * outFilePath)
{
    TFile outFile(outFilePath, "recreate");

    for (int ichip = 0; ichip < fNChips; ichip++)
    {
        TH2I * h2 = new TH2I(Form("chip%d", fChipIDs[ichip]), Form("Chip %d; x (position); y (position)", fChipIDs[ichip]), 1024, 0, 1024, 512, 0, 512);
        for (int irow = 0; irow < fInputData->getNRows(); irow++)
        {
            double * x = fInputData->getColumn("hit_x");
            double * y = fInputData->getColumn("hit_y");
            double * chip_id = fInputData->getColumn("chip_id(decimal)");

            if (chip_id[irow] == fChipIDs[ichip])   h2->Fill(x[irow], y[irow]);
        }
        h2->Write();
        delete h2;
    }

    TH2I * h2 = new TH2I("all_chips", "All Chips; x (position); y (position)", 1024, 0, 1024, 512, 0, 512);
    for (int irow = 0; irow < fInputData->getNRows(); irow++)
    {
        double * x = fInputData->getColumn("hit_x");
        double * y = fInputData->getColumn("hit_y");
        h2->Fill(x[irow], y[irow]);
    }

    h2->Write();
    delete h2;
}