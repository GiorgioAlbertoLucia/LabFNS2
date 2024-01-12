/*

*/

#include <Riostream.h>
#include <unordered_map>
#include <algorithm>

#include <TFile.h>
#include <TH2D.h>

#include "clusterer.hh"

Clusterer::Clusterer(const char * inFilePath, const char * outFilePath):
    fInputData(new FileManager(inFilePath)),
    fOutFile(NULL),
    fNEvents(0),
    fNChips(0),
    fChipIDs(NULL)
{
    fNEvents = fInputData->max("event_count");
    double * chipIDs = fInputData->unique("chip_id(decimal)", fNChips);
    
    fChipIDs = new int[fNChips];
    for (int i = 0; i < fNChips; i++)   fChipIDs[i] = (int) chipIDs[i];

    fOutFile = new TFile(outFilePath, "recreate");

    delete [] chipIDs;
}

Clusterer::~Clusterer() 
{   
    fOutFile->Close();
    delete fInputData;
}

void Clusterer::FindHits()
{
    std::unordered_map<int, std::vector<int>> event_indices;

    double *event = fInputData->getColumn("event_count");
    double *chip_id = fInputData->getColumn("chip_id(decimal)");
    double *x = fInputData->getColumn("hit_x");
    double *y = fInputData->getColumn("hit_y");

    for (int irow = 0; irow < fInputData->getNRows(); irow++)
    {
        event_indices[event[irow]].push_back(irow);
    }

    /*
    for (const auto &[ievent, indices] : event_indices)
    {
        std::cout << "Event " << ievent << std::endl;
        for (int i = 0; i < indices.size(); i++)
        {
            std::cout << indices[i] << " ";
        }
    }


    */

    for (auto it = event_indices.begin(); it != event_indices.end(); ++it)
    //for (auto it = event_indices.begin(); it != std::next(event_indices.begin(), 2); ++it)
    {
        int ievent = it->first;
        std::vector<int>& indices = it->second;

        std::cout << "Event " << ievent << std::endl;

        std::cout << "Indices size: " << indices.size() << std::endl;

        while (!indices.empty())
        {

            std::vector<int> next_indices;
            double mean_ix = x[indices[0]], mean_iy = y[indices[0]];
            int cluster_size = 1;

            for (int i = 1; i < indices.size(); i++)
            {
                int index = indices[i];

                if (chip_id[index] != chip_id[indices[0]]) continue;
                if (std::abs(x[index] - x[indices[0]]) != 1 && std::abs(y[index] - y[indices[0]]) != 1) continue;

                bool isValid = true;
                for (int j : next_indices)
                {
                    if (std::abs(x[index] - x[j]) != 1 && std::abs(y[index] - y[j]) != 1)
                    {
                        isValid = false;
                        break;
                    }
                }

                if (isValid)
                {
                    next_indices.push_back(index);
                    mean_ix += x[index];
                    mean_iy += y[index];
                    cluster_size++;
                }
            }

            mean_ix /= cluster_size;
            mean_iy /= cluster_size;

            Hit hit;
            hit.chipID = static_cast<int>(chip_id[indices[0]]);
            hit.event = static_cast<int>(ievent);
            hit.x = mean_ix;
            hit.y = mean_iy;
            hit.clusterSize = cluster_size;
            fHits.push_back(hit);

            // Remove processed indices from the vector
            indices.erase(std::remove_if(indices.begin(), indices.end(),
                [&next_indices](int i) { return std::find(next_indices.begin(), next_indices.end(), i) != next_indices.end(); }),
                indices.end());
            indices.erase(indices.begin());
        }
    }
}

/*

void Clusterer::FindHits()
{
    double * event = fInputData->getColumn("event_count");
    double * chip_id = fInputData->getColumn("chip_id(decimal)");
    double * x = fInputData->getColumn("hit_x");
    double * y = fInputData->getColumn("hit_y");

    
    for (int ievent = 0; ievent < fNEvents; ievent++)
    {
        std::cout << "Event " << ievent << std::endl;
        std::vector<double> ix, iy, ichip_id;
        for (int irow = 0; irow < fInputData->getNRows(); irow++)
        {
            if (event[irow] == ievent)
            {
                ix.push_back(x[irow]);
                iy.push_back(y[irow]);
                ichip_id.push_back(chip_id[irow]);
            }
        }

        int it = 0; 
        while (!ix.empty())
        {
            std::cout << "Iteration " << it++ << std::endl;

            std::vector<int> erase_position;
            std::vector<double> temp_ix, temp_iy;
            double mean_ix = ix[0], mean_iy = iy[0];
            int cluster_size = 1;

            for (int i = 1; i < ix.size(); i++)
            {

                if (ichip_id[i] != ichip_id[0])                                                 continue;
                if (std::abs(ix[i] - ix[0]) != 1 && std::abs(iy[i] - iy[0]) != 1)               continue;
                for (int j = 0; j < temp_ix.size(); j++)
                {
                    if (std::abs(ix[i] - temp_ix[j]) != 1 && std::abs(iy[i] - temp_iy[j]) != 1) continue;
                }
                
                erase_position.push_back(i);
                temp_ix.push_back(ix[i]);
                temp_iy.push_back(iy[i]);
                mean_ix += ix[i];
                mean_iy += iy[i];
                cluster_size++;
                
            }

            mean_ix /= cluster_size;
            mean_iy /= cluster_size;
            Hit hit;
            hit.chipID = (int)ichip_id[0];
            hit.event = (int)event[0];
            hit.x = mean_ix;
            hit.y = mean_iy;
            hit.clusterSize = cluster_size;
            fHits.push_back(hit);

            for (int i = erase_position.size() - 1; i >= 0; i--)
            {
                ix.erase(ix.begin() + erase_position[i]);
                iy.erase(iy.begin() + erase_position[i]);
                ichip_id.erase(ichip_id.begin() + erase_position[i]);
            }
            ix.erase(ix.begin());
            iy.erase(iy.begin());
            ichip_id.erase(ichip_id.begin());
        }
    }
}

*/ 

void Clusterer::VisualizeHits()
{

    for (int ichip = 0; ichip < fNChips; ichip++)
    {
        TH2D * h2 =  new TH2D(Form("chip%d", fChipIDs[ichip]), Form("Chip %d; x (position); y (position)", fChipIDs[ichip]), 1024, 0, 1024, 512, 0, 512);
        TH2D * h2hits = new TH2D(Form("chip%d_hits", fChipIDs[ichip]), Form("Chip %d; x (position); y (position)", fChipIDs[ichip]), 2048, 0, 1024, 1024, 0, 512);
        TH2D * h2cl1 = new TH2D(Form("chip%d_clsize1", fChipIDs[ichip]), Form("Chip %d - Cluster size = 1; x (position); y (position)", fChipIDs[ichip]), 1024, 0, 1024, 512, 0, 512);
        for (int ihit = 0; ihit < fHits.size(); ihit++)
        {
            if (fHits[ihit].chipID == fChipIDs[ichip])      
            {
                h2hits->Fill(fHits[ihit].x, fHits[ihit].y);
                if (fHits[ihit].clusterSize == 1)   h2cl1->Fill(fHits[ihit].x, fHits[ihit].y);
            }
        }

        double * chip_id = fInputData->getColumn("chip_id(decimal)");
        double * x = fInputData->getColumn("hit_x");
        double * y = fInputData->getColumn("hit_y");

        for (int irow = 0; irow < fInputData->getNRows(); irow++)
        {
            if (chip_id[irow] == fChipIDs[ichip])      h2->Fill(x[irow], y[irow]);
        }

        fOutFile->cd();
        h2->Write();
        h2hits->Write();
        h2cl1->Write();
        delete h2;
        delete h2hits;
        delete h2cl1;
    }

    TH2D * h2 = new TH2D("all_chips", "All Chips; x (position); y (position)", 2048, 0, 1024, 1024, 0, 512);
    TH2D * h2cl1 = new TH2D("all_chips_clsize1", "All Chips - Cluster size = 1; x (position); y (position)", 2048, 0, 1024, 1024, 0, 512);

    for (int ihit = 0; ihit < fHits.size(); ihit++)         
    {
        h2->Fill(fHits[ihit].x, fHits[ihit].y);
        if (fHits[ihit].clusterSize == 1)   h2cl1->Fill(fHits[ihit].x, fHits[ihit].y);
    }

    fOutFile->cd();
    h2->Write();
    h2cl1->Write();
    delete h2;
    delete h2cl1;
}

void Clusterer::VisualizeClusterSize()
{
    TH1D * h1 = new TH1D("clusterSize", "Cluster Size; Cluster Size; Counts", 10, 0, 10);
    for (int ihit = 0; ihit < fHits.size(); ihit++)         h1->Fill(fHits[ihit].clusterSize);

    fOutFile->cd();
    h1->Write();
    delete h1;
}