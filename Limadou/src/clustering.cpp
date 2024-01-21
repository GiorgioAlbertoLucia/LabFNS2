#include "lib/clusterer.hh"

void clustering()
{

    //const char * inFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/Gruppo1_2023/20231130_113112_DAQ/data.txt";
    //const char * inFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/Gruppo1_2023/20231130_115654_DAQ/data.txt";
    //const char * inFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/Gruppo1_2023/20231130_123218_DAQ/data.txt";
    const char * inFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/Gruppo1_2023/20231130_130551_DAQ/data.txt";

    //const char * outFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/output/background_113112.root";
    //const char * outFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/output/background_115654.root";
    //const char * outFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/output/source_123218.root";
    const char * outFile = "/Users/glucia/Projects/LabFNS2/Limadou/Data/output/source_130551.root";

    Clusterer clusterer(inFile, outFile);
    clusterer.FindHits();
    clusterer.VisualizeHits();
    clusterer.VisualizeClusterSize();
}