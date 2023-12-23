/*
    Script to do the full preprocessing of a file
*/

#include "preprocessor.hh"

void preprocessing()
{
    const char * inFilePath = "ITS3/Data/run175174828_230428174901.root";

    const double threshold = 1.;        // using 1 mV threshold
    Preprocessor p(inFilePath, threshold);
    p.BuildTree();
}
