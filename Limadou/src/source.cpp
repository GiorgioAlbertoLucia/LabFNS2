/*

*/

#include "lib/sourceAnalysis.hh"

void source()
{
    const int chipID = 113;
    const char * inFilePath = "Limadou/Data/output/source_123218.root";
    const char * cfgFilePath = "Limadou/src/config/cfgSourcePosition.yml";

    //const char * bkgFilePath1 = "Limadou/Data/output/background_113112.root";
    const char * bkgFilePath2 = "Limadou/Data/output/background_115654.root";

    SourceAnalysis sa(chipID, inFilePath);
    //sa.SubtractBackground(bkgFilePath1);    
    sa.SubtractBackground(bkgFilePath2);
    sa.Save();
    sa.FitHits(cfgFilePath);


}