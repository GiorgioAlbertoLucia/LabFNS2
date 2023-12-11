/*
    Script to do the full preprocessing of a file
*/

#include <TString.h>

#include "preprocessor.hh"

#include "../yaml/Yaml.hpp"

void preprocessing()
{
    const char * inFilePath = "ITS3/Data/run175174828_230428174901.root";
    const char * conversionLogPath = "ITS3/Data/output/calibrationOutput.log";

    const double threshold = 1.;        // using 1 mV threshold
    Preprocessor p(inFilePath, threshold);
    p.SetIgnorePoints(100);              // ignore the first 25 points of the signal

    const int nPixels = p.GetNPixels();
    double mV_to_electrons[nPixels][2];
    for (int iPixel = 0; iPixel < nPixels; iPixel++)
    {
        mV_to_electrons[iPixel][0] = 0.;
        mV_to_electrons[iPixel][1] = 0.;
    }

    Yaml::Node conversionLog;
    Yaml::Parse(conversionLog, conversionLogPath);
    for (int iPixel = 0; iPixel < nPixels; iPixel++)
    {
        mV_to_electrons[iPixel][0] = conversionLog[Form("pixel%d", iPixel)]["calibrationParsElectrons"][0].As<double>();
        mV_to_electrons[iPixel][1] = conversionLog[Form("pixel%d", iPixel)]["calibrationParsElectrons"][1].As<double>();
    }

    //p.DrawEvent(1);

    //p.UploadConversionValues(mV_to_electrons);
    p.BuildTree();
}
