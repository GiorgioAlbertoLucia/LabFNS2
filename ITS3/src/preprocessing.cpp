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
    p.SetIgnorePoints(100);              // ignore the first 100 points of the signal
    p.SetNSample(100);                   // use 100 points to calculate the mean and the RMS
    p.SetNDerivativePoints(40);          // use 40 points to calculate the derivative
    p.SetNSmoothingPoints(10);           // use 10 points to calculate the smoothing of the waveform for the derivative

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

    p.UploadConversionValues(mV_to_electrons);
    p.BuildTree();
}
