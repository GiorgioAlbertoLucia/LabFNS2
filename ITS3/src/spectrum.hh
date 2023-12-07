/*
    Class to study the spectrum of a radioactive isotope observed with a APTS sensor.
    The data is already available in the form of a ROOT file containing PixelData objects.
    The class is used to read the data from the ROOT file and to perform the analysis.
*/

#ifndef SPECTRUM_HH
#define SPECTRUM_HH

#include "../yaml/Yaml.hpp"

#include "pixelData.hh"

class Spectrum
{
    public:

        Spectrum(TH1D * inSpectrum, TFile * outFile);
        ~Spectrum();

        void EstimateBackground();
        void DrawSpectrumAndBackground() const;
        void SimpleFitSpectrum(const char * outLogPath, Yaml::Node & cfgFits);
        void FitSpectrum(const char * cfgPath, const int nFits = 1);
        void SaveSpectrum(const TString & name) const;

    private:
        

        // specifics for the spectrum histograms
        int fNBins;
        double fXmin;
        double fXmax;

        TH1D * fSpectrum;       // spectrum of the pixel signal
        TH1D * fBackground;     // spectrum of the background
        TH1D * fSignal;         // spectrum of the signal

        TFile * fOutFile;       // output file

        std::vector<TF1> fFits; // fits of the spectrum


};

#endif // SPECTRUM_HH