/*
    This is the header file for the preprocessor class.
    This class is used to read waveforms in input from a .root file and to process them.
    The waveforms are given as TGraph objects, one for each pixel.
    The preprocessor class is used to extract the following information from the waveforms:
        - baseline
        - minimum level
        - time of arrival at 10% of the amplitude
        - time of arrival at 90% of the amplitude
        - time of arrival at 50% of the amplitude
        - fall time
        - amplitude
        - RMS

    The preprocessor class is also used to build a new .root file containing the processed waveforms.
    The new .root file contains a TTree object with one entry for each pixel.
    The TTree object contains the following branches:
        nPixel PixelData object
        where PixelData is a struct containing the following information:
        - pixel: pixel number
        - baseline: baseline value
        - minLevel: minimum level of the waveform (basically baseline after the signal)
        - t10: time of arrival at 10% of the amplitude
        - t90: time of arrival at 90% of the amplitude
        - t50: time of arrival at 50% of the amplitude
        - fallTime: fall time of the waveform
        - amplitude: amplitude of the waveform
        - RMS: RMS of the waveform

    NOTE: The names of the TGraph in the input file are used to extract infomation from the file itself.
    In previous versions the name was "grEvXXXChanCYYYsampZZZ", where XXX is the event number, YYY is 
    the channel number and ZZZ is the sampling period in ps. This class is currently expecting 
    "grEvXXXPxYYYsampZZZ", where XXX is the event number, YYY is the pixel number and ZZZ is the
    sampling period in ps. For further changes in the TGraph input name please modifify the functions
    ReadInput() and ProcessEvent().
*/

#ifndef PREPROCESS_HH
#define PREPROCESS_HH

#include <TString.h>
#include <TFile.h>

/**
 * @brief 
 * 
 * 
 */
struct PixelData
{
    int pixel;          // pixel number
    int samplingPeriod; // sampling period in ps
    double baseline;    // baseline value
    double minLevel;    // minimum level of the waveform (basically baseline after the signal)
    double t10;         // time of arrival at 10% of the amplitude
    double t90;         // time of arrival at 90% of the amplitude
    double t50;         // time of arrival at 50% of the amplitude
    double fallTime;    // fall time of the waveform
    double amplitude;   // amplitude of the waveform
    double RMS;         // RMS of the waveform

    static TString GetBranchList()  { return "pixel/I:samplingPeriod/I:baseline/D:minLevel/D:t10/D:t90/D:t50/D:fallTime/D:amplitude/D:RMS/D"; }
};

class Preprocessor
{
    public:
        Preprocessor(const char * inFilePath, const double threshold = 10);
        ~Preprocessor();

        TString GetInFilePath() const { return fInFilePath; }
        int GetNPixels() const { return fNPixels; }
        int GetNEvents() const { return fNEvents; }
        int GetSamplingPeriod(const int pixel) const { return fSamplingPeriodDictionary[pixel]; }

        
        void BuildTree(const char * outFilePath = "default");
        bool ProcessEventScope(const int event, const int pixel, PixelData & pixelData, TFile * inFile);
        bool ProcessEventADC(const int event, const int pixel, PixelData & pixelData, TFile * inFile);
        void DrawEvent(const int event, const int pixel, const char * outFilePath = "default");

    protected:
        void ReadInput();
        void GenerateSamplingDictionary();

    private:

        TString fInFilePath;

        int fNPixels;
        int fNEvents;

        int * fSamplingPeriodDictionary;    // in ps
        
        double fThreshold;                  // in mV
};

#endif 