#ifndef PREPROCESS_HH
#define PREPROCESS_HH

#include <TString.h>

/**
 * @brief 
 * 
 * 
 */
struct PixelData
{
    int channel;
    double baseline;    // baseline value
    double minLevel;    // minimum level of the waveform (basically baseline after the signal)
    double t10;         // time of arrival at 10% of the amplitude
    double t90;         // time of arrival at 90% of the amplitude
    double t50;         // time of arrival at 50% of the amplitude
    double fallTime;    // fall time of the waveform
    double amplitude;   // amplitude of the waveform
    double RMS;         // RMS of the waveform

    static TString GetBranchList()  { return "channel/I:baseline/D:minLevel/D:t10/D:t90/D:t50/D:fallTime/D:amplitude/D:RMS/D"; }
};

class Preprocessor
{
    public:
        Preprocessor(const char * inFilePath, const double threshold = 10);
        ~Preprocessor();

        TString GetInFilePath() const { return fInFilePath; }
        int GetNChannels() const { return fNChannels; }
        int GetNEvents() const { return fNEvents; }
        int GetSamplingPeriod() const { return fSamplingPeriod; }

        void ReadInput();
        void BuildTree(const char * outFilePath = "default");
        bool ProcessEvent(const int event, const int channel, PixelData& channelData);
        void DrawEvent(const int event, const int channel, const char * outFilePath = "default");

    private:

        TString fInFilePath;

        int fNChannels;
        int fNEvents;
        int fSamplingPeriod;    // in ps
        
        double fThreshold;      // in mV
};

#endif 