#ifndef RISETIME_HH
#define RISETIME_HH

#include <Riostream.h>

#include <vector>
#include <algorithm>

/**
 * @brief Class to find the rise time of a signal. The rise time is defined as the time 
 *       between the x% and the 100-x% amplitude of the signal. x is the threshold.
 * 
 * @param threshold x% for the rise time calculation (number between 0 and 1)
 * @param window Window size for the derivative calculation
 * 
 */
class RiseTime
{
    public:
        RiseTime(double threshold, int window, std::vector<double> const & waveform, std::vector<double> const & timeFrame);
        ~RiseTime();

        double GetRiseTime() { return fRiseTime; }
        
        double findRiseTime();
        void drawWaveform(const char * outputPath, const int nEvent);

        void reset() { fRiseTime = 0.; fWaveform.clear(); fTimeFrame.clear(); };

    protected:
        std::vector<double> computeDerivative();
        std::pair<int, int> findRiseTimeIndices(const int begin);

    private:
        
        double fRiseTime;
        double fThreshold;                  // Threshold for rise time calculation in percent of max amplitude
        int fWindow;                        // Window for rise time calculation in samples (in timeframes)

        std::vector<double> fWaveform;      // Vector of the waveform
        std::vector<double> fTimeFrame;     // Vector of the timeframes

};


#endif