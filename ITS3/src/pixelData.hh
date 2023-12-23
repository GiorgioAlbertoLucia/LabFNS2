<<<<<<< HEAD
/*
    Struct to store all relevant information from a pixel.
    This struct is used to create a branch in the output tree.
*/

#ifndef PIXELDATA_HH
#define PIXELDATA_HH

#include <TString.h>

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

#endif // PIXELDATA_HH
=======
/*
    Struct to store all relevant information from a pixel.
    This struct is used to create a branch in the output tree.
*/

#ifndef PIXELDATA_HH
#define PIXELDATA_HH

#include <TString.h>

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
    double electrons;   // deposited charge in electrons
    double RMS;         // RMS of the waveform

    static TString GetBranchList()  { return "pixel/I:samplingPeriod/I:baseline/D:minLevel/D:t10/D:t90/D:t50/D:fallTime/D:amplitude/D:electrons/D:RMS/D"; }
};

#endif // PIXELDATA_HH
>>>>>>> 426078767d4e7cb224e63f70c50666d91a8e9623
