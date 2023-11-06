#include <Riostream.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <TLatex.h>

#include <algorithm>
#include <vector>

#include "rise_time.hh"


/**
 * @brief Constructs a RiseTime object with the given threshold and window size.
 * 
 * @param threshold The threshold used to find the rise time. It is a number between 0 and 1.
 * @param window The window size used to compute the derivative.
 * @param waveform The waveform to be analyzed.
 */
RiseTime::RiseTime(double threshold, int window, std::vector<double> const & waveform, std::vector<double> const & timeFrame):
    fWindow(window),
    fWaveform(waveform),
    fTimeFrame(timeFrame)
{
    if(threshold < 0 || threshold > 1)  throw std::invalid_argument("Threshold must be a value between 0 and 1.");
    fThreshold = threshold;
}

RiseTime::~RiseTime()
{

}

/*
    PUBLIC METHODS
*/

/**
 * @brief Find rise time searching for the signal peak using a window derivative.
 * 
 * @param derivative 
 * @return int 
 */
double RiseTime::findRiseTime()
{
    // Find the beginning of the signal
    auto derivative = RiseTime::computeDerivative();
    
    auto max = std::max_element(derivative.begin(), derivative.end());
    auto max_index = std::distance(derivative.begin(), max) * fWindow;  // Index of the maximum derivative in the waveform vector

    std::pair<int, int> riseTimeIndices = findRiseTimeIndices(max_index);
    fRiseTime = fTimeFrame.at(riseTimeIndices.second) - fTimeFrame.at(riseTimeIndices.first);
    
    return fRiseTime;
}

/**
 * @brief Draw waveform of the signal and print the rise time.
 * 
 */
void RiseTime::drawWaveform(const char * outputPath, const int nEvent)
{
    TCanvas *c1 = new TCanvas("c1","c1",2100,1000);
    TGraph *gr = new TGraph(fWaveform.size(), &fTimeFrame[0], &fWaveform[0]);
    gr->SetTitle(Form("Waveform - Event %d; Time (ns); Amplitude (mV)", nEvent));
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kBlue);
    gr->Draw("AL");

    auto latex = new TLatex();
    latex->SetTextSize(0.04);
    latex->SetTextFont(42);
    latex->SetNDC();
    latex->DrawLatex(0.2, 0.8, "^{90}Sr source acquisition");
    latex->DrawLatex(0.2, 0.75, Form("Event %d", nEvent));
    latex->DrawLatex(0.2, 0.7, Form("Rise time: %.1f to %.1f %", fThreshold*100, (1.0-fThreshold)*100));
    latex->DrawLatex(0.2, 0.65, Form("Rise Time = %f ns", fRiseTime));

    c1->SaveAs(outputPath);
}

/*
    PROTECTED METHODS
*/

std::pair<int, int> RiseTime::findRiseTimeIndices(const int begin)
{
    auto max = std::max_element(fWaveform.begin(), fWaveform.end());
    
    int firstIndex(0), lastIndex(0);
    for(unsigned long i = begin; i < fWaveform.size(); i++)
    {
        if(fWaveform[i] > fThreshold * (*max))
        {
            firstIndex = i;
            break;
        }
    }
    for(unsigned long i = firstIndex; i < fWaveform.size(); i++)
    {
        if(fWaveform[i] > (1.0 - fThreshold) * (*max))
        {
            lastIndex = i;
            break;
        }
    }

    return std::make_pair(firstIndex, lastIndex);
}

/**
 * @brief Computes the numeric derivative of a vector with given window size. Returns the vector of derivatives.
 * 
 * @return std::vector<double> 
 */
std::vector<double> RiseTime::computeDerivative()
{
    std::vector<double> derivative;
    derivative.reserve(int(fWaveform.size()/fWindow) + 1);

    for(unsigned long i = 0; i < fWaveform.size(); i+=fWindow)
    {
        if(i+fWindow > fWaveform.size())    derivative.push_back(0);
        else                                derivative.push_back((fWaveform[i+fWindow]-fWaveform[i])/fWindow);
    }
    
    return derivative;
}
