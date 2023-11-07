#include <Riostream.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLine.h>
#include <TAxis.h>

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
    std::pair<double, double> riseTimeExtremes = findRiseTimeExtremes();
    fRiseTime = riseTimeExtremes.second - riseTimeExtremes.first;
    
    return fRiseTime;
}

/**
 * @brief Draw waveform of the signal and print the rise time.
 * 
 */
void RiseTime::drawWaveform(const char * outputPath, const int nEvent)
{
    TCanvas c1("c1","c1",2100,1000);
    TGraph gr(fWaveform.size(), &fTimeFrame[0], &fWaveform[0]);
    gr.SetTitle(Form("Waveform - Event %d; Time (ns); Amplitude (mV)", nEvent));
    gr.SetLineColor(kBlue);
    gr.SetLineWidth(1);
    gr.SetMarkerColor(kBlue);
    gr.Draw("AL");

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.SetNDC();
    latex.DrawLatex(0.2, 0.8, "^{90}Sr source acquisition");
    latex.DrawLatex(0.2, 0.75, Form("Event %d", nEvent));
    latex.DrawLatex(0.2, 0.7, (Form("Rise time: from %d", int(fThreshold*100))+std::string("%")+Form(" to %d", int((1.0-fThreshold)*100))+std::string("%")).c_str());
    latex.DrawLatex(0.2, 0.65, Form("Rise Time = %.2f ns", fRiseTime));

    // add lines for rise time extremes
    std::pair<double, double> riseTimeExtremes = findRiseTimeExtremes();
    TLine line1(riseTimeExtremes.first, gr.GetYaxis()->GetXmin(), riseTimeExtremes.first, gr.GetYaxis()->GetXmax());
    TLine line2(riseTimeExtremes.second, gr.GetYaxis()->GetXmin(), riseTimeExtremes.second, gr.GetYaxis()->GetXmax());
    
    line1.SetLineColor(kRed);
    line1.SetLineWidth(1);
    line1.SetLineStyle(2);
    line1.Draw();
    line2.SetLineColor(kRed);
    line2.SetLineWidth(1);
    line2.SetLineStyle(2);
    line2.Draw();


    c1.SaveAs(outputPath);
}

/*
    PROTECTED METHODS
*/

/**
 * @brief Finds rise time extremes as the time at the x% and the 100-x% amplitude of the signal.
 * These points are found fitting the waveform with a line.
 * 
 * @param begin 
 * @return std::pair<int, int> 
 */
std::pair<double, double> RiseTime::findRiseTimeExtremes()
{
    // Find the beginning of the signal
    auto derivative = RiseTime::computeDerivative();
    
    auto max_der = std::max_element(derivative.begin(), derivative.end());
    auto begin = std::distance(derivative.begin(), max_der) * fWindow;  // Index of the maximum derivative in the waveform vector

    // find last point of the fit
    auto max = std::max_element(fWaveform.begin(), fWaveform.end());
    auto max_index = std::distance(fWaveform.begin(), max); 

    TGraph graph(fWaveform.size(), &fTimeFrame[0], &fWaveform[0]);
    auto fit = new TF1("fit", "[0]*x+[1]", fTimeFrame.at(begin), fTimeFrame.at(max_index));
    graph.Fit(fit, "RQC");

    fit->SetRange(fTimeFrame.at(0), fTimeFrame.at(fTimeFrame.size()-1));
    double firstExtreme = fit->GetX(fThreshold * *max);
    double lastExtreme = fit->GetX((1.0-fThreshold) * *max);

    delete fit;
    return std::make_pair(firstExtreme, lastExtreme);
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
