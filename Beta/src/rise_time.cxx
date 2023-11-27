#include <Riostream.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLine.h>
#include <TAxis.h>
#include <TH1D.h>

#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
#include <TLegend.h>

#include <algorithm>
#include <vector>

#include "rise_time.hh"

/*  USEFUL FUNCTIONS    */
std::vector<double> Smooth(std::vector<double>& y, const int nSmoothingPoints);
std::vector<double> DerivativeVector(std::vector<double>& y, std::vector<double>& x);
std::pair<double, double> GetMeanAndRMS(std::vector<double>& v, const int& begin, const int& end);
int FindChangingDerivative(std::vector<double>& yDerivative, const double& meanDer, const double& RMSDer, const int ignorePoints, const bool backwards = false);
int FindEdgeIndex(std::vector<double>& y, std::vector<double>& x, const char * whichEdge);




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
double RiseTime::findRiseTime(const double minAmplitude)
{
    std::pair<double, double> riseTimeExtremes = findRiseTimeExtremes(minAmplitude);
    fRiseTime = riseTimeExtremes.second - riseTimeExtremes.first;
    
    return fRiseTime;
}

/**
 * @brief Draw waveform of the signal and print the rise time.
 * 
 */
void RiseTime::drawWaveform(const char * outputPath, const int nEvent, const double minAmplitude)
{
    TCanvas c1("c1","c1",2100,1000);
    TGraph gr(fWaveform.size(), &fTimeFrame[0], &fWaveform[0]);
    gr.SetTitle(Form("Waveform - Event %d; Time (ns); Amplitude (mV)", nEvent));
    gr.SetLineColor(kBlue);
    gr.SetLineWidth(1);
    gr.SetMarkerColor(kBlue);
    gr.Draw("AL");

    // add lines for rise time extremes
    std::pair<double, double> riseTimeExtremes = findRiseTimeExtremes(minAmplitude);
    TLine line1(riseTimeExtremes.first, gr.GetYaxis()->GetXmin(), riseTimeExtremes.first, gr.GetYaxis()->GetXmax());
    TLine line2(riseTimeExtremes.second, gr.GetYaxis()->GetXmin(), riseTimeExtremes.second, gr.GetYaxis()->GetXmax());

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.SetNDC();
    latex.DrawLatex(0.2, 0.8, "^{90}Sr source acquisition");
    latex.DrawLatex(0.2, 0.75, Form("Event %d", nEvent));
    latex.DrawLatex(0.2, 0.7, (Form("Rise time: from %d", int(fThreshold*100))+std::string("%")+Form(" to %d", int((1.0-fThreshold)*100))+std::string("%")).c_str());
    latex.DrawLatex(0.2, 0.65, Form("Rise Time = %.2f ns", fRiseTime));
    
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

void RiseTime::drawLinearFit(const char * outputPath, const int nEvent, const double minAmplitude)
{
    gStyle->SetOptFit(0);

    TCanvas c1("c1","c1",2100,1000);
    TH1D hist("hist", Form("Waveform - Event %d; Time (ns); Amplitude (mV)", nEvent), fWaveform.size(), fTimeFrame.at(0), fTimeFrame.at(fTimeFrame.size()-1));
    hist.SetLineColor(kBlue);
    hist.SetLineWidth(1);

    for (unsigned long i = 0; i < fWaveform.size(); ++i)
    {
        hist.SetBinContent(i+1, fWaveform.at(i));
        hist.SetBinError(i+1, sqrt(fWaveform.at(i)));
    }

    // add lines for rise time extremes
    std::pair<double, double> riseTimeExtremes = findRiseTimeExtremes(minAmplitude);
    hist.GetXaxis()->SetRangeUser(riseTimeExtremes.first-5, riseTimeExtremes.second+5);

    // Find the beginning of the signal
    //auto derivative = RiseTime::computeDerivative();
    //auto max_der = std::max_element(derivative.begin(), derivative.end());
    //auto begin = std::distance(derivative.begin(), max_der) * fWindow;  // Index of the maximum derivative in the waveform vector
    auto begin = FindEdgeIndex(fWaveform, fTimeFrame, "left");

    // find last point of the fit
    auto max = std::max_element(fWaveform.begin(), fWaveform.end());
    auto max_index = std::distance(fWaveform.begin(), max);

    auto fit = new TF1("fit", "[0]*x+[1]", fTimeFrame.at(begin), fTimeFrame.at(max_index));
    fit->SetLineColor(kRed);
    fit->SetLineWidth(1);
    hist.Fit(fit, "RQC");

    // insert two orange points on the fit where the rise time extremes are
    TGraph graphExtremes(2);
    graphExtremes.SetPoint(0, riseTimeExtremes.first, fit->Eval(riseTimeExtremes.first));
    graphExtremes.SetPoint(1, riseTimeExtremes.second, fit->Eval(riseTimeExtremes.second));
    graphExtremes.SetMarkerStyle(20);
    graphExtremes.SetMarkerColor(kOrange);

    c1.cd();
    hist.Draw("hist");
    fit->Draw("same");
    graphExtremes.Draw("same p");

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.SetNDC();
    latex.DrawLatex(0.2, 0.8, "^{90}Sr source acquisition");
    latex.DrawLatex(0.2, 0.75, Form("Event %d", nEvent));
    latex.DrawLatex(0.2, 0.7, (Form("Rise time: from %d", int(fThreshold*100))+std::string("%")+Form(" to %d", int((1.0-fThreshold)*100))+std::string("%")).c_str());
    latex.DrawLatex(0.2, 0.65, Form("Rise Time = %.2f ns", fRiseTime));

    TLegend leg(0.2, 0.4, 0.5, 0.55);
    leg.AddEntry(&hist, "Waveform", "lf");
    leg.AddEntry(fit, "Linear fit", "l");
    leg.AddEntry(&graphExtremes, "Rise time extremes", "p");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.04);
    leg.SetFillStyle(0);
    leg.Draw("same");

    c1.SaveAs(outputPath);

    delete fit;
}

void RiseTime::drawDerivative(const char * outputPath, const int nEvent, const double minAmplitude)
{
    TCanvas c1("c1","c1",2100,1000);
    auto derivativeV = DerivativeVector(fWaveform, fTimeFrame);
    TGraph gr(fWaveform.size(), &fTimeFrame[0], &derivativeV[0]);
    gr.SetTitle(Form("Waveform derivative - Event %d; Time (ns); Derivative (mV ns^{-1})", nEvent));
    gr.SetMarkerStyle(20);
    gr.SetMarkerColor(kBlue);
    gr.SetLineColor(kBlue);
    gr.SetLineWidth(1);
    gr.GetXaxis()->SetRangeUser(-5, 5);
    gr.Draw("AL");

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.SetNDC();
    latex.DrawLatex(0.2, 0.75, "^{90}Sr source acquisition");
    latex.DrawLatex(0.2, 0.7, Form("Event %d", nEvent));

    // add graph with rise time extremes
    std::pair<double, double> riseTimeExtremes = findRiseTimeExtremes(minAmplitude);
    TGraph graphExtremes(2);
    graphExtremes.SetPoint(0, riseTimeExtremes.first, gr.GetYaxis()->GetXmin());
    graphExtremes.SetPoint(1, riseTimeExtremes.first, gr.GetYaxis()->GetXmax());
    graphExtremes.SetMarkerStyle(20);
    graphExtremes.SetMarkerColor(kOrange);
    graphExtremes.Draw("same p");

    TLegend leg(0.6, 0.6, 0.9, 0.8);
    leg.AddEntry(&gr, "Waveform derivative", "lf");
    leg.AddEntry(&graphExtremes, "Rise time extremes", "p");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.04);
    leg.SetFillStyle(0);
    leg.Draw("same");

    c1.SaveAs(outputPath);
}

/*
    PROTECTED METHODS
*/

/**
 * @brief Finds rise time extremes as the time at the x% and the 100-x% amplitude of the signal.
 * These points are found fitting the waveform with a line.
 * 
 * @param minAmplitude Minimum amplitude of the signal to be considered. If the maximum amplitude is below this value, the rise time is set to 1 s.
 * @return std::pair<int, int> 
 */
std::pair<double, double> RiseTime::findRiseTimeExtremes(const double minAmplitude)
{
    // Find the beginning of the signal
    int begin = FindEdgeIndex(fWaveform, fTimeFrame, "left");
    
    //auto derivative = RiseTime::computeDerivative();
    //auto max_der = std::max_element(derivative.begin(), derivative.end());
    //auto begin = std::distance(derivative.begin(), max_der) * fWindow;  // Index of the maximum derivative in the waveform vector

    //auto zero = std::find(fTimeFrame.begin(), fTimeFrame.end(), 0.);
    //auto begin = std::distance(fTimeFrame.begin(), zero);
    //const int begin = fTimeFrame.size() / 2;

    // find last point of the fit
    auto max = std::max_element(fWaveform.begin(), fWaveform.end());
    auto max_index = std::distance(fWaveform.begin(), max);

    // return 1 s if the max amplitude is below the minimum amplitude
    if(*max < minAmplitude) return std::make_pair(0, 1e9);

    TGraph graph(fWaveform.size(), &fTimeFrame[0], &fWaveform[0]);
    auto fit = new TF1("fit", "[0]*x+[1]", fTimeFrame.at(begin), fTimeFrame.at(max_index));
    graph.Fit(fit, "RQC");

    fit->SetRange(fTimeFrame.at(0), fTimeFrame.at(fTimeFrame.size()-1));
    double baseline = 0;
    for(int i = 20; i < 220; ++i)   baseline += fWaveform.at(i)/200;
    double firstExtreme = fit->GetX(fThreshold * (*max - baseline) + baseline);
    double lastExtreme = fit->GetX((1.0-fThreshold) * (*max - baseline) + baseline);

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


/*  USEFUL FUNCTIONS    */


/**
 * @brief Function to find the edge of the signal.
 * 
 * @param gr 
 * @return std::tuple<double, double, double> edgeLeft, edgeRight
 */
int FindEdgeIndex(std::vector<double>& y, std::vector<double>& x, const char * whichEdge = "left")
{
    const int nSmoothingPoints = 5;                         // number of points to consider for the smoothing
    auto ySmooth = Smooth(y, nSmoothingPoints);             // smooth the graph to find the edge more easily
    auto yDerivative = DerivativeVector(ySmooth, x);           // take the derivative of the smoothed graph
    
    const int nSamples = 80;                                // number of samples to consider for the plateau
    double meanDer{0.}, RMSDer{0.};
    auto resultDer = GetMeanAndRMS(yDerivative, nSmoothingPoints, nSmoothingPoints+nSamples);
    meanDer = resultDer.first;
    RMSDer = resultDer.second;

    if (std::string(whichEdge) == "left")       return FindChangingDerivative(yDerivative, meanDer, RMSDer, nSmoothingPoints);        
    else if (std::string(whichEdge) == "right") return FindChangingDerivative(yDerivative, meanDer, RMSDer, nSmoothingPoints, true);   
    else
    {
        std::cout << "Error: whichEdge must be either \"left\" or \"right\"." << std::endl;
        return 0.;
    }
}

/**
 * @brief Find the pont in which the derivative changes for three consecutive points with a value higher than meanDer + 3 * RMSDer.
 * 
 * @param grDer 
 * @param meanDer 
 * @param RMSDer 
 * @param ignorePoints number of points to ignore at the beginning and at the end of the TGraph (due to the smoothing)
 * @param backwards if true, the function loops through the TGraph from the last point to the first one
 * @return const int 
 */
int FindChangingDerivative(std::vector<double>& yDerivative, const double& meanDer, const double& RMSDer, const int ignorePoints, const bool backwards)
{
    const int nSamples = yDerivative.size();
    int changingPoint{0};
    
    if (!backwards)
    {
        for (int i = ignorePoints; i < nSamples - 2; ++i)
        {
            // positive derivative
            if (yDerivative[i] > meanDer + 3 * RMSDer && yDerivative[i + 1] > meanDer + 3 * RMSDer && yDerivative[i + 2] > meanDer + 3 * RMSDer)
            {
                changingPoint = i;
                break;
            }
            // negative derivative
            if (yDerivative[i] < meanDer - 3 * RMSDer && yDerivative[i + 1] < meanDer - 3 * RMSDer && yDerivative[i + 2] < meanDer - 3 * RMSDer)
            {
                changingPoint = i;
                break;
            }
        }
    }
    else
    {
        for (int i = nSamples - ignorePoints - 1; i > 1; --i)
        {
            // positive derivative
            if (yDerivative[i] > meanDer + 3 * RMSDer && yDerivative[i - 1] > meanDer + 3 * RMSDer && yDerivative[i - 2] > meanDer + 3 * RMSDer)
            {
                changingPoint = i;
                break;
            }
            // negative derivative  
            if (yDerivative[i] < meanDer - 3 * RMSDer && yDerivative[i - 1] < meanDer - 3 * RMSDer && yDerivative[i - 2] < meanDer - 3 * RMSDer)
            {
                changingPoint = i;
                break;
            }
        }
    }

    return changingPoint;
}

/**
 * @brief Get the Mean And RMS of a graph in a given range (points).
 * 
 * @param gr 
 * @param begin 
 * @param end 
 * @return std::pair<double, double> 
 */
std::pair<double, double> GetMeanAndRMS(std::vector<double>& v, const int& begin, const int& end)
{
    const int nSamples = v.size();

    double mean{0.};
    for (int i = begin; i < end; ++i)
    {
        mean += v[i];
    }
    mean /= (end - begin);

    double RMS{0.};
    for (int i = begin; i < end; ++i)
    {
        RMS += (v[i] - mean) * (v[i] - mean);
    }
    RMS /= (end - begin);
    RMS = std::sqrt(RMS);

    return std::make_pair(mean, RMS);
}

/**
 * @brief Function to obtain a smoothed version of a TGraph, by doing a running average over nSmoothingPoints points.
 * 
 * @param gr 
 * @param nSmoothingPoints 
 * @return TGraph* 
 */
std::vector<double> Smooth(std::vector<double>& y, const int nSmoothingPoints)
{

    const int nSamples = y.size();
    std::vector<double> ySmooth;
    ySmooth.reserve(nSamples);

    for (int i = 0; i < nSamples; ++i)
    {
        double sum{0.};
        for (int j = i - nSmoothingPoints; j < i + nSmoothingPoints; ++j)
        {
            if (j < 0 || j > nSamples) continue;
            sum += y[j];
        }
        ySmooth.push_back(sum / (2 * nSmoothingPoints + 1));
    }

    return ySmooth;
}

/**
 * @brief Function to obtain the derivative of a TGraph point by point.
 * 
 * @param gr 
 * @return TGraph* 
 */
std::vector<double> DerivativeVector(std::vector<double>& y, std::vector<double>& x)
{
    const int nSamples = y.size();
    std::vector<double> yDerivative;
    yDerivative.reserve(nSamples);

    for (int i = 0; i < nSamples; ++i)
    {
        double der{0.};
        for (int j = i - 1; j < i + 1; ++j)
        {
            if (j < 0 || j > nSamples) continue;
            der = (y[j + 1] - y[j - 1]) / (x[j + 1] - x[j - 1]);
        }
        yDerivative.push_back(der);
    }

    // set the first and last point of the derivative graph to the second and second to last point of the original graph
    yDerivative[0] = yDerivative[1];
    yDerivative[nSamples - 1] = yDerivative[nSamples - 2];

    return yDerivative;
}
