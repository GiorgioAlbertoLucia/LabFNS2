#define wfm_cxx
#include "wfm.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>

#include <algorithm> // for std::max_element
#include <numeric> // for std::accumulate, std::inner_product

// progress bar section
// -------------------------------------------------------------------------------------
#include <Riostream.h>
#include <chrono>
#include <thread>

/**
 * @brief Function to update the progress bar.
 * 
 * @param progress 
 * @param total 
 */
void updateProgressBar(int progress, int total, const std::chrono::steady_clock::time_point& startTime) 
{
    const int barWidth = 50;

    float percentage = static_cast<float>(progress) / total;
    int barLength = static_cast<int>(percentage * barWidth);

    auto currentTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < barLength) {
            std::cout << "=";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << std::setw(3) << static_cast<int>(percentage * 100.0) << "%";
    std::cout << "  Elapsed Time: " << elapsedTime << "s";
    std::cout.flush();
}
// -------------------------------------------------------------------------------------

// functions used in this file
// -------------------------------------------------------------------------------------
double FitToA(const double *t, const double *w, const int npoints, const double ToA_guess);
std::pair<double, double> findMaxAndToA(double *waveform, int size, double *time);
std::pair<double, double> calculateBaselineAndRMS(double *waveform, int begin, int end, double to_mV);
// -------------------------------------------------------------------------------------


void wfm::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
    }
}

/**
 * @brief Function to build a new tree with the variables of interest.
 * The tree contains the following variables:
 * - Amp2: amplitude of the waveform of the second channel (mV)
 * - Amp3: amplitude of the waveform of the third channel (mV)
 * - ToA2: time of arrival of the waveform of the second channel (ns)
 * - ToA3: time of arrival of the waveform of the third channel (ns)
 * - ToA2f: time of arrival of the waveform of the second channel (ns) (fitting)
 * - ToA3f: time of arrival of the waveform of the third channel (ns) (fitting)
 * - RMS2: RMS of the waveform of the second channel (mV)
 * - RMS3: RMS of the waveform of the third channel (mV)
 * - baseline2: baseline of the waveform of the second channel (mV)
 * - baseline3: baseline of the waveform of the third channel (mV)
 * 
 * @param filename 
 */
void wfm::BuildTree(const char * filename = "Beta/data/output/BetaOutput_clean.root")
{
    TFile outFile(filename, "RECREATE");
    TTree outTree("BetaTree", "BetaTree");

    // define variables to be saved in the tree
    double Amp2{0.}, Amp3{0.};
    double ToA2{0.}, ToA3{0.};
    double ToA2f{0.}, ToA3f{0.};
    double RMS2{0.}, RMS3{0.};
    double baseline2{0.}, baseline3{0.};       

    double amplitude_cut2{50.}, amplitude_cut3{60.};

    // define branches of the tree
    outTree.Branch("Amp2", &Amp2);   
    outTree.Branch("Amp3", &Amp3); 
    outTree.Branch("RMS2", &RMS2); 
    outTree.Branch("RMS3", &RMS3); 
    outTree.Branch("ToA2", &ToA2); 
    outTree.Branch("ToA3", &ToA3);   
    outTree.Branch("ToA2f", &ToA2f); 
    outTree.Branch("ToA3f", &ToA3f);      
    outTree.Branch("baseline2", &baseline2); 
    outTree.Branch("baseline3", &baseline3); 

    const double to_mV{1e3}, to_ns{1e9};
    const int begin{20}, end{220};                  // range of samples to use for baseline and RMS
    const auto startTime = std::chrono::steady_clock::now();
    
    
    for (int ientry = 0; ientry < fChain->GetEntriesFast(); ientry++)
    {
        updateProgressBar(ientry, fChain->GetEntriesFast(), startTime);

        Long64_t entry = LoadTree(ientry);
        if (entry < 0) break;
        fChain->GetEntry(ientry);

        // find the maximum of the waveform   
        auto maxToA2 = findMaxAndToA(w2, nw2, t2);
        const double max2 = maxToA2.first;
        ToA2 = maxToA2.second;
        auto maxToA3 = findMaxAndToA(w3, nw3, t3);
        const double max3 = maxToA3.first;
        ToA3 = maxToA3.second;

        // find the baseline of the waveform
        auto baselineRMS2 = calculateBaselineAndRMS(w2, begin, end, to_mV);
        baseline2 = baselineRMS2.first;
        RMS2 = baselineRMS2.second;
        auto baselineRMS3 = calculateBaselineAndRMS(w3, begin, end, to_mV);
        baseline3 = baselineRMS3.first;
        RMS3 = baselineRMS3.second;

        // find the amplitude of the waveform
        Amp2 = (max2 - baseline2);                                          // mV
        Amp3 = (max3 - baseline3);                                          // mV

        // find the time of arrival of the waveform (fitting)
        ToA2f = (Amp2 > amplitude_cut2) ? FitToA(t2, w2, nt2, ToA2) : ToA2; // ns
        ToA3f = (Amp3 > amplitude_cut3) ? FitToA(t3, w3, nt3, ToA3) : ToA3; // ns

        // fill the tree
        outTree.Fill();
        
    }

    // write the tree to the file
    outTree.Write();
    outFile.Close();
}


// functions used in this file

/**
 * @brief Function to fit the waveform with a gaussian and find the time of arrival.
 * If the fit fails, the function returns the ToA_guess.
 * 
 * @param t 
 * @param w 
 * @param npoints 
 * @param ToA_guess 
 * @return double 
 */
double FitToA(const double *t, const double *w, const int npoints, const double ToA_guess)
{
    TH1D hist("hist", "hist", npoints, t[0], t[npoints - 1]);
    for (int i = 0; i < npoints; i++)
    {
        hist.SetBinContent(i + 1, w[i]);
        hist.SetBinError(i + 1, 6/sqrt(12));
    }
    TF1 fit("fit", "gaus", ToA_guess - 7, ToA_guess + 7);
    hist.Fit(&fit, "RMLQ+");

    // check NaN values 
    if (fit.GetChisquare() != fit.GetChisquare())
    {
        std::cout << "NaN value in chisquare" << std::endl;
        return 0;
    }
    if (std::isnan(fit.GetParameter(1)))
    {
        std::cout << "NaN value in ToA" << std::endl;
        return ToA_guess;
    }
    
    return fit.GetParameter(1);
}

std::pair<double, double> findMaxAndToA(double *waveform, int size, double *time)
{
    double max_val = waveform[0] * 1e3; // Initialize max_val to the first element
    double ToA = time[0] * 1e9;          // Initialize ToA to the first element

    for (int i = 1; i < size; ++i)
    {
        if (waveform[i] > max_val)
        {
            max_val = waveform[i] * 1e3; // Update max_val if a higher value is found
            ToA = time[i] * 1e9;          // Update ToA accordingly
        }
    }

    return {max_val, ToA};
}

std::pair<double, double> calculateBaselineAndRMS(double *waveform, int begin, int end, double to_mV)
{
    double sum = 0.0, sumsq = 0.0;

    for (int i = begin; i < end; ++i)
    {
        sum += waveform[i];
        sumsq += waveform[i] * waveform[i];
    }

    double mean = sum / (end - begin);
    double RMS = sqrt((sumsq - 2 * sum * mean + (end - begin) * mean * mean) / (end - begin)) * to_mV;

    return {mean * to_mV, RMS};
}