// Implementation of Preprocessor class

#include "preprocessor.hh"
#include "colorTerminal.hh"

#include <iostream>
#include <fstream>
#include <streambuf>
#include <numeric>
#include <algorithm>

#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TKey.h>
#include <TTree.h>
#include <TList.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>


// -------------------------------------------------------------------------------------
// functions used in this file
std::pair<double, double> FindEdge(TGraph & gr, int nIgnorePoints = 0, int nChecks = 10, int nSample = 200, const int nDerivativePoints = 40, const int nSmoothingPoints = 10);
std::pair<int, int> FindEdgeIndex(TGraph & gr, int nIgnorePoints = 0, int nChecks = 10, int nSample = 200, const int nDerivativePoints = 40, const int nSmoothingPoints = 10);
int FindMinimumIndex(TGraph & gr);
int FindMaximumIndex(TGraph & gr);  
int FindChangingDerivative(TGraph & grDer, const double& meanDer, const double& RMSDer, const int ignorePoints, const int nChecks, const bool backwards = false);
std::pair<double, double> GetMeanAndRMS(TGraph & gr, const int& begin, const int& end);
TGraph* SmoothGraph(TGraph& gr, const int nSmoothingPoints = 2);
TGraph* DerivativeGraph(TGraph& gr, const int nDerivativePoints = 40);
int FindPointIndex(TGraph& gr, const double& x, const int nIgnorePoints = 0);
int FindClosestPointIndex(TGraph& gr, const double y, const int nIgnorePoints = 0);

// -------------------------------------------------------------------------------------


// progress bar section
// -------------------------------------------------------------------------------------
#include <chrono>
#include <thread>
void updateProgressBar(int progress, int total, const std::chrono::steady_clock::time_point& startTime);
// -------------------------------------------------------------------------------------


Preprocessor::Preprocessor(const char * inFilePath, const double threshold): 
    fInFilePath(inFilePath), 
    fNPixels(0), 
    fNEvents(0), 
    fSamplingPeriodDictionary(NULL), 
    fThreshold(threshold), 
    fIgnorePoints(0),
    fChecks(10),
    fNSample(200),
    fNDerivativePoints(40),
    fNSmoothingPoints(10)
{
    Preprocessor::ReadInput();   
    Preprocessor::GenerateSamplingDictionary();

    fmVToElectrons = new double*[fNPixels];
    for (int i = 0; i < fNPixels; i++)
    {
        fmVToElectrons[i] = new double[2];
        fmVToElectrons[i][0] = 0.;
        fmVToElectrons[i][1] = 0.;
    }
}

Preprocessor::~Preprocessor()
{
    delete fSamplingPeriodDictionary;
    delete fmVToElectrons;
}

/*  PROTECTED */

/**
 * @brief Loops throught the whole file and finds information such as number of channels and number of events.
 * 
 */
void Preprocessor::ReadInput()
{
    auto inFile = TFile::Open(fInFilePath.Data());
    const int nKeys = inFile->GetNkeys();
    TIter next(inFile->GetListOfKeys());

    int event{0}, pixel{0}, samplingPeriod{0}, index{0};
    const auto startTime = std::chrono::steady_clock::now();
    TKey *key;

    std::cout << "Reading input file: " << BLUE << UNDERLINE << fInFilePath << RESET << std::endl;
    while ((key = (TKey*)next())) 
    {
        updateProgressBar(index, nKeys, startTime);

        if (strncmp(key->GetClassName(), "TGraph", 6) == 0)  
        {
            sscanf(key->GetName(), "grEv%dPx%dsamp%d", &event, &pixel, &samplingPeriod);
            fNEvents = (event > fNEvents) ? event : fNEvents;
            fNPixels = (pixel > fNPixels) ? pixel : fNPixels;
        }

        index++;
    }
    std::cout << std::endl;

    fNPixels++;
}

/**
 * @brief Function to generate a dictionary which associates the sampling period to each pixel index.
 * Assume that the correct information can be found in event 0.
 * 
 */
void Preprocessor::GenerateSamplingDictionary()
{
    auto inFile = TFile::Open(fInFilePath.Data());
    auto listKeys = inFile->GetListOfKeys();
    long int nKeys = inFile->GetNkeys();

    int event{0}, pixel{0}, samplingPeriod{0};
    fSamplingPeriodDictionary = new int[fNPixels];

    std::cout << "Generating sampling period dictionary... " << std::endl;

    for (long int ipixel = 0; ipixel < fNPixels; ipixel++) 
    {
        fSamplingPeriodDictionary[ipixel] = 0;

        for (long int i = 0; i < nKeys; i++)             // get the information looping through the names of the objects in the file
        {   
            TKey* key = (TKey*)listKeys->At(i);             // get the ith object
            if(!key) 
            {
                std::cerr << "\tError: key " << i << " not found." << std::endl;
                continue;              
            }            
            TString className = key->GetClassName();        // get the type of the object
            TString objectName = key->GetName();            // get the name of the object

            if (className == "TGraph")  
            {
                sscanf(objectName.Data(), "grEv%dPx%dsamp%d", &event, &pixel, &samplingPeriod);
                if (event == 0 && pixel == ipixel) 
                {
                    fSamplingPeriodDictionary[pixel] = samplingPeriod;
                    break;    
                }
            }
        }
    }

    inFile->Close();
}


/*  PUBLIC  */

/**
 * @brief Upload a dictionary to convert mV to electrons.
 * 
 * @param mV_to_electrons 
 */
void Preprocessor::UploadConversionValues(double (* mVToElectrons)[2])
{
    for (int i = 0; i < fNPixels; i++)
    {
        fmVToElectrons[i][0] = mVToElectrons[i][0];
        fmVToElectrons[i][1] = mVToElectrons[i][1];
    }
}

/**
 * @brief Function to build the tree with the preprocessed data. All events in the input files will be processed.
 * 
 * @param outFilePath 
 */
void Preprocessor::BuildTree(const char * outFilePath)
{
    TString sOutFilePath(outFilePath);
    if (sOutFilePath.EqualTo("default")) 
    {
        sOutFilePath = fInFilePath;
        sOutFilePath.ReplaceAll(".root", "_preprocessed.root");
    }

    int nPixels = fNPixels;
    PixelData pixelData[nPixels];

    auto inFile = TFile::Open(fInFilePath.Data());
    TFile outFile(sOutFilePath.Data(), "RECREATE");
    TTree outTree("PreprocessedData", "outTree");

    int event{0};
    outTree.Branch("Event", &event, "Event/I");
    outTree.Branch("nPixels", &nPixels, "nPixels/I");
    for(int i = 0; i < fNPixels; i++) outTree.Branch(Form("pixel%d", i), &pixelData[i]);//, PixelData::GetBranchList().Data());
    
    std::cout << "Writing output file: " << BLUE << UNDERLINE << sOutFilePath << RESET << std::endl;
    const auto startTime = std::chrono::steady_clock::now();

    for (int ievent = 0; ievent < fNEvents; ++ievent)
    {
        updateProgressBar(ievent, fNEvents, startTime);
        
        event = ievent;
        for (int ipixel = 0; ipixel < fNPixels; ipixel++)   
        {
            if (fSamplingPeriodDictionary[ipixel] == 25)            Preprocessor::ProcessEventScope(ievent, ipixel, pixelData[ipixel], inFile);
            else if (fSamplingPeriodDictionary[ipixel] == 250000)   Preprocessor::ProcessEventADC(ievent, ipixel, pixelData[ipixel], inFile);
        }

        outTree.Fill();
    }

    inFile->Close();

    std::cout << std::endl;
    outFile.cd();
    outTree.Write();
    outFile.Close();
}

/**
 * @brief Function to process a single event read by the oscilloscope. 
 * It reads the graph from the input file and computes the variables of interest.
 * If no signal is found (i.e. the signal amplitude is below the set threshold) the time variables are set to -9999.
 * If a graph is not found, the function returns false.
 * 
 * @param event 
 * @param pixel 
 * @param pixelData 
 * @return true 
 * @return false 
 */
bool Preprocessor::ProcessEventScope(const int event, const int pixel, PixelData & pixelData, TFile * inFile)
{
    pixelData.pixel = pixel;  
    pixelData.samplingPeriod = fSamplingPeriodDictionary[pixel];
    pixelData.baseline = -9999.;
    pixelData.minLevel = -9999.;
    pixelData.t10 = -9999.;
    pixelData.t90 = -9999.;
    pixelData.t50 = -9999.;
    pixelData.fallTime = -9999.;
    pixelData.amplitude = -9999.;
    pixelData.electrons = -9999.;
    pixelData.RMS = -9999.;

    TString grName = Form("grEv%dPx%dsamp%d", event, pixel, fSamplingPeriodDictionary[pixel]);
    auto gr = (TGraph*)inFile->Get(grName.Data());

    if(!gr) 
    {
        std::cerr << std::endl << "Error: graph " << grName << " not found." << std::endl;
        return 0;
    }

    const int minimumIndex = FindMinimumIndex(*gr);
    const double minimum = gr->GetPointY(minimumIndex);

    auto edgesIdx = FindEdgeIndex(*gr, fIgnorePoints, fChecks, fNSample, fNDerivativePoints, fNSmoothingPoints);
    const double edgeLeft = gr->GetPointX(edgesIdx.first);
    const double edgeRight = gr->GetPointX(edgesIdx.second);

    // auto edges = FindEdge(*gr, fIgnorePoints, fChecks, fNSample, fNDerivativePoints, fNSmoothingPoints);
    // if edgeLeft is too little, it means that the signal is not present in the graph
    // therefore, in order to have a meaningful value for the baseline, we take the maximum of the graph as the edgeLeft
    // an event like this will have fallTime = -9999. (etc.)
    // const bool isSignalPresent = (edges.first > 1e-5);
    // const double edgeLeft = isSignalPresent ? gr->GetPointX(FindMaximumIndex(*gr)) : edges.first;
    // const double edgeRight = edges.second;
    

    int oldVerbosity = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;                         // disable implicit multithreading

    TF1 baselineFit(Form("baselineFit%d", pixel), "[0]");
    gr->Fit(&baselineFit, "Q+", "", gr->GetPointX(fIgnorePoints), edgeLeft);
    pixelData.baseline = baselineFit.GetParameter(0);

    TF1 minLevelFit(Form("minLevelFit%d", pixel), "[0]");
    gr->Fit(&minLevelFit, "Q+", "", edgeRight, gr->GetPointX(gr->GetN()-1-fIgnorePoints));
    pixelData.minLevel = minLevelFit.GetParameter(0);

    if (pixelData.minLevel  > pixelData.baseline)   pixelData.amplitude = -9999.;
    else                                            pixelData.amplitude = abs(pixelData.minLevel - pixelData.baseline);
    if (fmVToElectrons) 
    {
      if (pixelData.amplitude < 0.) pixelData.electrons = -9999.;
      else                          pixelData.electrons = (pixelData.amplitude - fmVToElectrons[pixel][0]) / fmVToElectrons[pixel][1];
    }

    if (edgeLeft > edgeRight || pixelData.minLevel > pixelData.baseline) 
    {
        pixelData.amplitude = -9999.;
        pixelData.t10 = -9999.;
        pixelData.t50 = -9999.;
        pixelData.t90 = -9999.;
        pixelData.fallTime = -9999.;
    }
    else
    {
        pixelData.t10 = gr->GetPointX(FindClosestPointIndex(*gr, pixelData.baseline - 0.1 * pixelData.amplitude, fIgnorePoints));
        pixelData.t50 = gr->GetPointX(FindClosestPointIndex(*gr, pixelData.baseline - 0.5 * pixelData.amplitude, fIgnorePoints));
        pixelData.t90 = gr->GetPointX(FindClosestPointIndex(*gr, pixelData.baseline - 0.9 * pixelData.amplitude, fIgnorePoints));
        pixelData.fallTime = pixelData.t50 - pixelData.t10;

        //TF1 fit(Form("fit%d", pixel), "pol1", edgeLeft, edgeRight);
        //gr->Fit(&fit, "RMQ+");
        //pixelData.t10 = fit.GetX(pixelData.baseline - 0.1 * pixelData.amplitude);
        //pixelData.t50 = fit.GetX(pixelData.baseline - 0.5 * pixelData.amplitude);
        //pixelData.t90 = fit.GetX(pixelData.baseline - 0.9 * pixelData.amplitude);
        //pixelData.fallTime = pixelData.t50 - pixelData.t10;
    }

    gErrorIgnoreLevel = oldVerbosity;                   // re-enable implicit multithreading

    auto grMeanAndRMS = GetMeanAndRMS(*gr, fIgnorePoints, fIgnorePoints + fNSample);
    pixelData.RMS = grMeanAndRMS.second;

    // Debugging session
    if (pixel == 9 and event == 5)
    {
        TFile outFile("ITS3/Data/Event5Pixel9.root", "recreate");
        gr->Write();
        TGraph * grPoints = new TGraph(0);
        grPoints->SetPoint(0, pixelData.t10, pixelData.baseline - 0.1 * pixelData.amplitude);
        grPoints->SetPoint(1, pixelData.t50, pixelData.baseline - 0.5 * pixelData.amplitude);
        grPoints->SetPoint(2, pixelData.t90, pixelData.baseline - 0.9 * pixelData.amplitude);
        grPoints->Write();

        auto canvas = new TCanvas("c", "", 1000, 1000);
        gr->Draw("AL");
        grPoints->SetMarkerSize(2);
        grPoints->SetMarkerStyle(20);
        grPoints->SetMarkerColor(kRed);
        grPoints->Draw("same");
        outFile.cd();
        canvas->Write();

        outFile.Close();

        std::cout << "Event 5, Pixel 9" << std::endl;
        std::cout << "Amplitude: " << pixelData.amplitude << std::endl;
        std::cout << "t10: " << pixelData.t10 << std::endl;
        std::cout << "t50: " << pixelData.t50 << std::endl;
        std::cout << "t90: " << pixelData.t90 << std::endl; 
        std::cout << "FallTime: " << pixelData.fallTime << std::endl;
        std::cout << std::endl;
    }
    /*
    if (pixel == 1)
    {
        std::cout << std::endl;
        std::cout << GREEN << "Event: " << event << RESET << std::endl;
        std::cout << "pixel: " << pixel << std::endl;
        std::cout << "MinimumPos: " << gr->GetPointX(minimumIndex) << std::endl;
        std::cout << "Minimum: " << minimum << std::endl;
        std::cout << "EdgeLeft: " << edgeLeft << std::endl;
        std::cout << "EdgeRight: " << edgeRight << std::endl;
        std::cout << "Baseline: " << pixelData.baseline << std::endl;
        std::cout << "MinLevel: " << pixelData.minLevel << std::endl;
        std::cout << "Amplitude: " << pixelData.amplitude << std::endl;
        std::cout << "t10: " << pixelData.t10 << std::endl;
        std::cout << "t50: " << pixelData.t50 << std::endl;
        std::cout << "t90: " << pixelData.t90 << std::endl;
        std::cout << "FallTime: " << pixelData.fallTime << std::endl;
        std::cout << "RMS: " << pixelData.RMS << std::endl;
    }
    */

    delete gr;
    return 1;
}

/**
 * @brief Function to process a single event read by the oscilloscope. 
 * It reads the graph from the input file and computes the variables of interest.
 * If no signal is found (i.e. the signal amplitude is below the set threshold) the time variables are set to -9999.
 * If a graph is not found, the function returns false.
 * 
 * @param event 
 * @param pixel 
 * @param pixelData 
 * @return true 
 * @return false 
 */
bool Preprocessor::ProcessEventADC(const int event, const int pixel, PixelData & pixelData, TFile * inFile)
{
    pixelData.pixel = pixel;  
    pixelData.samplingPeriod = fSamplingPeriodDictionary[pixel];
    pixelData.baseline = -9999.;
    pixelData.minLevel = -9999.;
    pixelData.t10 = -9999.;
    pixelData.t90 = -9999.;
    pixelData.t50 = -9999.;
    pixelData.fallTime = -9999.;
    pixelData.amplitude = -9999.;
    pixelData.electrons = -9999.;
    pixelData.RMS = -9999.;

    TString grName = Form("grEv%dPx%dsamp%d", event, pixel, fSamplingPeriodDictionary[pixel]);
    auto gr = (TGraph*)inFile->Get(grName.Data());

    if(!gr) 
    {
        std::cerr << std::endl << "Error: graph " << grName << " not found." << std::endl;
        return 0;
    }

    const int minimumIndex = FindMinimumIndex(*gr);
    const double minimum = gr->GetPointY(minimumIndex);
    pixelData.minLevel = minimum;
    auto edges = FindEdge(*gr);
    const double edgeLeft = gr->GetPointX(99);            // take the 100th point as the edgeLeft

    int oldVerbosity = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;                             // disable implicit multithreading

    TF1 baselineFit(Form("baselineFit%d", pixel), "[0]");
    gr->Fit(&baselineFit, "Q+", "", gr->GetPointX(0), gr->GetPointX(0));
    pixelData.baseline = baselineFit.GetParameter(0);

    gErrorIgnoreLevel = oldVerbosity;                       // re-enable implicit multithreading


    if (pixelData.baseline - pixelData.minLevel < fThreshold)   pixelData.amplitude = -9999.;
    else                                                        pixelData.amplitude = abs(pixelData.minLevel - pixelData.baseline);
    if (fmVToElectrons) pixelData.electrons = (pixelData.amplitude - fmVToElectrons[pixel][0]) / fmVToElectrons[pixel][1];

    // for the ADC data, time variables cannot be computed. They are set to -9999.
    pixelData.t10 = -9999.;
    pixelData.t50 = -9999.;
    pixelData.t90 = -9999.;
    pixelData.fallTime = -9999.;

    auto grMeanAndRMS = GetMeanAndRMS(*gr, fIgnorePoints, fIgnorePoints + fNSample);
    pixelData.RMS = grMeanAndRMS.second;

    // Debugging session
    /*
    if (pixel == 1)
    {
        std::cout << std::endl;
        std::cout << GREEN << "Event: " << event << RESET << std::endl;
        std::cout << "pixel: " << pixel << std::endl;
        std::cout << "MinimumPos: " << gr->GetPointX(minimumIndex) << std::endl;
        std::cout << "Minimum: " << minimum << std::endl;
        std::cout << "EdgeLeft: " << edgeLeft << std::endl;
        std::cout << "EdgeRight: " << edgeRight << std::endl;
        std::cout << "Baseline: " << pixelData.baseline << std::endl;
        std::cout << "MinLevel: " << pixelData.minLevel << std::endl;
        std::cout << "Amplitude: " << pixelData.amplitude << std::endl;
        std::cout << "t10: " << pixelData.t10 << std::endl;
        std::cout << "t50: " << pixelData.t50 << std::endl;
        std::cout << "t90: " << pixelData.t90 << std::endl;
        std::cout << "FallTime: " << pixelData.fallTime << std::endl;
        std::cout << "RMS: " << pixelData.RMS << std::endl;
    }
    */

    delete gr;
    return 1;
}

/**
 * @brief For a single event and pixel, draws the waveform and its derivative.
 * 
 * @param event 
 * @param pixel 
 * @param outFilePath 
 */
void Preprocessor::DrawPixel(const int event, const int pixel, const char * outFilePath)
{
    TString sOutFilePath(outFilePath);
    if (sOutFilePath.EqualTo("default"))    sOutFilePath = Form("ITS3/Data/Event%dPixel%d.root", event, pixel);
    
    if (pixel > fNPixels || event > fNEvents)
    {
        std::cerr << RED << "Error: pixel or event out of range." << RESET << std::endl;
        return;
    }

    auto inFile = TFile::Open(fInFilePath.Data());
    TString grName = Form("grEv%dPx%dsamp%d", event, pixel, fSamplingPeriodDictionary[pixel]);
    auto gr = (TGraph*)inFile->Get(grName.Data());
    auto grDerivative = DerivativeGraph(*gr);

    auto c = new TCanvas("c", "c", 1000, 1000);
    c->Divide(2);
    c->cd(1);
    gr->Draw("AL");
    c->cd(2);
    grDerivative->Draw("AL");
    //c->SaveAs(sOutFilePath.Data());
    c->Draw();

    auto outFile = TFile::Open(sOutFilePath.Data(), "RECREATE");
    gr->Write();
    grDerivative->Write();
    outFile->Close();
}

/**
 * @brief For a single event, draws the waveform and its derivative for all the pixels.
 * 
 * @param event 
 * @param outFilePath 
 * @param preprocessed 
 */
void Preprocessor::DrawEvent(const int event, const char * outFilePath)
{
    TString sOutFilePath(outFilePath);
    if (sOutFilePath.EqualTo("default"))    sOutFilePath = Form("ITS3/Data/Event%d.root", event);
    
    if (event > fNEvents)
    {
        std::cerr << RED << "Error: event out of range." << RESET << std::endl;
        return;
    }

    auto inFile = TFile::Open(fInFilePath.Data());
    auto outFile = TFile::Open(sOutFilePath.Data(), "RECREATE");
    std::cout << "Writing output file: " << BLUE << UNDERLINE << sOutFilePath << RESET << std::endl;

    for (int i = 0; i < fNPixels; i++)
    {
        TString grName = Form("grEv%dPx%dsamp%d", event, i, fSamplingPeriodDictionary[i]);
        
        auto canvas = new TCanvas(Form("cEv%dPx%dsamp%d", event, i, fSamplingPeriodDictionary[i]), "", 1000, 1000);
        auto canvasClean = new TCanvas(Form("cEv%dPx%dsamp%d_clean", event, i, fSamplingPeriodDictionary[i]), "", 1000, 1000);

        auto gr = (TGraph*)inFile->Get(grName.Data());
        auto grDerivative = DerivativeGraph(*gr);
        grDerivative->SetTitle(Form("Signal derivative - Event %d, Pixel %d, Sampling Period %d ps; Time (ns); Derivative (a.u.)", event, i, fSamplingPeriodDictionary[i]));
        auto grSecondDerivative = DerivativeGraph(*grDerivative);
        auto grSmooth = SmoothGraph(*gr, 10);
        canvas->Divide(2);

        auto edgesIdx = FindEdgeIndex(*gr, fIgnorePoints, fChecks, fNSample, fNDerivativePoints, fNSmoothingPoints);
        auto grEdge = new TGraph(2);
        grEdge->SetPoint(0, gr->GetPointX(edgesIdx.first), gr->GetPointY(edgesIdx.first));
        grEdge->SetPoint(1, gr->GetPointX(edgesIdx.second), gr->GetPointY(edgesIdx.second));
        grEdge->SetMarkerStyle(20);
        grEdge->SetMarkerColor(kOrange-3);

        auto baselineFit = new TF1(Form("baselineFit%d", i), "[0]");
        baselineFit->SetLineColor(kRed);
        gr->Fit(baselineFit, "Q+", "", gr->GetPointX(fIgnorePoints), gr->GetPointX(edgesIdx.first));
        auto minLevelFit = new TF1(Form("minLevelFit%d", i), "[0]");
        minLevelFit->SetLineColor(kRed);
        gr->Fit(minLevelFit, "Q+", "", gr->GetPointX(edgesIdx.second), gr->GetPointX(gr->GetN()-1-fIgnorePoints));

        gStyle->SetOptFit(0);

        // recreate gr without the first and last nIgnorePoints points
        auto grClean = new TGraph(gr->GetN() - 2 * fIgnorePoints);
        for (int iPoint = fIgnorePoints; iPoint < gr->GetN() - fIgnorePoints; ++iPoint) grClean->SetPoint(iPoint - fIgnorePoints, gr->GetPointX(iPoint), gr->GetPointY(iPoint));
        grClean->SetName(Form("grEv%dPx%dsamp%d_clean", event, i, fSamplingPeriodDictionary[i]));
        grClean->SetTitle(Form("Event %d, Pixel %d, Sampling Period %d ps; Time (ns); Amplitude (mV)", event, i, fSamplingPeriodDictionary[i]));

        gr->Write();
        grClean->Write();
        grDerivative->Write();
        grSecondDerivative->Write();
        grSmooth->Write();

        canvas->cd(1);
        grClean->Draw("AL");
        grEdge->Draw("same");
        baselineFit->Draw("same");
        minLevelFit->Draw("same");
        canvas->cd(2);
        grDerivative->Draw("same");
        canvas->Write();

        canvasClean->cd();
        grClean->Draw("AL");
        baselineFit->SetRange(grClean->GetPointX(0), grClean->GetPointX(edgesIdx.first-fIgnorePoints));
        minLevelFit->SetRange(grClean->GetPointX(edgesIdx.second+fIgnorePoints), grClean->GetPointX(grClean->GetN()-1));
        baselineFit->Draw("same");
        minLevelFit->Draw("same");
        grEdge->Draw("same");
        canvasClean->Write();

        delete gr;
        delete grClean;
        delete grDerivative;
        delete canvas;
        delete baselineFit;
        delete minLevelFit;
        delete grEdge;
        delete canvasClean;
    }

    outFile->Close();
    inFile->Close();    
}

// -------------------------------------------------------------------------------------
// functions used in this file

/**
 * @brief Function to find the edge of the signal.
 * 
 * @param gr 
 * @param nIgnorePoints number of points to ignore at the beginning and at the end of the TGraph (due to the smoothing and/or 
 * peculiarity of the waveform)
 * @return std::tuple<double, double, double> edgeLeft, edgeRight
 */
std::pair<double, double> FindEdge(TGraph & gr, int nIgnorePoints, int nChecks, int nSample, const int nDerivativePoints, const int nSmoothingPoints)
{
    auto edgesIdx = FindEdgeIndex(gr, nIgnorePoints, nChecks, nSample, nDerivativePoints, nSmoothingPoints);
    return std::make_pair(gr.GetPointX(edgesIdx.first), gr.GetPointX(edgesIdx.second));
}

/**
 * @brief Function to find the edge of the signal.
 * 
 * @param gr 
 * @return std::tuple<double, double, double> edgeLeft, edgeRight
 */
std::pair<int, int> FindEdgeIndex(TGraph & gr, int nIgnorePoints, int nChecks, int nSample, const int nDerivativePoints, const int nSmoothingPoints)
{
    auto grSmooth = SmoothGraph(gr, nSmoothingPoints);                  // smooth the graph to find the edge more easily
    auto grDerivative = DerivativeGraph(*grSmooth, nDerivativePoints);  // take the derivative of the smoothed graph

    double meanDer{0.}, RMSDer{0.};
    if (nSample > gr.GetN())  nSample = int(gr.GetN()/2);
    if (nIgnorePoints == 0)  nIgnorePoints = nSmoothingPoints;
    auto resultDer = GetMeanAndRMS(*grDerivative, nIgnorePoints, nIgnorePoints+nSample);
    meanDer = resultDer.first;
    RMSDer = resultDer.second;

    double edgeLeft{0.}, edgeRight{0.};    
    int edgeLeftIndex = FindChangingDerivative(*grDerivative, meanDer, RMSDer, nIgnorePoints, nChecks);          // correct for smoothing
    int edgeRightIndex = FindChangingDerivative(*grDerivative, meanDer, RMSDer, nIgnorePoints, nChecks, true);   // correct for smoothing

    delete grSmooth;
    delete grDerivative;
    
    return std::make_pair(edgeLeftIndex, edgeRightIndex);
}

/**
 * @brief Function to find the minimum value in a graph.
 * 
 * @param gr 
 * @return double 
 */
int FindMinimumIndex(TGraph & gr)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    double minimum{999999.};
    int minimumIndex{0};
    for (int i = 0; i < nSamples; ++i)
    {
        if (y[i] < minimum)     
        {
            minimum = y[i];
            minimumIndex = i;
        }   
    }

    return minimumIndex;
}

/**
 * @brief Function to find the maximum value in a graph.
 * 
 * @param gr 
 * @return double 
 */
int FindMaximumIndex(TGraph & gr)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    double maximum{-9999999.};
    int maximumIndex{0};
    for (int i = 0; i < nSamples; ++i)
    {
        if (y[i] > maximum)     
        {
            maximum = y[i];
            maximumIndex = i;
        }   
    }

    return maximumIndex;
}

/**
 * @brief Find the point at which the derivative changes for n consecutive points with a value higher or lower than meanDer ± 3 * RMSDer.
 * 
 * @param grDer 
 * @param meanDer 
 * @param RMSDer 
 * @param ignorePoints Number of points to ignore at the beginning and at the end of the TGraph (due to smoothing).
 * @param nChecks Number of consecutive points with derivative beyond meanDer ± 3 * RMSDer to consider a change.
 * @param backwards If true, the function loops through the TGraph from the last point to the first one.
 * @return const int 
 */
int FindChangingDerivative(TGraph & grDer, const double & meanDer, const double & RMSDer, const int ignorePoints, const int nChecks, const bool backwards)
{
    const int nSamples = grDer.GetN();
    double * x = grDer.GetX();
    double * y = grDer.GetY();

    int changingPoint{0};
    int consecutiveCount{0};

    if (!backwards)
    {
        for (int i = ignorePoints; i < nSamples - nChecks + 1; ++i)
        {
            // Check if the derivative of n consecutive points is beyond meanDer ± 3 * RMSDer
            bool isChange = true;
            for (int j = 0; j < nChecks; ++j)
            {
                if (y[i + j] <= meanDer + 5 * RMSDer && y[i + j] >= meanDer - 5 * RMSDer)
                {
                    isChange = false;
                    break;
                }
            }

            if (isChange)
            {
                changingPoint = i;
                consecutiveCount = nChecks;
                break;
            }
        }
    }
    else
    {
        for (int i = nSamples - ignorePoints - 1; i > nChecks - 2; --i)
        {
            // Check if the derivative of n consecutive points is beyond meanDer ± 3 * RMSDer
            bool isChange = true;
            for (int j = 0; j < nChecks; ++j)
            {
                if (y[i - j] <= meanDer + 5 * RMSDer && y[i - j] >= meanDer - 5 * RMSDer)
                {
                    isChange = false;
                    break;
                }
            }

            if (isChange)
            {
                changingPoint = i;
                consecutiveCount = nChecks;
                break;
            }
        }
    }

    if (consecutiveCount < nChecks)
    {
        // If no change is found, return 0
        changingPoint = 0;
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
std::pair<double, double> GetMeanAndRMS(TGraph & gr, const int & begin, const int & end)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    double mean{0.};
    for (int i = begin; i < end; ++i)
    {
        mean += y[i];
    }
    mean /= (end - begin);

    double RMS{0.};
    for (int i = begin; i < end; ++i)
    {
        RMS += (y[i] - mean) * (y[i] - mean);
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
TGraph* SmoothGraph(TGraph& gr, const int nSmoothingPoints)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    auto grSmooth = new TGraph(0);
    grSmooth->SetName(Form("%s_smoothed", gr.GetName()));
    for (int i = 0; i < nSamples; ++i)
    {
        double sum{0.};
        for (int j = i - nSmoothingPoints; j < i + nSmoothingPoints; ++j)
        {
            if (j < 0 || j > nSamples) continue;
            sum += y[j];
        }
        grSmooth->SetPoint(i, x[i], sum / (2 * nSmoothingPoints + 1));
    }

    return grSmooth;
}

/**
 * @brief Function to obtain the derivative of a TGraph point by point. The derivative is computed as the average of the derivative done on 
 * nDerivativePoints points before and after the point in which the derivative is evaluated.
 * 
 * @param gr 
 * @param nDerivativePoints
 * @return TGraph* 
 */
TGraph* DerivativeGraph(TGraph& gr, const int nDerivativePoints)
{
    // compute the derivative of a TGraph point by point considering nDerivativePoints points before and after the point in which the derivative is evaluated
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    auto grDerivative = new TGraph(0);
    grDerivative->SetName(Form("%s_derivative", gr.GetName()));
    
    for (int i = 1; i < nSamples - nDerivativePoints - 1; ++i)
    {
        double der{0.};
        for (int j = i; j < i + nDerivativePoints; j++)   der += (y[j + 1] - y[j - 1]) / (x[j + 1] - x[j - 1]);
        der /= nDerivativePoints;
        grDerivative->SetPoint(i, x[i], der);
    }
    grDerivative->SetPoint(0, x[0], grDerivative->GetPointY(1));
    grDerivative->SetPoint(nSamples - 1, x[nSamples - 1], grDerivative->GetPointY(nSamples - 2));

    // set the first and last point of the derivative graph to the second and second to last point of the original graph
    grDerivative->SetPoint(0, x[0], grDerivative->GetPointY(1));
    grDerivative->SetPoint(nSamples - 1, x[nSamples - 1], grDerivative->GetPointY(nSamples - 2));

    return grDerivative;
}

/**
 * @brief Find the first point in the TGraph with given value. The first skipPoints points are ignored.
 * 
 * @param gr 
 * @param x 
 * @param skipPoints 
 * @return int 
 */
int FindPointIndex(TGraph & gr, const double & x, const int nIgnorePoints)
{
    const int nSamples = gr.GetN();
    double * xGraph = gr.GetX();
    double * yGraph = gr.GetY();

    for (int i = nIgnorePoints; i < nSamples; ++i)  if (xGraph[i] == x) return i;
    return 0;
}

/**
 * @brief Find the point in TGraph with value closest to the target value.
 * 
 * @param gr
 * @param target
 * @return int
*/
int FindClosestPointIndex(TGraph & gr, const double target, const int nIgnorePoints)
{
    const int size = gr.GetN();
    double * array = gr.GetY();

    size_t effectiveSize = (size > nIgnorePoints) ? size - nIgnorePoints : 0;
    array += nIgnorePoints;

    // Create a vector of indices and sort it based on the values in the array
    size_t * indices = new size_t[effectiveSize];
    std::iota(indices, indices + effectiveSize, 0);
    std::sort(indices, indices + effectiveSize,
              [&array](size_t i1, size_t i2) { return array[i1] < array[i2]; });

    // Use binary search on the sorted indices
    auto it = std::lower_bound(indices, indices + effectiveSize, target,
                               [&array](size_t i, double value) { return array[i] < value; });

    if (it == indices) {
        // Target is less than or equal to the first element
        size_t result = *it;
        delete[] indices;
        return result;
    } else if (it == indices + effectiveSize) {
        // Target is greater than or equal to the last element
        size_t result = *(--it);
        delete[] indices;
        return result;
    }

    // Check the closest value between the current iterator and the one before
    size_t index1 = *it;
    size_t index2 = *(--it);

    size_t result = std::abs(array[index1] - target) < std::abs(array[index2] - target) ? index1 : index2;
    result += nIgnorePoints;

    delete[] indices;
    return result;
}

// -------------------------------------------------------------------------------------
// progress bar section
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
    std::cout << "] " << static_cast<int>(percentage * 100.0) << "%";
    std::cout << "  Elapsed Time: " << elapsedTime << "s";
    std::cout.flush();
}
// -------------------------------------------------------------------------------------
