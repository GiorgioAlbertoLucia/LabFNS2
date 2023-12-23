// Implementation of Preprocessor class

#include "preprocessor.hh"
#include "colorTerminal.hh"

#include <iostream>
#include <fstream>
#include <streambuf>

#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TKey.h>
#include <TTree.h>
#include <TList.h>
#include <TCanvas.h>


// -------------------------------------------------------------------------------------
// functions used in this file

std::pair<double, double> FindEdge(TGraph & gr);
int FindMinimumIndex(TGraph & gr);
int FindMaximumIndex(TGraph & gr);
int FindChangingDerivative(TGraph & grDer, const double& meanDer, const double& RMSDer, const int ignorePoints, const bool backwards = false);
std::pair<double, double> GetMeanAndRMS(TGraph & gr, const int& begin, const int& end);
TGraph* SmoothGraph(TGraph& gr, const int nSmoothingPoints = 2);
TGraph* DerivativeGraph(TGraph& gr);

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
    fThreshold(threshold)
{
    Preprocessor::ReadInput();   
    Preprocessor::GenerateSamplingDictionary();
}

Preprocessor::~Preprocessor()
{
    delete fSamplingPeriodDictionary;
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
 * If no signal is found (i.e. the signal amplitude is below the set threshold) the time variables are set to -1.
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
    pixelData.baseline = -1.;
    pixelData.minLevel = -1.;
    pixelData.t10 = -1.;
    pixelData.t90 = -1.;
    pixelData.t50 = -1.;
    pixelData.fallTime = -1.;
    pixelData.amplitude = -1.;
    pixelData.RMS = -1.;

    TString grName = Form("grEv%dPx%dsamp%d", event, pixel, fSamplingPeriodDictionary[pixel]);
    auto gr = (TGraph*)inFile->Get(grName.Data());

    if(!gr) 
    {
        std::cerr << std::endl << "Error: graph " << grName << " not found." << std::endl;
        return 0;
    }

    const int minimumIndex = FindMinimumIndex(*gr);
    const double minimum = gr->GetPointY(minimumIndex);
    auto edges = FindEdge(*gr);
    const double edgeLeft = (edges.first < 1e-5) ? gr->GetPointX(FindMaximumIndex(*gr)) : edges.first;
    const double edgeRight = edges.second;
    // if edgeLeft is too little, it means that the signal is not present in the graph
    // therefore, in order to have a meaningful value for the baseline, we take the maximum of the graph as the edgeLeft

    // avoid automatic output due to TF1s not being able to fit the graph
    std::streambuf* coutBuffer = std::cout.rdbuf();
    std::streambuf* cerrBuffer = std::cerr.rdbuf();
    std::ofstream nullStream("/dev/null"); 
    std::cout.rdbuf(nullStream.rdbuf());
    std::cerr.rdbuf(nullStream.rdbuf());

    TF1 baselineFit(Form("baselineFit%d", pixel), "[0]");
    gr->Fit(&baselineFit, "Q", "", gr->GetPointX(0), edgeLeft);
    pixelData.baseline = baselineFit.GetParameter(0);

    TF1 minLevelFit(Form("minLevelFit%d", pixel), "[0]");
    gr->Fit(&minLevelFit, "Q", "", edgeRight, gr->GetPointX(minimumIndex));
    pixelData.minLevel = minLevelFit.GetParameter(0);

    std::cout.rdbuf(coutBuffer);        // restore the original state of std::cout
    std::cerr.rdbuf(cerrBuffer);        // restore the original state of std::cerr

    pixelData.amplitude = abs(pixelData.minLevel - pixelData.baseline);

    if (pixelData.amplitude < fThreshold) 
    {
        pixelData.t10 = -1.;
        pixelData.t50 = -1.;
        pixelData.t90 = -1.;
        pixelData.fallTime = -1.;
    }
    else
    {
        TF1 fit(Form("fit%d", pixel), "pol1", edgeLeft, edgeRight);
        gr->Fit(&fit, "RMLQ+");
        pixelData.t10 = fit.GetX(pixelData.baseline - 0.1 * pixelData.amplitude);
        pixelData.t50 = fit.GetX(pixelData.baseline - 0.5 * pixelData.amplitude);
        pixelData.t90 = fit.GetX(pixelData.baseline - 0.9 * pixelData.amplitude);
        pixelData.fallTime = pixelData.t50 - pixelData.t10;
    }

    const int nSamples = 80;
    auto grMeanAndRMS = GetMeanAndRMS(*gr, 0, nSamples);
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
 * @brief Function to process a single event read by the oscilloscope. 
 * It reads the graph from the input file and computes the variables of interest.
 * If no signal is found (i.e. the signal amplitude is below the set threshold) the time variables are set to -1.
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
    pixelData.baseline = -1.;
    pixelData.minLevel = -1.;
    pixelData.t10 = -2.;
    pixelData.t90 = -2.;
    pixelData.t50 = -2.;
    pixelData.fallTime = -2.;
    pixelData.amplitude = -1.;
    pixelData.RMS = -1.;

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
    //const double edgeLeft = gr->GetPointX(99);              // take the 100th point as the edgeLeft
    const double edgeLeft = 250*99;

    // avoid automatic output due to TF1s not being able to fit the graph
    std::streambuf* coutBuffer = std::cout.rdbuf();
    std::streambuf* cerrBuffer = std::cerr.rdbuf();
    std::ofstream nullStream("/dev/null"); 
    std::cout.rdbuf(nullStream.rdbuf());
    std::cerr.rdbuf(nullStream.rdbuf());

    TF1 baselineFit(Form("baselineFit%d", pixel), "[0]");
    gr->Fit(&baselineFit, "Q", "", gr->GetPointX(0), edgeLeft);
    pixelData.baseline = baselineFit.GetParameter(0);

    std::cout.rdbuf(coutBuffer);        // restore the original state of std::cout
    std::cerr.rdbuf(cerrBuffer);        // restore the original state of std::cerr

    pixelData.amplitude = abs(gr->GetPointY(FindMinimumIndex(*gr)) - pixelData.baseline);

    // for the ADC data, time variables cannot be computed. They are set to -2.
    pixelData.t10 = -2.;
    pixelData.t50 = -2.;
    pixelData.t90 = -2.;
    pixelData.fallTime = -2.;

    const int nSamples = 80;
    auto grMeanAndRMS = GetMeanAndRMS(*gr, 0, nSamples);
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
 * @brief For a single event, draws the waveform and its derivative.
 * 
 * @param event 
 * @param pixel 
 * @param outFilePath 
 */
void Preprocessor::DrawEvent(const int event, const int pixel, const char * outFilePath)
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

// -------------------------------------------------------------------------------------
// functions used in this file

/**
 * @brief Function to find the edge of the signal.
 * 
 * @param gr 
 * @return std::tuple<double, double, double> edgeLeft, edgeRight
 */
std::pair<double, double> FindEdge(TGraph & gr)
{
    const int nSmoothingPoints = 5;                         // number of points to consider for the smoothing
    auto grSmooth = SmoothGraph(gr, nSmoothingPoints);      // smooth the graph to find the edge more easily
    auto grDerivative = DerivativeGraph(*grSmooth);         // take the derivative of the smoothed graph

    const int nSamples = 80;                                // number of samples to consider for the plateau
    double meanDer{0.}, RMSDer{0.};
    auto resultDer = GetMeanAndRMS(*grDerivative, nSmoothingPoints, nSmoothingPoints+nSamples);
    meanDer = resultDer.first;
    RMSDer = resultDer.second;

    double edgeLeft{0.}, edgeRight{0.};
    int edgeLeftIndex = FindChangingDerivative(*grDerivative, meanDer, RMSDer, nSmoothingPoints);          // correct for smoothing
    int edgeRightIndex = FindChangingDerivative(*grDerivative, meanDer, RMSDer, nSmoothingPoints, true);   // correct for smoothing
    edgeLeft = gr.GetPointX(edgeLeftIndex);
    edgeRight = gr.GetPointX(edgeRightIndex);

    delete grSmooth;
    delete grDerivative;
    
    return std::make_pair(edgeLeft, edgeRight);
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
 * @brief Find the pont in which the derivative changes for three consecutive points with a value higher than meanDer + 3 * RMSDer.
 * 
 * @param grDer 
 * @param meanDer 
 * @param RMSDer 
 * @param ignorePoints number of points to ignore at the beginning and at the end of the TGraph (due to the smoothing)
 * @param backwards if true, the function loops through the TGraph from the last point to the first one
 * @return const int 
 */
int FindChangingDerivative(TGraph & grDer, const double& meanDer, const double& RMSDer, const int ignorePoints, const bool backwards)
{
    const int nSamples = grDer.GetN();
    double * x = grDer.GetX();
    double * y = grDer.GetY();

    int changingPoint{0};
    if (!backwards)
    {
        for (int i = ignorePoints; i < nSamples - 2; ++i)
        {
            // positive derivative
            if (y[i] > meanDer + 3 * RMSDer && y[i + 1] > meanDer + 3 * RMSDer && y[i + 2] > meanDer + 3 * RMSDer)
            {
                changingPoint = i;
                break;
            }
            // negative derivative
            if (y[i] < meanDer - 3 * RMSDer && y[i + 1] < meanDer - 3 * RMSDer && y[i + 2] < meanDer - 3 * RMSDer)
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
            if (y[i] > meanDer + 3 * RMSDer && y[i - 1] > meanDer + 3 * RMSDer && y[i - 2] > meanDer + 3 * RMSDer)
            {
                changingPoint = i;
                break;
            }
            // negative derivative  
            if (y[i] < meanDer - 3 * RMSDer && y[i - 1] < meanDer - 3 * RMSDer && y[i - 2] < meanDer - 3 * RMSDer)
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
std::pair<double, double> GetMeanAndRMS(TGraph & gr, const int& begin, const int& end)
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
 * @brief Function to obtain the derivative of a TGraph point by point.
 * 
 * @param gr 
 * @return TGraph* 
 */
TGraph* DerivativeGraph(TGraph& gr)
{
    // compute the derivative of a TGraph point by point considering nDerivativePoints points before and after the point in which the derivative is evaluated
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    auto grDerivative = new TGraph(0);
    grDerivative->SetName(Form("%s_derivative", gr.GetName()));
    for (int i = 0; i < nSamples; ++i)
    {
        double der{0.};
        for (int j = i - 1; j < i + 1; ++j)
        {
            if (j < 0 || j > nSamples) continue;
            der = (y[j + 1] - y[j - 1]) / (x[j + 1] - x[j - 1]);
        }
        grDerivative->SetPoint(i, x[i], der);
    }

    // set the first and last point of the derivative graph to the second and second to last point of the original graph
    grDerivative->SetPoint(0, x[0], grDerivative->GetPointY(1));
    grDerivative->SetPoint(nSamples - 1, x[nSamples - 1], grDerivative->GetPointY(nSamples - 2));

    return grDerivative;
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
    //std::cout << "] " << std::setw(3) << static_cast<int>(percentage * 100.0) << "%";
    std::cout << "  Elapsed Time: " << elapsedTime << "s";
    std::cout.flush();
}
// -------------------------------------------------------------------------------------
