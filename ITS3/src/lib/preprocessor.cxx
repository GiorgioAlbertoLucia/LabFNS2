// Implementation of Preprocessor class

#include "preprocessor.hh"
#include "graphUtilities.hh"
#include "progressBar.hh"
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
    fNSmoothingPoints(10),
    fClusterThreshold(10)
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
    /*
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
    */
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
 * @brief Function to process a single evet and to find the seed and the cluster.
 * 
 * @param clusterSize
 * @param seedData
 * @param clusterData
 * @param pixelData
 * 
 * @return true
*/
bool Preprocessor::GenerateSeedAndCluster(int & clusterSize, PixelData & seedData, PixelData & clusterData, PixelData * pixelData)
{
    clusterSize = 0;

    seedData.pixel = -fNPixels;  
    seedData.samplingPeriod = 0;
    seedData.baseline = 0.;
    seedData.minLevel = 0.;
    seedData.t10 = 0.;
    seedData.t90 = 0.;
    seedData.t50 = 0.;
    seedData.fallTime = -9999.;
    seedData.amplitude = -9999.;
    seedData.electrons = -9999.;
    seedData.RMS = 0.;

    clusterData.pixel = 0;  
    clusterData.samplingPeriod = 0;
    clusterData.baseline = 0.;
    clusterData.minLevel = 0.;
    clusterData.t10 = 0.;
    clusterData.t90 = 0.;
    clusterData.t50 = 0.;
    clusterData.fallTime = 0.;
    clusterData.amplitude = 0.;
    clusterData.electrons = 0.;
    clusterData.RMS = 0.;


    for (int ipixel = 0; ipixel < fNPixels; ipixel++)
    {
        std::cout << "ampl pixel " << ipixel << ": " << (pixelData+ipixel)->amplitude << std::endl;
        if ((pixelData+ipixel)->amplitude > seedData.amplitude)
        {
            seedData.pixel = ipixel;
            seedData.amplitude = (pixelData+ipixel)->amplitude;
            seedData.electrons = (pixelData+ipixel)->electrons;
            seedData.fallTime= (pixelData+ipixel)->fallTime;
        }
        if ((pixelData+ipixel)->amplitude > fClusterThreshold)
        {
            clusterSize++;
            clusterData.amplitude += (pixelData+ipixel)->amplitude;
            clusterData.electrons += (pixelData+ipixel)->electrons;
        } 
    }
    if (clusterSize == 0) 
    {
        clusterData.amplitude = -9999.;
        clusterData.electrons = -9999.;
    }

    return 1;
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
        grDerivative->SetTitle(Form("Signal derivative - Event %d, Pixel %d, Sampling Period %d ps; Time (ps); Derivative (a.u.)", event, i, fSamplingPeriodDictionary[i]));
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
        grClean->SetTitle(Form("Event %d, Pixel %d, Sampling Period %d ps; Time (ps); Amplitude (mV)", event, i, fSamplingPeriodDictionary[i]));

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

/**
 * @brief Function to generate the seed and cluster tree from tree generated in BuildTree().
*/
void Preprocessor::BuildSeedAndClusterTree(const char * inFilePath, const char * outFilePath)
{
    TString sInFilePath(inFilePath);
    if (sInFilePath.EqualTo("default")) 
    {
        sInFilePath = fInFilePath;
        sInFilePath.ReplaceAll(".root", "_preprocessed.root");
    }
    TString sOutFilePath(outFilePath);
    if (sOutFilePath.EqualTo("default")) 
    {
        sOutFilePath = fInFilePath;
        sOutFilePath.ReplaceAll(".root", "_seed_and_cluster.root");
    }

    auto inFile = TFile::Open(sInFilePath.Data());
    std::cout << "Reading input file: " << BLUE << UNDERLINE << sInFilePath.Data() << RESET << std::endl;
    auto inTree = (TTree*)inFile->Get("PreprocessedData");

    const int nPixels = fNPixels;
    PixelData * pixelData[nPixels];
    for (PixelData *& pxDt: pixelData)  pxDt = new PixelData();
    TBranch * bPixelData[nPixels];
    for (int ipixel = 0; ipixel < nPixels; ipixel++)    inTree->SetBranchAddress(Form("pixel%d", ipixel), pixelData+ipixel, &bPixelData[ipixel]);

    TFile outFile(sOutFilePath.Data(), "RECREATE");
    TTree outTree("SeedAndCluster", "SeedAndCluster");
    std::cout << "Writing output file: " << BLUE << UNDERLINE << sOutFilePath << RESET << std::endl;

    int event{0}, clusterSize{0};
    PixelData seedData, clusterData;

    outTree.Branch("Event", &event, "Event/I");
    outTree.Branch("ClusterSize", &clusterSize,"ClusterSize/I");
    outTree.Branch("Seed", &seedData);
    outTree.Branch("Cluster", &clusterData);

    const auto startTime = std::chrono::steady_clock::now();
    for (int ievent = 0; ievent < fNEvents; ++ievent)
    //for (int ievent = 0; ievent < 3; ++ievent)
    {
        updateProgressBar(ievent, fNEvents, startTime);
        
        event = ievent;
        inTree->GetEntry(ievent);
        for (int ipixel = 0; ipixel < nPixels; ipixel++)    std::cout << "ampl pixel " << ipixel << ": " << pixelData[ipixel]->amplitude << std::endl;
        Preprocessor::GenerateSeedAndCluster(clusterSize, seedData, clusterData, *pixelData);

        outTree.Fill();
    }

    inFile->Close();
    std::cout << std::endl;

    outFile.cd();
    outTree.Write();
    outFile.Close();
}
