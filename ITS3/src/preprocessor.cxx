// Implementation of Preprocessor class

#include "preprocessor.hh"
#include "colorTerminal.hh"

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
double FindMinimumIndex(TGraph & gr);
int FindChangingDerivative(TGraph & grDer, const double& meanDer, const double& RMSDer, const int ignorePoints, const bool backwards = false);
std::pair<double, double> GetMeanAndRMS(TGraph & gr, const int& begin, const int& end);
TGraph* SmoothGraph(TGraph& gr, const int nSmoothingPoints = 2);
TGraph* DerivativeGraph(TGraph& gr);

// -------------------------------------------------------------------------------------


// progress bar section
// -------------------------------------------------------------------------------------
#include <Riostream.h>
#include <chrono>
#include <thread>
void updateProgressBar(int progress, int total, const std::chrono::steady_clock::time_point& startTime);
// -------------------------------------------------------------------------------------


Preprocessor::Preprocessor(const char * inFilePath, const double threshold): 
    fInFilePath(inFilePath), 
    fNChannels(0), 
    fNEvents(0), 
    fSamplingPeriod(0), 
    fThreshold(threshold)
{
    Preprocessor::ReadInput();   
}

Preprocessor::~Preprocessor()
{

}

/**
 * @brief Loops throught the whole file and finds information such as number of channels and number of events.
 * 
 */
void Preprocessor::ReadInput()
{
    auto inFile = TFile::Open(fInFilePath.Data());
    long int nKeys = inFile->GetNkeys();
    auto listKeys = inFile->GetListOfKeys();
    int event{0}, channel{0};
    const auto startTime = std::chrono::steady_clock::now();

    std::cout << "Reading input file: " << BLUE << UNDERLINE << fInFilePath << RESET << std::endl;
    for (long int i = 0; i < nKeys; i++)             // get the information looping through the names of the objects in the file
    {   
        updateProgressBar(i, nKeys, startTime);

        
        TKey* key = (TKey*)listKeys->At(i);             // get the ith object
        if(!key) 
        {
            //std::cerr << "\tError: key " << i << " not found." << std::endl;
            continue;              
        }            
        TString className = key->GetClassName();        // get the type of the object
        TString objectName = key->GetName();            // get the name of the object

        if (className == "TGraph")  sscanf(objectName.Data(), "grEv%dChanC%dsamp%d", &event, &channel, &fSamplingPeriod);
        if (event > fNEvents) fNEvents = event;    
        if (channel > fNChannels) fNChannels = channel;

        delete key;
    }

    std::cout << std::endl;
    fNEvents++;
    fNChannels++;

    //delete listKeys;
    inFile->Close();
    
}

void Preprocessor::BuildTree(const char * outFilePath)
{
    TString sOutFilePath(outFilePath);
    if (sOutFilePath.EqualTo("default")) 
    {
        sOutFilePath = fInFilePath;
        sOutFilePath.ReplaceAll(".root", "_preprocessed.root");
    }

    const int nChannels = fNChannels;
    PixelData channelData[nChannels];

    TFile outFile(sOutFilePath.Data(), "RECREATE");
    TTree outTree("PreprocessedData", "outTree");

    int event{0};
    outTree.Branch("Event", &event, "Event/I");
    for(int i = 1; i < fNChannels + 1; i++) outTree.Branch(Form("Channel%d", i), &channelData[i-1]);//, PixelData::GetBranchList().Data());
    
    const auto startTime = std::chrono::steady_clock::now();

    for (int ievent = 0; ievent < fNEvents; ++ievent)
    {
        updateProgressBar(ievent, fNEvents, startTime);
        
        event = ievent;
        for (int ichannel = 1; ichannel < fNChannels + 1; ichannel++)   // channels start from 1
        {
            bool outcome = Preprocessor::ProcessEvent(ievent, ichannel, channelData[ichannel-1]);
            if (!outcome) continue;                                     // if the event is not found, skip it
        }
        
        outTree.Fill();
    }

    std::cout << std::endl;
    outFile.cd();
    outTree.Write();
    outFile.Close();
}

/**
 * @brief Function to process a single event. It reads the graph from the input file and computes the variables of interest.
 * If no signal is found (i.e. the signal amplitude is below the set threshold) the time variables are set to -1.
 * If a graph is not found, the function returns false.
 * 
 * @param event 
 * @param channel 
 * @param channelData 
 * @return true 
 * @return false 
 */
bool Preprocessor::ProcessEvent(const int event, const int channel, PixelData& channelData)
{
    channelData.channel = channel;  
    channelData.baseline = 0.;
    channelData.minLevel = 0.;
    channelData.t10 = 0.;
    channelData.t90 = 0.;
    channelData.t50 = 0.;
    channelData.fallTime = 0.;
    channelData.amplitude = 0.;
    channelData.RMS = 0.;

    auto inFile = TFile::Open(fInFilePath.Data());
    TString grName = Form("grEv%dChanC%dsamp%d", event, channel, fSamplingPeriod);
    auto gr = (TGraph*)inFile->Get(grName.Data());

    if(!gr) 
    {
        std::cerr << "Error: graph " << grName << " not found." << std::endl;
        return 0;
    }

    const int minimumIndex = FindMinimumIndex(*gr);
    const double minimum = gr->GetPointY(minimumIndex);
    auto edges = FindEdge(*gr);
    const double edgeLeft = edges.first;
    const double edgeRight = edges.second;

    // Debugging session

    TF1 baselineFit("baselineFit", "[0]");
    gr->Fit(&baselineFit, "Q", "", gr->GetPointX(0), edgeLeft);
    channelData.baseline = baselineFit.GetParameter(0);

    TF1 minLevelFit("minLevelFit", "[0]");
    gr->Fit(&minLevelFit, "Q", "", edgeRight, gr->GetPointX(minimumIndex));
    channelData.minLevel = minLevelFit.GetParameter(0);

    channelData.amplitude = abs(channelData.minLevel - channelData.baseline);

    if (channelData.amplitude < fThreshold) 
    {
        channelData.t10 = -1.;
        channelData.t50 = -1.;
        channelData.t90 = -1.;
        channelData.fallTime = -1.;
    }
    else
    {
        TF1 fit("fit", "pol1", edgeLeft, edgeRight);
        gr->Fit(&fit, "RMLQ+");
        channelData.t10 = fit.GetX(channelData.baseline - 0.1 * channelData.amplitude);
        channelData.t50 = fit.GetX(channelData.baseline - 0.5 * channelData.amplitude);
        channelData.t90 = fit.GetX(channelData.baseline - 0.9 * channelData.amplitude);
        channelData.fallTime = channelData.t50 - channelData.t10;
    }

    const int nSamples = 80;
    auto grMeanAndRMS = GetMeanAndRMS(*gr, 0, nSamples);
    channelData.RMS = grMeanAndRMS.second;

    // Debugging session
    /*
    if (channel == 1)
    {
        std::cout << std::endl;
        std::cout << GREEN << "Event: " << event << RESET << std::endl;
        std::cout << "Channel: " << channel << std::endl;
        std::cout << "MinimumPos: " << gr->GetPointX(minimumIndex) << std::endl;
        std::cout << "Minimum: " << minimum << std::endl;
        std::cout << "EdgeLeft: " << edgeLeft << std::endl;
        std::cout << "EdgeRight: " << edgeRight << std::endl;
        std::cout << "Baseline: " << channelData.baseline << std::endl;
        std::cout << "MinLevel: " << channelData.minLevel << std::endl;
        std::cout << "Amplitude: " << channelData.amplitude << std::endl;
        std::cout << "t10: " << channelData.t10 << std::endl;
        std::cout << "t50: " << channelData.t50 << std::endl;
        std::cout << "t90: " << channelData.t90 << std::endl;
        std::cout << "FallTime: " << channelData.fallTime << std::endl;
        std::cout << "RMS: " << channelData.RMS << std::endl;
    }
    */

    delete gr;
    inFile->Close();

    return 1;
}

void Preprocessor::DrawEvent(const int event, const int channel, const char * outFilePath)
{
    TString sOutFilePath(outFilePath);
    if (sOutFilePath.EqualTo("default"))    sOutFilePath = Form("ITS3/Data/Event%dChannel%d.root", event, channel);
    
    if (channel > fNChannels || event > fNEvents)
    {
        std::cerr << RED << "Error: channel or event out of range." << RESET << std::endl;
        return;
    }

    auto inFile = TFile::Open(fInFilePath.Data());
    TString grName = Form("grEv%dChanC%dsamp%d", event, channel, fSamplingPeriod);
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
double FindMinimumIndex(TGraph & gr)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    double minimum{0.};
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
                std::cout << "meanDer: " << meanDer << ", RMSDer: " << RMSDer << std::endl;
                std::cout << "y["<<i<<"]: " << y[i] << ", y["<<i<<"+1]: " << y[i+1] << ", y["<<i<<"+2]: " << y[i+2] << std::endl;
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
    std::cout << "] " << std::setw(3) << static_cast<int>(percentage * 100.0) << "%";
    std::cout << "  Elapsed Time: " << elapsedTime << "s";
    std::cout.flush();
}
// -------------------------------------------------------------------------------------
