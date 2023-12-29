/*
    Script to do the energy calibration of the APTS sensor.
    It converts amplitude in mV to deposited charge in electrons.

    The script takes as input the peak position of known structures of a radioactive isotope.
    These informations are stored in a log file that can be read by mini-yaml.
    The peak position is given in mV and it is converted to electrons using the energy calibration.

    The energy calibration is performed using the following formula:
        E = (V - [0]) / [1]
    where:
        E is the energy in eV
        V is the voltage in V
        [1] and [0] are the calibration parameters

    The calibration parameters are obtained by fitting the peak positions of the radioactive isotope.
    A point at zero is assumed. The validity of the assumption is then checked with a Gauss test.
*/

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TString.h>

#include "../yaml/Yaml.hpp"

void calibration()
{
    const char * inLogPath = "ITS3/Data/output/spectrumOutput.log";
    const char * outLogPath = "ITS3/Data/output/calibrationOutput.log";
    const char * outFilePath = "ITS3/Data/output/calibrationOutput.root";

    const char * inFilePath = "ITS3/Data/run175174828_230428174901_preprocessed.root";
    const char * inTreeName = "PreprocessedData";

    auto inFile = TFile::Open(inFilePath);
    auto inTree = (TTree*)inFile->Get(inTreeName);

    // read the log file
    Yaml::Node inLog;
    Yaml::Parse(inLog, inLogPath);

    // create the output file
    std::ofstream createOutput(outLogPath);
    createOutput.close();

    auto outFile = new TFile(outFilePath, "recreate");

    const int nPixels = 4;
    int pixels[] = {5, 6, 9, 10};
    double keV_to_electrons = 3.6e3;
    double peakPositions[] = {5.89875, 6.49045};                                                // in electrons
    double peakPositionsElectrons[] = {5.89875 * keV_to_electrons, 6.49045 * keV_to_electrons}; // in electrons
    double peakPositionsError[] = {0., 0.};                                                     // in electrons

    // loop over the pixels
    for (int i = 0; i < nPixels; i++)
    {
        const int iPixel = pixels[i];
        
        auto g = new TGraphErrors();
        g->SetName(Form("calibrationPx%d", iPixel));
        g->SetTitle(Form("Calibration for pixel %d", iPixel));
        g->GetXaxis()->SetTitle("Energy (keV)");
        g->GetYaxis()->SetTitle("Amplitude (mV)");

        auto ge = new TGraphErrors();
        ge->SetName(Form("calibrationPx%dElectrons", iPixel));
        ge->SetTitle(Form("Calibration for pixel %d in electrons", iPixel));
        ge->GetXaxis()->SetTitle("Charge (e^{-})");
        ge->GetYaxis()->SetTitle("Amplitude (mV)");

        TH1D * hRMS = new TH1D(Form("RMSpx%d", iPixel), "", 100, 0, 1);
        inTree->Draw(Form("pixel%d.RMS>>RMSpx%d", iPixel, iPixel));
        const double RMS = hRMS->GetMean();

        g->SetPoint(0, 0, 0);           // add the point at zero
        g->SetPointError(0, 0, RMS);      // add the point at zero
        ge->SetPoint(0, 0, 0);          // add the point at zero
        ge->SetPointError(0, 0, RMS);     // add the point at zero

        std::cout << "pixel" << iPixel << ":" << std::endl;
        for (int iPeak = 0; iPeak < inLog[Form("pixel%d", iPixel)]["nFits"].As<int>(); iPeak++)
        {
            Yaml::Node cfgPeak = inLog[Form("pixel%d", iPixel)]["fits"][iPeak];
            double peak = cfgPeak["peak"].As<double>();
            double peakError = cfgPeak["peakerror"].As<double>();
            
            g->SetPoint(iPeak + 1, peakPositions[iPeak], peak);
            g->SetPointError(iPeak + 1, peakPositionsError[iPeak], peakError);
            ge->SetPoint(iPeak + 1, peakPositionsElectrons[iPeak], peak);
            ge->SetPointError(iPeak + 1, peakPositionsError[iPeak] * keV_to_electrons, peakError);


            std::cout << "  peak" << iPeak << ":" << std::endl;
            std::cout << "    -  " << peak << " +/- " << peakError << std::endl;
            std::cout << "    -  " << peakPositions[iPeak] << " +/- " << peakPositionsError[iPeak] << std::endl;
            std::cout << "    -  " << peakPositionsElectrons[iPeak] << " +/- " << peakPositionsError[iPeak] * keV_to_electrons << std::endl;
            std::cout << std::endl;
        }

        outFile->cd();
        auto f = new TF1(Form("calibrationPx%dFit", iPixel), "[0] + [1]*x", 0, 10*keV_to_electrons);
        g->Fit(f, "RMQ+");
        g->Write();
        f->Write();

        auto fe = new TF1(Form("calibrationPx%dFitElectrons", iPixel), "[0] + [1]*x", 0, 10*keV_to_electrons);
        ge->Fit(fe, "RMQ+");
        ge->Write();
        fe->Write();

        // write fit results to log file
        std::streambuf* originalCoutBuffer = std::cout.rdbuf();
        std::ofstream outputFile(outLogPath, std::ios_base::app);
        std::cout.rdbuf(outputFile.rdbuf());

        std::cout << "pixel" << iPixel << ":" << std::endl;
        std::cout << "  calibrationPars:" << std::endl;
        std::cout << "    -  " << f->GetParameter(0) << std::endl;
        std::cout << "    -  " << f->GetParameter(1) << std::endl;
        std::cout << "  calibrationParsError:" << std::endl;
        std::cout << "    -  " << f->GetParError(0) << std::endl;
        std::cout << "    -  " << f->GetParError(1) << std::endl;
        std::cout << "  calibrationParsElectrons:" << std::endl;
        std::cout << "    -  " << fe->GetParameter(0) << std::endl;
        std::cout << "    -  " << fe->GetParameter(1) << std::endl;
        std::cout << "  calibrationParsErrorElectrons:" << std::endl;
        std::cout << "    -  " << fe->GetParError(0) << std::endl;
        std::cout << "    -  " << fe->GetParError(1) << std::endl;
        std::cout << "  chi2: " << f->GetChisquare() << std::endl;
        std::cout << "  NDF: " << f->GetNDF() << std::endl;
        std::cout << std::endl;

        std::cout.rdbuf(originalCoutBuffer);

        delete f;
        delete fe;
        delete g;
        delete ge;
        delete hRMS;
    }

    outFile->Close();
    inFile->Close();
}