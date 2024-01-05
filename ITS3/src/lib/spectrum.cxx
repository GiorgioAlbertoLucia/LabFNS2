/*
    Class implementation of spectrum.hh
*/

#include <iostream>
#include <fstream>
#include <streambuf>

#include "../../../yaml/Yaml.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCanvas.h>

#include <RooRealVar.h> 
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHistPdf.h>
#include <RooPolynomial.h>

#include "spectrum.hh"

/**
 * @brief Construct a new Spectrum:: Spectrum object.
 * 
 */
Spectrum::Spectrum(TH1D * inSpectrum, TFile * outFile):
    fSpectrum(inSpectrum),
    fBackground(nullptr),
    fSignal(nullptr),
    fOutFile(outFile)
{
    fNBins = fSpectrum->GetNbinsX();
    fXmin = fSpectrum->GetXaxis()->GetXmin();
    fXmax = fSpectrum->GetXaxis()->GetXmax();
    
    //fFits = new std::vector<TF1>();
}

Spectrum::~Spectrum()
{
    delete fSpectrum;
    delete fBackground;
    delete fSignal;
}

void Spectrum::EstimateBackground()
{
    const int nBins = fSpectrum->GetNbinsX();
    const double xMin = fSpectrum->GetXaxis()->GetXmin();
    const double xMax = fSpectrum->GetXaxis()->GetXmax();
    
    double source[nBins];
    for (int iBin = 1; iBin < nBins; iBin++)    source[iBin] = fSpectrum->GetBinContent(iBin);

    auto spectrum = new TSpectrum();
    spectrum->Background(source, nBins, 15, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, true, TSpectrum::kBackSmoothing3, false);
    fBackground = new TH1D("background", "background", nBins, xMin, xMax);
    for (int iBin = 0; iBin < nBins; iBin++)
    {
        fBackground->SetBinContent(iBin, source[iBin]);
        fBackground->SetBinError(iBin, sqrt(source[iBin]));
    }

    fSignal = new TH1D("signal", "signal", nBins, xMin, xMax);
    fSignal->Add(fSpectrum, fBackground, 1, -1);
}

void Spectrum::SimpleFitSpectrum(const char * outLogPath, Yaml::Node & cfgFits)
{
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::ofstream outputFile(outLogPath, std::ios_base::app);
    std::cout.rdbuf(outputFile.rdbuf());

    std::cout << "  nFits: " << cfgFits["nFits"].As<int>() << std::endl;
    std::cout << "  fits:" << std::endl;

    for (int iFit = 0; iFit < cfgFits["nFits"].As<int>(); iFit++)
    {
        auto cfgFit = cfgFits["fits"][iFit];
        
        TF1 fit(cfgFit["name"].As<std::string>().c_str(), cfgFit["func"].As<std::string>().c_str(), cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());
        for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    fit.SetParameter(iParam, cfgFit["param"][iParam].As<double>());
        fSpectrum->Fit(&fit, "RM+", "", cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());
        
        // write fit results to log file
        std::cout << "    - params:" << std::endl;
        for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    std::cout << "        - " << fit.GetParameter(iParam) << std::endl;
        std::cout << "      parerrors:" << std::endl;
        for (int iParam = 0; iParam < cfgFit["nParams"].As<int>(); iParam++)    std::cout << "        - " << fit.GetParError(iParam) << std::endl;
        std::cout << "      chi2: " << fit.GetChisquare() << std::endl;
        std::cout << "      NDF: " << fit.GetNDF() << std::endl;
        std::cout << "      peak: " << fit.GetParameter(1) << std::endl;
        std::cout << "      peakerror: " << fit.GetParError(1) << std::endl; 
        std::cout << std::endl;

        fit.SetLineColor(cfgFit["color"].As<int>());
        fit.SetRange(cfgFit["range"][0].As<double>(), cfgFit["range"][1].As<double>());
        fFits.push_back(fit);
    }

    std::cout.rdbuf(originalCoutBuffer);

    // save spectrum to file 
    fOutFile->cd();
    fSpectrum->Write();
}

void Spectrum::FitSpectrum(const char * cfgPath, const int nFits)
{

    RooRealVar x("x", "x", fXmin, fXmax);
    x.setRange("fitRange", 62, 72);

    RooDataHist data("data", "data", RooArgList(x), fSpectrum);
    RooDataHist dataBkg("dataBkg", "dataBkg", RooArgList(x), fBackground);

    RooRealVar c0("c0", "c0", 1, -40, 40);
    RooRealVar c1("c1", "c1", 1, -40, 40);
    RooRealVar c2("c2", "c2", 1, -40, 40);

    RooArgList coeffs(c0, c1, c2);
    RooPolynomial bkg("bkg", "bkg", x, coeffs);

    bkg.fitTo(dataBkg, RooFit::Range("fitRange"), RooFit::SumW2Error(true));

    RooPlot* frame = x.frame();
    dataBkg.plotOn(frame);
    bkg.plotOn(frame);

    TCanvas canvas("canvas", "Background Fit Result");
    frame->Draw();
    canvas.SaveAs("background_fit_result.png");

    // SIGNAL FIT

    RooRealVar mean("mean", "mean of Gaussian", 68, 60, 80);
    RooRealVar sigma("sigma", "width of Gaussian", 2, 0.1, 10);
    RooGaussian signal("signal", "Gaussian Signal", x, mean, sigma);

    // Create a RooRealVar to represent the fraction of signal in the combined fit
    RooRealVar nsig("nsig", "number of signal events", 100, 0, 10000);
    RooRealVar nbkg("nbkg", "number of background events", 500, 0, 10000);

    // Create a RooAddPdf for the combined fit
    RooAddPdf model_combined("model_combined", "Combined Signal + Background Model", RooArgList(signal, bkg), RooArgList(nsig, nbkg));

    // Perform the combined fit
    model_combined.fitTo(data, RooFit::SumW2Error(true));

    // Plot the result
    RooPlot* frame_combined = x.frame();
    data.plotOn(frame_combined);
    model_combined.plotOn(frame_combined);
    model_combined.plotOn(frame_combined, RooFit::Components("signal"), RooFit::LineStyle(kDashed));
    model_combined.plotOn(frame_combined, RooFit::Components("bkg"), RooFit::LineStyle(kDotted));

    TCanvas canvas1("canvas1", "Fit Result");
    frame_combined->Draw();
    canvas1.SaveAs("combined_fit_result.png");


    /*    
    RooRealVar x("x", "x", fXmin, fXmax);

    RooDataHist data("data", "data", RooArgList(x), fSpectrum);
    RooDataHist dataBkg("dataBkg", "dataBkg", RooArgList(x), fBackground);


    for (int iFit = 0; iFit < nFits; iFit++)
    {
        x.setRange("fitRange", 62, 72);

        // definition of RooVars
        RooRealVar mean("mean", "mean", 69.5, 68, 70);
        RooRealVar sigma("sigma", "sigma", 3, 0, 10);
        
        RooRealVar c0("c0", "c0", 1, 0, 10);
        RooRealVar c1("c1", "c1", 1, 0, 10);
        RooRealVar c2("c2", "c2", 1, 0, 10);

        RooRealVar nsig("nsig", "nsig", 100, 0, 1000);
        RooRealVar nbkg("nbkg", "nbkg", 100, 0, 1000);

        RooGaussian signal("signal", "signal", x, mean, sigma);
        RooArgList coeffs(c0, c1, c2);
        RooPolynomial bkg("bkg", "bkg", x, coeffs);

        auto bkgFitRes = bkg.fitTo(dataBkg, RooFit::Extended(true), RooFit::Range("fitRange"));
        bkgFitRes->Print();

        RooHistPdf bkgPdf("bkgPdf", "bkgPdf", RooArgSet(x), dataBkg, 0);
        //RooAddPdf signalAndBkg("signalAndBkg", "signalAndBkg", RooArgList(signal, bkgPdf), RooArgList(nsig, nbkg));
//
        //auto fitRes = signalAndBkg.fitTo(data, RooFit::Extended(true), RooFit::Range("fitRange"));
        //fitRes->Print();

        RooPlot* frame = x.frame();
        data.plotOn(frame);
        //signalAndBkg.plotOn(frame);
        //signalAndBkg.plotOn(frame, RooFit::Components("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
        bkgPdf.plotOn(frame);
        bkgPdf.plotOn(frame, RooFit::Components("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));

        TCanvas canvas("canvas", "Fit Result");
        frame->Draw();
        canvas.SaveAs("fit_result.png");

        //fFits.push_back(model);
    }

    */
    
}

void Spectrum::DrawSpectrumAndBackground() const
{
    
    fOutFile->cd();

    fSpectrum->SetLineColor(kOrange-3);
    fBackground->SetLineColor(kBlue-3);
    fSignal->SetLineColor(kRed-3);

    fSpectrum->Write();
    fBackground->Write();
    fSignal->Write();

    TCanvas canvas("spectrum&bkg", "Spectrum and Background", 800, 600);
    fSpectrum->Draw("hist");
    fBackground->Draw("hist same");
    canvas.Write();
}