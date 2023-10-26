// Analysis class generated for Lab2 beta setup
#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void analysis::Loop()
{

   if (fChain == 0) return;
   TH1D * histoAmpli2 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
   TH1D * histoAmpli3 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
   TH1D * histoToA2 = new TH1D("histoToA2","histoToA2",200,0.,2.5); //ns
   TH1D * histoToA3 = new TH1D("histoToA3","histoToA3",200,0.,2.5);
   TH1D * histoRMS2 = new TH1D("histoRMS2","histoRMS2",200,0.,1); //mV
   TH1D * histoRMS3 = new TH1D("histoRMS3","histoRMS3",200,0.,1);
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   

      // Il vostro codice va qui
      const int numPoints=w2->size();
      double baseLine2{0.}, baseLine3{0.};
      double RMS2{0.}, RMS3{0.};
      TH1D* histoNoise2 = new TH1D("histoNoise2","histoNoise2",200,-1.,1.);
      TH1D* histoNoise3 = new TH1D("histoNoise3","histoNoise3",200,-1.,1.);
      double amplitude2{0.}, amplitude3{0.};
      double toA2{0.}, toA3{0.};
      for(int ii=0;ii<numPoints;ii++)
      {
         if(ii>20 && ii<220)
         {
            histoNoise2->Fill(w2->at(ii)*1000);
            //cout<<"cout3"<<endl;
            histoNoise3->Fill(w3->at(ii)*1000);
         }
         if(w2->at(ii)*1000>amplitude2)
         {
            amplitude2=w2->at(ii)*1000;
            toA2=t2->at(ii)*1000000000;
         }
             if(w3->at(ii)*1000>amplitude3)
         {
            amplitude3=w3->at(ii)*1000;
            toA3=t3->at(ii)*1000000000;
         }
      }
      
      baseLine2=histoNoise2->GetMean();
      baseLine3=histoNoise3->GetMean();
      RMS2=histoNoise2->GetRMS();
      RMS3=histoNoise3->GetRMS();
      histoAmpli2->Fill(amplitude2);
      histoAmpli3->Fill(amplitude3);
      histoToA2->Fill(toA2);
      histoToA3->Fill(toA3);
      histoRMS2->Fill(RMS2);
      histoRMS3->Fill(RMS3);
      delete histoNoise2;
      delete histoNoise3;

      nbytes += nb;
   }

// Gli istogrammi vanno fittati e plottati qui
TCanvas * canvasA = new TCanvas("canvasA","canvasA",1000,1000);
canvasA->Divide(2,1);
canvasA->cd(1);
histoAmpli2->Draw("hist");
histoAmpli2->GetXaxis()->SetTitle("Amplitude (mV)"); 
histoAmpli2->GetYaxis()->SetTitle("Entries (a.u.)");
histoAmpli2->SetTitle("Sensor A - Channel 2");
canvasA->cd(2);
histoAmpli3->Draw("hist");
histoAmpli3->SetLineColor(kRed);
histoAmpli3->GetXaxis()->SetTitle("Amplitude (mV)"); 
histoAmpli3->GetYaxis()->SetTitle("Entries (a.u.)");
histoAmpli3->SetTitle("Sensor B - Channel 3");

TCanvas * canvasT = new TCanvas("canvasT","canvasT",1000,1000);
canvasT->Divide(2,1);
canvasT->cd(1);
histoToA2->Draw("hist");
histoToA2->GetXaxis()->SetTitle("Time of arrival (ns)"); 
histoToA2->GetYaxis()->SetTitle("Entries (a.u.)");
histoToA2->SetTitle("Sensor A - Channel 2");
canvasT->cd(2);
histoToA3->Draw("hist");
histoToA3->SetLineColor(kRed);
histoToA3->GetXaxis()->SetTitle("Time of arrival (ns)");
histoToA3->GetYaxis()->SetTitle("Entries (a.u.)");
histoToA3->SetTitle("Sensor B - Channel 3");

TCanvas * canvasRMS = new TCanvas("canvasRMS","canvasRMS",1000,1000);
canvasRMS->Divide(2,1);
canvasRMS->cd(1);
histoRMS2->Draw("hist");
histoRMS2->GetXaxis()->SetTitle("Baseline RMS (mV)"); 
histoRMS2->GetYaxis()->SetTitle("Entries (a.u.)");
histoRMS2->SetTitle("Sensor A - Channel 2");
canvasRMS->cd(2);
histoRMS3->Draw("hist");
histoRMS3->SetLineColor(kRed);
histoRMS3->GetXaxis()->SetTitle("Baseline RMS (mV)"); 
histoRMS3->GetYaxis()->SetTitle("Entries (a.u.)");
histoRMS3->SetTitle("Sensor B - Channel 3");


}
