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
   TH1D* histoAmpli2 = new TH1D("histoAmpli2","histoAmpli2",200,0.,1.); //mV, to check
   TH1D* histoAmpli3 = new TH1D("histoAmpli3","histoAmpli3",200,0.,1.);
   TH1D* histoToA2 = new TH1D("histoToA2","histoToA2",200,0.,3.); //ns, to check
   TH1D* histoToA3 = new TH1D("histoToA3","histoToA3",200,0.,3.);
   TH1D* histoRMS2 = new TH1D("histoRMS2","histoRMS2",200,0.,1.); //mV, to check
   TH1D* histoRMS3 = new TH1D("histoRMS3","histoRMS3",200,0.,1.);
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   
      cout<<"jentry="<<jentry<<endl;

      // Il vostro codice va qui
      const int numPoints=w2[jentry].size();
      cout<<"numPoints="<<numPoints<<endl;
      double baseLine2{0.}, baseLine3{0.};
      double RMS2{0.}, RMS3{0.};
      cout<<"cout5"<<endl;
      TH1D* histoNoise2 = new TH1D("histoNoise2","histoNoise2",200,-1.,1.);
      TH1D* histoNoise3 = new TH1D("histoNoise3","histoNoise3",200,-1.,1.);
      cout<<"cout45"<<endl;
      double amplitude2{0.}, amplitude3{0.};
      double toA2{0.}, toA3{0.};
      for(int ii=0;ii<numPoints;ii++)
      {
         //cout<<"ii="<<ii<<endl;
         if(ii>20 && ii<220)
         {
            histoNoise2->Fill(w2[jentry][ii]);
            //cout<<"cout3"<<endl;
            histoNoise3->Fill(w3[jentry][ii]);
         }
         if(w2[jentry][ii]>amplitude2)
         {
            amplitude2=w2[jentry][ii];
            toA2=t2[jentry][ii];
         }
             if(w3[jentry][ii]>amplitude3)
         {
            amplitude3=w3[jentry][ii];
            toA3=t3[jentry][ii];
         }
      }
      cout<<"cout2"<<endl;
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


      nbytes += nb;
      
   }

// Gli istogrammi vanno fittati e plottati qui

}
