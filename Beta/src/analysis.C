// Analysis class generated for Lab2 beta setup
#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>
#include <iostream>

void SetStyle(bool graypalette)
{
  //gStyle->Reset("Plain");
  //gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette)
    gStyle->SetPalette(8,0);
  else
    gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetHistLineWidth(1);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.1,"x");
  //gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //gStyle->SetTickLength(0.04,"X")
  //gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}


void analysis::Loop()
{
   SetStyle(true);
   if (fChain == 0) return;
   TH1D * histoEvent100_ch2; //mV
   TH1D * histoEvent100_ch3; //mV
   double Amp2{0.}, Amp3{0.};
   double ToA2{0.}, ToA3{0.};
   double RMS2{0.}, RMS3{0.};
   double baseline2{0.}, baseline3{0.};  
   bool FillTTree = true;                 
   
   double amplitude_cut2{50.}, amplitude_cut3{60.};

   TFile outFile("Beta/data/output/BetaOutput.root", "recreate");
   TTree * tree;
   if (FillTTree)
   {
      tree = new TTree("BetaTree","BetaTree");
      tree->Branch("Amp2", &Amp2);   
      tree->Branch("Amp3", &Amp3); 
      tree->Branch("RMS2", &RMS2); 
      tree->Branch("RMS3", &RMS3); 
      tree->Branch("ToA2", &ToA2); 
      tree->Branch("ToA3", &ToA3);       
      tree->Branch("baseline2", &baseline2); 
      tree->Branch("baseline3", &baseline3); 
   }
   

   TH1D * histoAmpli2 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
   TH1D * histoAmpli3 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
   TH1D * histoAmpli2_sel = new TH1D("histoAmpli2_sel","histoAmpli2_Sel",200,0.,400); //mV
   TH1D * histoAmpli3_sel = new TH1D("histoAmpli3_sel","histoAmpli3_Sel",200,0.,400);
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
      Amp2=0.;
      Amp3=0.;
      ToA2=0.;
      ToA3=0.;
      RMS2=0.;
      RMS3=0.;
      baseline2=0.;
      baseline3=0.;
      
      // Il vostro codice va qui
      const int numPoints = w2->size();
      TH1D* histoNoise2 = new TH1D("histoNoise2","histoNoise2",200,-1.,1.);
      TH1D* histoNoise3 = new TH1D("histoNoise3","histoNoise3",200,-1.,1.);
      if (jentry==100)
      {
         histoEvent100_ch2 = new TH1D("histoEvent100_ch2","histoEvent100_ch2",t2->size(),t2->front(),t2->back());
         histoEvent100_ch3 = new TH1D("histoEvent100_ch3","histoEvent100_ch3",t3->size(),t3->front(),t3->back());
      }
      for(int ii=0;ii<numPoints;ii++)
      {
         if (jentry==100)
         {
            histoEvent100_ch2->Fill(t2->at(ii),w2->at(ii));
            histoEvent100_ch3->Fill(t3->at(ii),w3->at(ii));
         }

         if(ii>20 && ii<220)
         {
            histoNoise2->Fill(w2->at(ii)*1000);
            histoNoise3->Fill(w3->at(ii)*1000);
         }
         if(w2->at(ii)*1000>Amp2)
         {
            Amp2=w2->at(ii)*1000;
            ToA2=t2->at(ii)*1000000000;
         }
         if(w3->at(ii)*1000>Amp3)
         {
            Amp3=w3->at(ii)*1000;
            ToA3=t3->at(ii)*1000000000;
         }
      }
      
      if(Amp2>amplitude_cut2 && Amp3>amplitude_cut3)
      {
         histoAmpli2_sel->Fill(Amp2);
         histoAmpli3_sel->Fill(Amp3);
         //Add code here
      }
      
      baseline2=histoNoise2->GetMean();
      baseline3=histoNoise3->GetMean();
      RMS2=histoNoise2->GetRMS();
      RMS3=histoNoise3->GetRMS();
      histoAmpli2->Fill(Amp2);
      histoAmpli3->Fill(Amp3);
      histoToA2->Fill(ToA2);
      histoToA3->Fill(ToA3);
      histoRMS2->Fill(RMS2);
      histoRMS3->Fill(RMS3);
      delete histoNoise2;
      delete histoNoise3;
      if (FillTTree)
         tree->Fill();
      nbytes += nb;
   }

   if (FillTTree)
   {
      tree->Write();
      outFile.Close();
   }
    
    
// Gli istogrammi vanno fittati e plottati qui
//TCanvas * canvasA = new TCanvas("canvasA","canvasA",2000,1000);
//canvasA->Divide(2,1);
//canvasA->cd(1);
//gPad->SetLogy();
//histoAmpli2->Draw("hist");
//histoAmpli2->GetXaxis()->SetTitle("Amplitude (mV)"); 
//histoAmpli2->GetYaxis()->SetTitle("Entries (a.u.)");
//histoAmpli2->SetTitle("Sensor A - Channel 2");
//canvasA->cd(2);
//gPad->SetLogy();
//histoAmpli3->Draw("hist");
//histoAmpli3->SetLineColor(kRed);
//histoAmpli3->GetXaxis()->SetTitle("Amplitude (mV)"); 
//histoAmpli3->GetYaxis()->SetTitle("Entries (a.u.)");
//histoAmpli3->SetTitle("Sensor B - Channel 3");
//canvasA->SaveAs("Beta/data/output/Amplitude.pdf");
//
//TCanvas * canvasA_Sel = new TCanvas("canvasA_Sel","canvasA_Sel",2000,1000);
//canvasA_Sel->Divide(2,1);
//canvasA_Sel->cd(1);
//gPad->SetLogy();
//histoAmpli2_sel->Draw("hist");
//histoAmpli2_sel->GetXaxis()->SetTitle("Amplitude (mV)"); 
//histoAmpli2_sel->GetYaxis()->SetTitle("Entries (a.u.)");
//histoAmpli2_sel->SetTitle("Sensor A - Channel 2");
//canvasA_Sel->cd(2);
//gPad->SetLogy();
//histoAmpli3_sel->Draw("hist");
//histoAmpli3_sel->SetLineColor(kRed);
//histoAmpli3_sel->GetXaxis()->SetTitle("Amplitude (mV)"); 
//histoAmpli3_sel->GetYaxis()->SetTitle("Entries (a.u.)");
//histoAmpli3_sel->SetTitle("Sensor B - Channel 3");
//canvasA_Sel->SaveAs("Beta/data/output/Amplitude_selected.pdf");
//
//TCanvas * canvasT = new TCanvas("canvasT","canvasT",1000,1000);
//canvasT->Divide(2,1);
//canvasT->cd(1);
//histoToA2->Draw("hist");
//histoToA2->GetXaxis()->SetTitle("Time of arrival (ns)"); 
//histoToA2->GetYaxis()->SetTitle("Entries (a.u.)");
//histoToA2->SetTitle("Sensor A - Channel 2");
//canvasT->cd(2);
//histoToA3->Draw("hist");
//histoToA3->SetLineColor(kRed);
//histoToA3->GetXaxis()->SetTitle("Time of arrival (ns)");
//histoToA3->GetYaxis()->SetTitle("Entries (a.u.)");
//histoToA3->SetTitle("Sensor B - Channel 3");
//
//TCanvas * canvasRMS = new TCanvas("canvasRMS","canvasRMS",1000,1000);
//canvasRMS->Divide(2,1);
//canvasRMS->cd(1);
//histoRMS2->Draw("hist");
//histoRMS2->GetXaxis()->SetTitle("Baseline RMS (mV)"); 
//histoRMS2->GetYaxis()->SetTitle("Entries (a.u.)");
//histoRMS2->SetTitle("Sensor A - Channel 2");
//canvasRMS->cd(2);
//histoRMS3->Draw("hist");
//histoRMS3->SetLineColor(kRed);
//histoRMS3->GetXaxis()->SetTitle("Baseline RMS (mV)"); 
//histoRMS3->GetYaxis()->SetTitle("Entries (a.u.)");
//histoRMS3->SetTitle("Sensor B - Channel 3");
//
//TCanvas * canvas100 = new TCanvas("canvas100","canvas100",1000,1000);
//double max = histoEvent100_ch2->GetMaximum();
//double min = histoEvent100_ch2->GetMinimum();
//histoEvent100_ch3->GetMaximum()>max ? max=histoEvent100_ch3->GetMaximum() : max=max;
//histoEvent100_ch3->GetMinimum()<min ? min=histoEvent100_ch3->GetMinimum() : min=min;
//double Xmin = histoEvent100_ch2->GetXaxis()->GetXmin();
//double Xmax = histoEvent100_ch2->GetXaxis()->GetXmax();
//histoEvent100_ch3->GetXaxis()->GetXmin()<Xmin ? Xmin=histoEvent100_ch3->GetXaxis()->GetXmin() : Xmin=Xmin;
//histoEvent100_ch3->GetXaxis()->GetXmax()>Xmax ? Xmax=histoEvent100_ch3->GetXaxis()->GetXmax() : Xmax=Xmax;
//canvas100->DrawFrame(Xmin,min*1.1,Xmax,max*1.1,"Event 100;t (s);V (mV)");
//histoEvent100_ch2->Draw("hist,same");
//histoEvent100_ch2->GetXaxis()->SetTitle("Baseline 100 (mV)"); 
//histoEvent100_ch2->GetYaxis()->SetTitle("Entries (a.u.)");
//histoEvent100_ch2->SetTitle("Sensor A - Channel 2");
//histoEvent100_ch3->Draw("hist,same");
//histoEvent100_ch3->SetLineColor(kRed);
//histoEvent100_ch3->GetXaxis()->SetTitle("Baseline 100 (mV)"); 
//histoEvent100_ch3->GetYaxis()->SetTitle("Entries (a.u.)");
//histoEvent100_ch3->SetTitle("Sensor B - Channel 3");
//TLegend * legend = new TLegend(0.2,0.7,0.4,0.8);
//legend->AddEntry(histoEvent100_ch2,"Sensor A - Channel 2","l");
//legend->AddEntry(histoEvent100_ch3,"Sensor B - Channel 3","l");
//legend->SetTextSize(0.03);
//legend->Draw();
//canvas100->Modified();
//canvas100->Update();
//canvas100->SaveAs("Beta/data/output/Event_100.pdf");
}