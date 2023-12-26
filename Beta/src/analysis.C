// Analysis class generated for Lab2 beta setup
#define analysis_cxx

#include "analysis.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>
#include <TF1.h>
#include <THStack.h>
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

void setHistos(TH1D * histo, const char * title, const char * xtitle, const char * ytitle, const int color)
{
    histo->SetTitle(title);
    histo->GetXaxis()->SetTitle(xtitle);
    histo->GetYaxis()->SetTitle(ytitle);
    histo->SetFillColor(color);
    histo->SetLineColor(color);
}


void analysis::Loop(const bool FillTTree = true)
{
    SetStyle(true);
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    TH1D * histoEvent100_ch2; //mV
    TH1D * histoEvent100_ch3; //mV
    double Amp2{0.}, Amp3{0.};
    double ToA2{0.}, ToA3{0.};
    double ToA2f{0.}, ToA3f{0.};
    double RMS2{0.}, RMS3{0.};
    double baseline2{0.}, baseline3{0.};       
    int xx{0}, yy{0};           
    
    double amplitude_cut2{50.}, amplitude_cut3{60.};

    TH1D * histoAmpli2 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
    TH1D * histoAmpli3 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
    TH1D * histoAmpli2A0 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
    TH1D * histoAmpli3A0 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
    TH1D * histoAmpli2A1 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
    TH1D * histoAmpli3A1 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
    TH1D * histoAmpli2A2 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
    TH1D * histoAmpli3A2 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
    TH1D * histoAmpli2A3 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
    TH1D * histoAmpli3A3 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
    TH1D * histoAmpli2A4 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
    TH1D * histoAmpli3A4 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
    TH1D * histoAmpli2A5 = new TH1D("histoAmpli2","histoAmpli2",200,0.,400); //mV
    TH1D * histoAmpli3A5 = new TH1D("histoAmpli3","histoAmpli3",200,0.,400);
    TH1D * histoAmpli2_sel = new TH1D("histoAmpli2_sel","histoAmpli2_Sel",200,0.,400); //mV
    TH1D * histoAmpli3_sel = new TH1D("histoAmpli3_sel","histoAmpli3_Sel",200,0.,400);
    TH1D * histoToA2 = new TH1D("histoToA2","histoToA2",25,-0.025,1.225); //ns
    TH1D * histoToA3 = new TH1D("histoToA3","histoToA3",25,-0.025,1.225);
    TH1D * histoRMS2 = new TH1D("histoRMS2","histoRMS2",150,0.3,0.8); //mV
    TH1D * histoRMS3 = new TH1D("histoRMS3","histoRMS3",150,0.3,0.8);
    TH1D * histoChiSquare2 = new TH1D("histoChiSquare2","histoChiSquare2",nentries,0.,nentries);
    TH1D * histoChiSquare3 = new TH1D("histoChiSquare3","histoChiSquare3",nentries,0.,nentries);

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
        tree->Branch("ToA2f", &ToA2f); 
        tree->Branch("ToA3f", &ToA3f);      
        tree->Branch("baseline2", &baseline2); 
        tree->Branch("baseline3", &baseline3); 
    }
   

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   
        Amp2=0.;
        Amp3=0.;
        ToA2=0.;
        ToA3=0.;
        ToA2f=0.;
        ToA3f=0.;
        RMS2=0.;
        RMS3=0.;
        baseline2=0.;
        baseline3=0.;
        xx=0;
        yy=0;

        
        const int numPoints = w2->size();
        TH1D* histoNoise2 = new TH1D("histoNoise2","histoNoise2",200,-1.,1.);   // histograms to measure rms
        TH1D* histoNoise3 = new TH1D("histoNoise3","histoNoise3",200,-1.,1.);
        TH1D* histoTime2 = new TH1D("histoTime2","histoTime2",500,-1.5,2.5);     // histograms to measure ToA
        TH1D* histoTime3 = new TH1D("histoTime3","histoTime3",500,-1.5,2.5);

        if (jentry==100)
        {
            histoEvent100_ch2 = new TH1D("histoEvent100_ch2","histoEvent100_ch2",t2->size(),t2->front(),t2->back());
            histoEvent100_ch3 = new TH1D("histoEvent100_ch3","histoEvent100_ch3",t3->size(),t3->front(),t3->back());
        }
        for(int ii=0;ii<numPoints;ii++)
        {
            histoTime2->Fill(1.e9*(t2->at(ii)),1000*(w2->at(ii)));
            histoTime3->Fill(1.e9*(t3->at(ii)),1000*(w3->at(ii)));
            
            if (jentry==100)
            {
                histoEvent100_ch2->Fill(t2->at(ii),w2->at(ii)*1000);
                histoEvent100_ch3->Fill(t3->at(ii),w3->at(ii)*1000);
            }

            if(ii>20 && ii<220)
            {
                histoNoise2->Fill(w2->at(ii)*1000);
                histoNoise3->Fill(w3->at(ii)*1000);
            }
        }

        xx = findIndexOfMaxElement(w2);
        yy = findIndexOfMaxElement(w3);
        Amp2 = w2->at(xx)*1000;
        Amp3 = w3->at(yy)*1000;

        ToA2 = t2->at(xx)*1.e9;
        ToA3 = t3->at(yy)*1.e9;
        if(jentry==100 || jentry==200 || jentry==300) cout<<"tempo 2: "<<ToA2<<" tempo 3: "<<ToA3<<endl;
        if(jentry==100 || jentry==200 || jentry==300) cout<<"tempo 2bis: "<<t2->at(xx)<<" tempo 3bis: "<<t3->at(yy)<<endl;

        if(Amp2>amplitude_cut2 && Amp3>amplitude_cut3)
        {
            histoAmpli2_sel->Fill(Amp2);
            histoAmpli3_sel->Fill(Amp3);
            
            TF1* fit2 = new TF1("fit2","gaus",1.e9*(t2->at(xx-10)),1.e9*(t2->at(xx+10)));
            TF1* fit3 = new TF1("fit3","gaus",1.e9*(t3->at(yy-10)),1.e9*(t3->at(yy+10)));
            FitSignal(histoTime2, histoTime3, ToA2f, ToA3f, fit2, fit3);
            if (fit2->GetChisquare()/fit2->GetNDF() == fit2->GetChisquare()/fit2->GetNDF())
                histoChiSquare2->Fill(jentry, fit2->GetChisquare()/fit2->GetNDF());
            if (fit3->GetChisquare()/fit3->GetNDF() == fit3->GetChisquare()/fit3->GetNDF())
                histoChiSquare3->Fill(jentry, fit3->GetChisquare()/fit3->GetNDF());
                
            if(jentry==100 ){
                TCanvas * canvas1 = new TCanvas("canvas1","canvas1",1000,700);
                canvas1->Divide(2,1);
                canvas1->cd(1);
                histoTime2->SetLineColor(kBlue -10);
                histoTime2->SetTitle("Sensor A - Channel 2");
                histoTime2->GetXaxis()->SetTitle("Time (ns)");
                histoTime2->GetYaxis()->SetTitle("Amplitude (mV)");
                histoTime2->Draw("hist");
                fit2->Draw("same");
                fit2->SetLineColor(kRed);
                fit2->SetLineWidth(2);
                gStyle->SetOptFit(11111111);
                canvas1->cd(2);
                histoTime3->SetLineColor(kRed -10);
                histoTime3->SetTitle("Sensor B - Channel 3");
                histoTime3->GetXaxis()->SetTitle("Time (ns)");
                histoTime3->GetYaxis()->SetTitle("Amplitude (mV)");
                histoTime3->Draw("hist");
                fit3->Draw("same");
                fit3->SetLineColor(kRed);
                fit3->SetLineWidth(2);
                gStyle->SetOptFit(11111111);
                cout<<"chi2/dof fit2: "<<fit2->GetChisquare()<<"/"<<fit2->GetNDF()<<endl;
                cout<<"chi2/dof fit3: "<<fit3->GetChisquare()<<"/"<<fit3->GetNDF()<<endl;
                canvas1->SaveAs("Beta/data/output/checkfit.pdf");
            }
            delete fit2;
            delete fit3;
        }
        
        else
        {
            ToA2f = ToA2;
            ToA3f = ToA3;
        }
        if(jentry==100 || jentry==200 || jentry==300) cout<<"tempo 2: "<<ToA2<<" tempo 2 con fit: "<<ToA2f<<endl;
        if(jentry==100 || jentry==200 || jentry==300) cout<<"tempo 3: "<<ToA3<<" tempo 3 con fit: "<<ToA3f<<endl;

        baseline2=histoNoise2->GetMean();
        baseline3=histoNoise3->GetMean();
        RMS2=histoNoise2->GetRMS();
        RMS3=histoNoise3->GetRMS();
        histoAmpli2->Fill(Amp2);
        histoAmpli3->Fill(Amp3);
        if(jentry<nentries/6)
        { 
            if(Amp2>15) histoAmpli2A0->Fill(Amp2);
            if(Amp3>15) histoAmpli3A0->Fill(Amp3);
        }
        if(jentry>nentries/6 && jentry<2*nentries/6)
        { 
            if(Amp2>15) histoAmpli2A1->Fill(Amp2);
            if(Amp3>15) histoAmpli3A1->Fill(Amp3);
        }
        if(jentry>2*nentries/6 && jentry<3*nentries/6)
        { 
            if(Amp2>15) histoAmpli2A2->Fill(Amp2);
            if(Amp3>15) histoAmpli3A2->Fill(Amp3);
        }
        if(jentry>3*nentries/6 && jentry<4*nentries/6)
        { 
            if(Amp2>15) histoAmpli2A3->Fill(Amp2);
            if(Amp3>15) histoAmpli3A3->Fill(Amp3);
        }
        if(jentry>4*nentries/6 && jentry<5*nentries/6)
        { 
            if(Amp2>15) histoAmpli2A4->Fill(Amp2);
            if(Amp3>15) histoAmpli3A4->Fill(Amp3);
        }
        if(jentry>5*nentries/6 && jentry<6*nentries/6)
        { 
            if(Amp2>15) histoAmpli2A5->Fill(Amp2);
            if(Amp3>15) histoAmpli3A5->Fill(Amp3);
        }
        histoToA2->Fill(ToA2);
        histoToA3->Fill(ToA3);
        histoRMS2->Fill(RMS2);
        histoRMS3->Fill(RMS3);

        delete histoNoise2;
        delete histoNoise3;
        delete histoTime2;
        delete histoTime3;

        if (FillTTree)  tree->Fill();
        nbytes += nb;
   }

    histoChiSquare2->Write();
    histoChiSquare3->Write();

   if (FillTTree)
   {
      tree->Write();
      outFile.Close();
   }
    
   // -----------------------------------------
   //             plotting section
   if(!FillTTree)
   {
      TCanvas * canvasA = new TCanvas("canvasA","canvasA",2000,1000);
      canvasA->Divide(2,1);
      canvasA->cd(1);
      gPad->SetLogy();
      histoAmpli2->Draw("hist");
      histoAmpli2->GetXaxis()->SetTitle("Amplitude (mV)"); 
      histoAmpli2->GetYaxis()->SetTitle("Entries (a.u.)");
      histoAmpli2->SetTitle("Sensor A - Channel 2");
      canvasA->cd(2);
      gPad->SetLogy();
      histoAmpli3->Draw("hist");
      histoAmpli3->SetLineColor(kRed);
      histoAmpli3->GetXaxis()->SetTitle("Amplitude (mV)"); 
      histoAmpli3->GetYaxis()->SetTitle("Entries (a.u.)");
      histoAmpli3->SetTitle("Sensor B - Channel 3");
      canvasA->SaveAs("Beta/data/output/Amplitude.pdf");

      TCanvas * canvasA_Sel = new TCanvas("canvasA_Sel","canvasA_Sel",2000,1000);
      canvasA_Sel->Divide(2,1);
      canvasA_Sel->cd(1);
      gPad->SetLogy();
      histoAmpli2_sel->Draw("hist");
      histoAmpli2_sel->GetXaxis()->SetTitle("Amplitude (mV)"); 
      histoAmpli2_sel->GetYaxis()->SetTitle("Entries (a.u.)");
      histoAmpli2_sel->SetTitle("Sensor A - Channel 2");
      canvasA_Sel->cd(2);
      gPad->SetLogy();
      histoAmpli3_sel->Draw("hist");
      histoAmpli3_sel->SetLineColor(kRed);
      histoAmpli3_sel->GetXaxis()->SetTitle("Amplitude (mV)"); 
      histoAmpli3_sel->GetYaxis()->SetTitle("Entries (a.u.)");
      histoAmpli3_sel->SetTitle("Sensor B - Channel 3");
      canvasA_Sel->SaveAs("Beta/data/output/Amplitude_selected.pdf");

      TCanvas * canvasT = new TCanvas("canvasT","canvasT",1000,1000);
      canvasT->Divide(2,1);
      canvasT->cd(1);
      histoToA2->Draw("hist");
      histoToA2->GetXaxis()->SetTitle("Time of arrival (ns)"); 
      histoToA2->GetYaxis()->SetTitle("Entries (a.u.)");
      histoToA2->SetTitle("Sensor A - Channel 2");
      TF1 * fitToA2 = new TF1("fitToA2","gaus",0.5,0.9);
      histoToA2->Fit(fitToA2,"rm+");
      fitToA2->SetLineColor(kBlue+2);
      fitToA2->Draw("same");
      gStyle->SetOptFit(1111);
      TLegend * legend1 = new TLegend(0.2,0.7,0.4,0.8);
      legend1->AddEntry(histoToA2,"Sensor A - Channel 2","l");
      legend1->AddEntry(fitToA2,"p0*exp(#frac{-(x-p1)^{2}}{2*p2^{2}})","l");
      legend1->SetTextSize(0.03);
      legend1->Draw("same");
      canvasT->cd(2);
      histoToA3->Draw("hist");
      histoToA3->SetLineColor(kRed);
      histoToA3->GetXaxis()->SetTitle("Time of arrival (ns)");
      histoToA3->GetYaxis()->SetTitle("Entries (a.u.)");
      histoToA3->SetTitle("Sensor B - Channel 3");
      TF1 * fitToA3 = new TF1("fitToA3","gaus",0.4,0.8);
      histoToA3->Fit(fitToA3,"rm+");
      fitToA3->SetLineColor(kRed+2);
      fitToA3->Draw("same");
      gStyle->SetOptFit(1111);
      TLegend * legend2 = new TLegend(0.15,0.7,0.35,0.8);
      legend2->AddEntry(histoToA3,"Sensor B - Channel 3","l");
      legend2->AddEntry(fitToA3,"p0*exp(#frac{-(x-p1)^{2}}{2*p2^{2}})","l");
      legend2->SetTextSize(0.03);
      legend2->Draw("same");
      cout<<"chi2/dof fit2: "<<fitToA2->GetChisquare()<<"/"<<fitToA2->GetNDF()<<endl;
      cout<<"chi2/dof fit3: "<<fitToA3->GetChisquare()<<"/"<<fitToA3->GetNDF()<<endl;
      canvasT->SaveAs("Beta/data/output/time.pdf");

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

      TCanvas * canvas100 = new TCanvas("canvas100","canvas100",1000,1000);
      double max = histoEvent100_ch2->GetMaximum();
      double min = histoEvent100_ch2->GetMinimum();
      histoEvent100_ch3->GetMaximum()>max ? max=histoEvent100_ch3->GetMaximum() : max=max;
      histoEvent100_ch3->GetMinimum()<min ? min=histoEvent100_ch3->GetMinimum() : min=min;
      double Xmin = histoEvent100_ch2->GetXaxis()->GetXmin();
      double Xmax = histoEvent100_ch2->GetXaxis()->GetXmax();
      histoEvent100_ch3->GetXaxis()->GetXmin()<Xmin ? Xmin=histoEvent100_ch3->GetXaxis()->GetXmin() : Xmin=Xmin;
      histoEvent100_ch3->GetXaxis()->GetXmax()>Xmax ? Xmax=histoEvent100_ch3->GetXaxis()->GetXmax() : Xmax=Xmax;
      canvas100->DrawFrame(Xmin,min*1.1,Xmax,max*1.1,"Event 100;t (s);V (mV)");
      histoEvent100_ch2->Draw("hist,same");
      histoEvent100_ch2->GetXaxis()->SetTitle("Baseline 100 (mV)"); 
      histoEvent100_ch2->GetYaxis()->SetTitle("Entries (a.u.)");
      histoEvent100_ch2->SetTitle("Sensor A - Channel 2");
      histoEvent100_ch3->Draw("hist,same");
      histoEvent100_ch3->SetLineColor(kRed);
      histoEvent100_ch3->GetXaxis()->SetTitle("Baseline 100 (mV)"); 
      histoEvent100_ch3->GetYaxis()->SetTitle("Entries (a.u.)");
      histoEvent100_ch3->SetTitle("Sensor B - Channel 3");
      TLegend * legend = new TLegend(0.6,0.65,0.8,0.8);
      legend->AddEntry(histoEvent100_ch2,"Sensor A - Channel 2","l");
      legend->AddEntry(histoEvent100_ch3,"Sensor B - Channel 3","l");
      legend->SetTextSize(0.05);
      legend->Draw();
      canvas100->Modified();
      canvas100->Update();
      canvas100->SaveAs("Beta/data/output/Event_100.pdf");

        TCanvas * canvasAmpliScan = new TCanvas("canvasAmpliScan","canvasAmpliScan",2000,1000);
        canvasAmpliScan->Divide(2,1);
        canvasAmpliScan->cd(1);
        cout<<" aaaaa2"<<endl;
        setHistos(histoAmpli2A0,"Sensor A - Channel 2","Amplitude (mV)","Entries (a.u.)",kBlue-10);
        setHistos(histoAmpli2A1,"Sensor A - Channel 2","Amplitude (mV)","Entries (a.u.)",kBlue-7);
        setHistos(histoAmpli2A2,"Sensor A - Channel 2","Amplitude (mV)","Entries (a.u.)",kBlue-4);
        setHistos(histoAmpli2A3,"Sensor A - Channel 2","Amplitude (mV)","Entries (a.u.)",kBlue-1);
        setHistos(histoAmpli2A4,"Sensor A - Channel 2","Amplitude (mV)","Entries (a.u.)",kBlue+2);
        setHistos(histoAmpli2A5,"Sensor A - Channel 2","Amplitude (mV)","Entries (a.u.)",kBlue+4);
        cout<<" aaaaa3"<<endl;
        auto hs1 = new THStack("hs1","Sensor A - Channel 2");
        cout<<" aaaaa4"<<endl;
        hs1->Add(histoAmpli2A0);
        hs1->Add(histoAmpli2A1);
        hs1->Add(histoAmpli2A2);
        hs1->Add(histoAmpli2A3);
        hs1->Add(histoAmpli2A4);
        hs1->Add(histoAmpli2A5);
        //hs1->GetXaxis()->SetTitle("Amplitude (mV)");
        //hs1->GetYaxis()->SetTitle("Entries (a.u.)");
        //hs1->GetXaxis()->SetRangeUser(20,400);
        hs1->Draw("hist");
        //hs1->GetXaxis()->SetTitle("Amplitude (mV)");  
        //hs1->GetYaxis()->SetTitle("Entries (a.u.)"); 
        cout<<" aaaaa"<<endl;
        /*histoAmpli2A0->Draw("hist");
        histoAmpli2A1->Draw("hist,same");
        histoAmpli2A2->Draw("hist,same");
        histoAmpli2A3->Draw("hist,same");
        histoAmpli2A4->Draw("hist,same");
        histoAmpli2A5->Draw("hist,same");*/
        TLegend * legendA = new TLegend(0.5,0.5,0.8,0.8);
        legendA->AddEntry(histoAmpli2A0,"Events 0 - 5385","f");
        legendA->AddEntry(histoAmpli2A1,"Events 5385 - 10771","f");
        legendA->AddEntry(histoAmpli2A2,"Events 10771 - 16154 ","f");
        legendA->AddEntry(histoAmpli2A3,"Events 16154 - 21539","f");
        legendA->AddEntry(histoAmpli2A4,"Events 21539 - 26992 ","f");
        legendA->AddEntry(histoAmpli2A5,"Events 26992 - 32315 ","f");
        legendA->SetTextSize(0.03);
        legendA->Draw();
        cout<<" aaaaadddd"<<endl;
        canvasAmpliScan->cd(2);
        setHistos(histoAmpli3A0,"Sensor B - Channel 3","Amplitude (mV)","Entries (a.u.)",kRed-10);
        setHistos(histoAmpli3A1,"Sensor B - Channel 3","Amplitude (mV)","Entries (a.u.)",kRed-7);
        setHistos(histoAmpli3A2,"Sensor B - Channel 3","Amplitude (mV)","Entries (a.u.)",kRed-4);
        setHistos(histoAmpli3A3,"Sensor B - Channel 3","Amplitude (mV)","Entries (a.u.)",kRed-1);
        setHistos(histoAmpli3A4,"Sensor B - Channel 3","Amplitude (mV)","Entries (a.u.)",kRed+2);
        setHistos(histoAmpli3A5,"Sensor B - Channel 3","Amplitude (mV)","Entries (a.u.)",kRed+4);
        //histoAmpli3A0->Draw("hist");
        /*histoAmpli3A1->Draw("hist,same");
        histoAmpli3A2->Draw("hist,same");
        histoAmpli3A3->Draw("hist,same");
        histoAmpli3A4->Draw("hist,same");
        histoAmpli3A5->Draw("hist,same");*/
        auto hs2 = new THStack("hs2","Sensor B - Channel 3");
        hs2->Add(histoAmpli3A0);
        hs2->Add(histoAmpli3A1);
        hs2->Add(histoAmpli3A2);
        hs2->Add(histoAmpli3A3);
        hs2->Add(histoAmpli3A4);
        hs2->Add(histoAmpli3A5);
        //hs2->GetXaxis()->SetTitle("Amplitude (mV)");
        //hs2->GetYaxis()->SetTitle("Entries (a.u.)");
        hs2->Draw("hist");
        //hs2->GetXaxis()->SetTitle("Amplitude (mV)");  
        //hs2->GetYaxis()->SetTitle("Entries (a.u.)"); 
        //hs2->GetXaxis()->SetRangeUser(20,400);
        TLegend * legendB = new TLegend(0.5,0.55,0.8,0.85);
        legendB->AddEntry(histoAmpli3A0,"Events 0 - 5385","f");
        legendB->AddEntry(histoAmpli3A1,"Events 5385 - 10771","f");
        legendB->AddEntry(histoAmpli3A2,"Events 10771 - 16154","f");
        legendB->AddEntry(histoAmpli3A3,"Events 16154 - 21539","f");
        legendB->AddEntry(histoAmpli3A4,"Events 21539 - 26992","f");
        legendB->AddEntry(histoAmpli3A5,"Events 26992 - 32315","f");
        legendB->SetTextSize(0.03);
        legendB->Draw();
        cout<<" aaaddffddfsda"<<endl;
        canvasAmpliScan->SaveAs("Beta/data/output/AmplitudeScan.pdf");
        cout<<" aafgsfhdaaa"<<endl;
   }    
   
}