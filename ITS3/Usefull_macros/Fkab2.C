#include <TFile.h>
#include "Riostream.h"
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TAttLine.h>
#include <TStyle.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <TString.h>
#include <TPaveStats.h>
using namespace std;
void NormalizeHistos(TH1F* histo){
  /*double tot=histo->Integral();
      if(tot>0){
	histo->Scale(1./tot);
	histo->GetYaxis()->SetTitle("Entries (a.u.)");
      }*/
      Double_t factor = 1.;
histo->Scale(factor/histo->Integral("width"));
    }

void DoFit(TH1F* h, TF1* funcKa, TF1* funcKb){
  double roughmean=h->GetBinCenter(h->GetMaximumBin());
  double roughpeak=h->GetBinContent(h->GetMaximumBin());
  double roughbase=h->GetBinContent(h->GetXaxis()->FindBin(65.));
  funcKa->SetParameters(roughbase,roughpeak,roughmean,2);
  h->Fit(funcKa, "RM+");
  int peakb=-1;
  double cntpeakb=0;
  for(int jb=1; jb<h->GetNbinsX(); jb++){
    double x=h->GetBinCenter(jb);
    double y=h->GetBinContent(jb);
    if(x>funcKa->GetParameter(2)+3*funcKa->GetParameter(3) && y>cntpeakb){
      peakb=jb;
      cntpeakb=y;
    }
  }
  funcKb->SetParameters(h->GetBinContent(peakb),h->GetBinCenter(peakb),funcKa->GetParameter(3));
  h->Fit(funcKb,"","",funcKa->GetParameter(2)+3*funcKa->GetParameter(3),100);
  return;
}

void ComputeLinearParams(TGraphErrors* g, double& slop, double& cnst){
  double x1,y1,x2,y2;
  g->GetPoint(0,x1,y1);
  g->GetPoint(1,x2,y2);
  slop=(y2-y1)/(x2-x1);
  cnst=y1-slop*x1;
  return;
}

void Fkab2(TString filename="APTS10_Vbb-4V_WaltConf_trgOR_20220311_TTree.root",string filename2="DerivativeCalibration_APTS06.txt", double maxFallTime=1000.){//quando si esegue va messo vb con l'intero del Vbb
/*double gain[5][16] =   {{0.6022, 0.5691, 0.5864, 0.6044, 0.5737, 0.5342, 0.5454, 0.5819, 0.5750, 0.5163, 0.5501, 0.5872, 0.5653, 0.5973, 0.5962, 0.5667},
                        {0.6583, 0.6230, 0.6435, 0.6555, 0.6261, 0.5750, 0.5916, 0.6335, 0.6289, 0.5628, 0.5977, 0.6444, 0.6232, 0.6533, 0.6492, 0.6280},
                        {0.6782, 0.6471, 0.6635, 0.6787, 0.6496, 0.5932, 0.6123, 0.6585, 0.6520, 0.5801, 0.6213, 0.6681, 0.6532, 0.6739, 0.6742, 0.6582},
                        {0.6927, 0.6571, 0.6710, 0.6829, 0.6578, 0.605338, 0.62469, 0.6686, 0.6617, 0.588382, 0.624109, 0.6761, 0.6581, 0.6805, 0.6852, 0.6709},
                        {0.6903, 0.6584, 0.6716, 0.6823, 0.6616, 0.6015, 0.6121, 0.6707, 0.6641, 0.5855, 0.6280, 0.6825, 0.6639, 0.6826, 0.6895, 0.6762}};*/
  double eA=((5.89875*100)+(5.88765*51))/151;// energie dei picchi, per la Ka uso la media pesata dei due
  double eB=6.49045;
  double SikA = 1.74; //Silicon KAlpha peak
  double FeEscape = 4.15;
  
  
  TH1F* hAmplChan[4];
  TH1F* hAmplChaneq[4];
  TH1F* hAmplChanCluSiz1[4];
  TH1F* hAmplChanFastCluSiz1[4];
  TH1F* hAmplChantot;
  TH1F* hAmplChantoten;


  TFile* f=new TFile(filename.Data());
    //double p[5][4];
 /*for(int ii=0;ii<5;ii++)
 {
  for(int gg=0;gg<16;gg++)
  {
    if(gg==5) p[ii][gg]=gain[ii][gg];
    if(gg==6) p[ii][gg]=gain[ii][gg];
    if(gg==9) p[ii][gg]=gain[ii][gg];
    if(gg==10) p[ii][gg]=gain[ii][gg];
    
  }
 }
    for(int t=0;t<16;t++)
    {
      if(t==5) cout<<"gain= "<<gain[vb][t]<<endl;
      if(t==6) cout<<"gain= "<<gain[vb][t]<<endl;
      if(t==9) cout<<"gain= "<<gain[vb][t]<<endl;
      if(t==10) cout<<"gain= "<<gain[vb][t]<<endl;
      if(t<4) cout<<"p= "<<p[vb][t]<<endl;
      
      }*/
   //fattore di calibrazione per ogni canale
  double q[4];
  double pratio[4];
  double qratio[4];
  double p_TOTAL=0;
  double q_TOTAL=0;
  double Ampleq[4];
  double p[4]; //fattore di calibrazione per ogni canale

  double alf=0.,bet=0.;  stringstream infile(filename2.c_str());
   ifstream input(infile.str().c_str(), ifstream::in);
  int k=0;
  string line;

  //ciclo per la lettura del file e memorizzo i dati 
  while(getline(input, line)){
    stringstream read(line);
    read>>p[k]>>q[k]>>p_TOTAL>>q_TOTAL;
    cout<<p[k]<<" "<<q[k]<<" "<<p_TOTAL<<" "<<q_TOTAL<<endl; 
    k++;
    }
  
  TTree* tree=(TTree*)f->Get("treeParams");
  double amplVec[4];
  double fallTimeVec[4];
  double t50[4];
  double t10[4];

   
   hAmplChantot=new TH1F ("hAmplChantot"," All clusters ; Signal Amplitude tot (mV) ; Entries",200.,0.,140.);//isto di ampl tot dopo equalizzazione
   hAmplChantoten=new TH1F ("hAmplChantoten"," All clusters ; Energy (KeV) ; Entries",100,0.,10.);//isto in energia
   int w=0;
  for(int k=0; k<4; k++){
    if(k==0) w=5;
    if(k==1) w=6;
    if(k==2) w=9;
    if(k==3) w=10;
    tree->SetBranchAddress(Form("FallTimePx%d",w),&fallTimeVec[k]);
    tree->SetBranchAddress(Form("SignalAmplPx%d",w),&amplVec[k]);
    tree->SetBranchAddress(Form("t50Px%d",w),&t50[k]);
    tree->SetBranchAddress(Form("t10Px%d",w),&t10[k]);
    hAmplChan[k]=new TH1F(Form("hAmplChan%d",k+1),Form(" All clusters ; Signal Amplitude Px%d (mV) ; Entries",w),200,8.,130.);//ampiezza per canale non eq
    hAmplChaneq[k]=new TH1F(Form("hAmplChaneq%d",k+1),Form(" All clusters ; Signal Amplitude equalized Px%d (mV) ; Entries",w),200,10.,150.);//ampiezze per canale dopo eq

    //hAmplChanCluSiz1[k]=new TH1F(Form("hAmplChanCluSiz1%d",k+1),Form(" Cluster Size = 1 ; Signal Amplitude Ch%d (mV) ; Entries",k+1),80,60.,100.);
    //hAmplChanFastCluSiz1[k]=new TH1F(Form("hAmplChanFastCluSiz1%d",k+1),Form("Cluster Size = 1, FallTime < %.0f; Signal Amplitude Ch%d (mV) ; Entries",maxFallTime,k+1),80,60.,100.);
  }
  for(int ient=0; ient<tree->GetEntriesFast(); ient++){
    tree->GetEvent(ient);
    //int clusiz=0;
    double maxsig=-999.;
    int maxsigpix=-1;
    
    for(int k=0; k<4; k++){//prima era tutto 
    if(k==0) w=5;
    if(k==1) w=6;
    if(k==2) w=9;
    if(k==3) w=10;
      if(amplVec[k]>8. && (t50[k]-t10[k])*2<1000.){
        hAmplChantot->Fill(amplVec[k]/p[k]);
  	//clusiz+=1;
//cout<<amplVec[k]<<endl;
	hAmplChan[k]->Fill(amplVec[k]);
  hAmplChaneq[k]->Fill(amplVec[k]/p[k]);
 
	if(amplVec[k]>maxsig){
	  maxsig=amplVec[k];
	  maxsigpix=k;
	}
      }
   
    /*if(clusiz==1){
      hAmplChanCluSiz1[maxsigpix]->Fill(amplVec[maxsigpix]*1000.);
      if(fallTimeVec[maxsigpix]<maxFallTime) hAmplChanFastCluSiz1[maxsigpix]->Fill(amplVec[maxsigpix]);
    }*/
  }
  }

  //TF1* funcKa=new TF1("funcKa","[0]*(1-TMath::Erf((x-[2])/[3]))+[1]*TMath::Exp(-0.5*(x-[2])*(x-[2])/[3]/[3])",60,90.);
  //TF1* funcKb=new TF1("funcKb","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",80,100.);
 
  TF1* funcKa=new TF1("funcKa","[0]*(1-TMath::Erf((x-[2])/[3]))+[1]*TMath::Exp(-0.5*(x-[2])*(x-[2])/[3]/[3])",45.,50.);
  //dipende dal Vbb e cambia molto da pixel a pixel tantissimo, conviene farlo dal canvas
  TF1* funcKb=new TF1("funcKb","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",50.,55.);//dipende dal Vbb e cambia molto da pixel a pixel tantissimo, conviene farlo dal canvas
   
  TF1* funcKaeq=new TF1("funcKaeq","[0]*(1-TMath::Erf((x-[2])/[3]))+[1]*TMath::Exp(-0.5*(x-[2])*(x-[2])/[3]/[3])",40.,45.);
  //Vbb:0=40,45| 1,2=75, 82| |3.6= 100,120
  TF1* funcKbeq=new TF1("funcKbeq","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",45.,50.);
  //vbb 0=45,50| 1.2=| |3.6=120,130 
  TF1* funcKaeqt=new TF1("funcKaeqt","[0]*(1-TMath::Erf((x-[2])/[3]))+[1]*TMath::Exp(-0.5*(x-[2])*(x-[2])/[3]/[3])",100.,120.);
  //Vbb: 1,2=75, 82| |3.6= 100,120

  TF1* funcKbeqt=new TF1("funcKbeqt","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",120.,132.);
  //vbb 1.2=| |3.6=120,130 
  TF1* pikSi1=new TF1("piKsi1","[3]*x+[4]+[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",18.,30.);// questi due non vengono se non si taglia su ft circa a 320 ps
  //Vbb= 3.6->23,40
  TF1* pikSi2=new TF1("piKsi2","[3]*x+[4]+[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",48.,61.);
//Vbb= 3.6->75,85
  pikSi1->SetParameter(1,28.);//Vbb: 3.6->28
  pikSi1->SetParameter(2,0.8);//Vbb: 3.6->0.8
  
  pikSi2->SetParameter(1,80.);//Vbb: 3.6->80
  pikSi2->SetParameter(2,0.5);//Vbb: 3.6->0.5
  funcKa->SetLineWidth(3);
  funcKa->SetLineColor(kRed+1);
  funcKb->SetLineWidth(3);
  funcKb->SetLineColor(4);

  funcKaeq->SetLineWidth(2);
  funcKaeq->SetLineColor(kRed+1);
  funcKbeq->SetLineWidth(3);
  funcKbeq->SetLineColor(4);

  funcKaeqt->SetLineWidth(2);
  funcKaeqt->SetLineColor(kRed+1);
  //funcKbeqt->SetParameter(1,125.);
  //funcKbeqt->SetParameter(2,1.);
  

  funcKbeqt->SetLineWidth(3);
  funcKbeqt->SetLineColor(4);

  pikSi1->SetLineWidth(3);
  pikSi1->SetLineColor(kGreen+1);
 
  pikSi2->SetLineWidth(3);
  pikSi2->SetLineColor(kBlack);
 


  
  TGraphErrors* gAall[4];
  TGraphErrors* gAalleq[4];
  //TGraphErrors* gAclu1[4];
  //TGraphErrors* gAfast[4];
  TGraphErrors* calEn= new TGraphErrors(3);
  TGraphErrors* calEn4= new TGraphErrors(5);
  int cols[4]={1,kRed+1,kGreen+2,kBlue};
  for(int k=0; k<4; k++){
    gAall[k]=new TGraphErrors(0);
    gAalleq[k]=new TGraphErrors(0);
    //gAclu1[k]=new TGraphErrors(0);
    //gAfast[k]=new TGraphErrors(0);
    gAall[k]->SetMarkerColor(cols[k]);
    gAall[k]->SetLineColor(cols[k]);
    gAall[k]->SetMarkerStyle(25);
    gAall[k]->SetMarkerSize(1.2);
    gAalleq[k]->SetMarkerColor(cols[k]);
    gAalleq[k]->SetLineColor(cols[k]);
    gAalleq[k]->SetMarkerStyle(25);
    gAalleq[k]->SetMarkerSize(1.2);
    //gAclu1[k]->SetMarkerColor(cols[k]);
    //gAclu1[k]->SetLineColor(cols[k]);
    //gAclu1[k]->SetMarkerStyle(20);
    //gAfast[k]->SetMarkerColor(cols[k]);
   // gAfast[k]->SetLineColor(cols[k]);
   // gAfast[k]->SetMarkerStyle(47);
  }
  
  double slopeAll[4]; //questo non so cosa faccia,ma l'ho lasciato
  double constAll[4];
  double slopeAlleq[4];  
  double constAlleq[4];
  //double slopeClu1[4];  
  //double constClu1[4];
  //double slopeFast[4];  
  //double constFast[4];
    for(int k=0;k<4;k++)
  {
  NormalizeHistos(hAmplChan[k]);
  NormalizeHistos(hAmplChaneq[k]);
  }
  NormalizeHistos(hAmplChantot);
   NormalizeHistos(hAmplChantoten);
  TCanvas* cA = new TCanvas("cA","Amplitudes",1650,900);
  cA->Divide(2,2);
  for(int k=0; k<4; k++){
    cA->cd(k+1);
    DoFit(hAmplChan[k],funcKa,funcKb);
     cout<<"*************** non eq"<<endl;
    hAmplChan[k]->Draw("HIST");
    funcKa->DrawClone("same");
    funcKb->DrawClone("same");
    gAall[k]->SetPoint(0,eA,funcKa->GetParameter(2));
    gAall[k]->SetPointError(0,0.,funcKa->GetParError(2));
    gAall[k]->SetPoint(1,eB,funcKb->GetParameter(1));
    gAall[k]->SetPointError(1,0.,funcKb->GetParError(1));
    ComputeLinearParams(gAall[k],slopeAll[k],constAll[k]);
  }
  
  TCanvas* cAeq=new TCanvas("cAeq","Amplitudes equalizated",1650,900);
  cAeq->Divide(2,2);
  for(int k=0; k<4; k++){
    cAeq->cd(k+1);
    DoFit(hAmplChaneq[k],funcKaeq,funcKbeq);
    cout<<"*************  eq"<<endl;
    hAmplChaneq[k]->Draw("HIST");
    funcKaeq->DrawClone("same");
    funcKbeq->DrawClone("same");
     gAalleq[k]->SetPoint(0,eA,funcKaeq->GetParameter(2));
    gAalleq[k]->SetPointError(0,0.,funcKaeq->GetParError(2));
    gAalleq[k]->SetPoint(1,eB,funcKbeq->GetParameter(1));
    gAalleq[k]->SetPointError(1,0.,funcKbeq->GetParError(1));
    ComputeLinearParams(gAalleq[k],slopeAlleq[k],constAlleq[k]);
    cout<<"valori medi dei picchi:alfa "<<funcKaeq->GetParameter(2)<<" beta: "<<funcKbeq->GetParameter(1)<<endl;
  }


   TCanvas* cAeqt=new TCanvas("cAeqt","Amplitudes equalizated Tot",1650,900);
 
  
    cAeqt->cd();
    DoFit(hAmplChantot,funcKaeqt,funcKbeqt);
    cout<<"**************** eq tot"<<endl;
    hAmplChantot->Draw("HIST");
    //funcKaeqt->DrawClone("same");
    //funcKbeqt->DrawClone("same");
    cout<<"**************  silicio"<<endl;
    //hAmplChantot->Fit(pikSi1,"RM+");
   // hAmplChantot->Fit(pikSi2,"RM+");
    //pikSi1->DrawClone("same");
    //pikSi2->DrawClone("same");

    calEn->SetPoint(0,0.,0.);
    calEn->SetPoint(1,funcKaeqt->GetParameter(2),eA);
    calEn->SetPointError(1,funcKaeqt->GetParError(2),0.);
    calEn->SetPoint(2,funcKbeqt->GetParameter(1),eB);
    calEn->SetPointError(2,funcKbeqt->GetParError(1),0.);
    calEn->SetTitle("Calibrazione in energia");

    calEn->GetXaxis()->SetTitle("Amplitude [mV]");
    calEn->GetYaxis()->SetTitle("Energy [KeV]");
    calEn->SetMarkerSize(0.8);
    calEn->SetMarkerStyle(21);
     TCanvas* cr=new TCanvas("cr","retta di calibrazioen in energia",1650,900);
    cr->cd();
    calEn->Draw("AP");


    TF1* retcal=new TF1("retcal","[0]*x+[1]",0.,130.);
    retcal->SetLineColor(94);
    retcal->SetParameter(0,0.06);
    calEn->Fit(retcal,"RM+");

    //calEn4->SetPoint(0,0.,0.);
    calEn4->SetPoint(0,funcKaeqt->GetParameter(2),eA);
    calEn4->SetPointError(0,funcKaeqt->GetParError(2),0.);
    calEn4->SetPoint(2,funcKbeqt->GetParameter(1),eB);
    calEn4->SetPointError(2,funcKbeqt->GetParError(1),0.);
    calEn4->SetPoint(3,pikSi1->GetParameter(1),SikA);
    calEn4->SetPointError(3,pikSi1->GetParError(1),0.);
    calEn4->SetPoint(1,pikSi2->GetParameter(1),FeEscape);
    calEn4->SetPointError(1,pikSi2->GetParError(1),0.);
    calEn4->SetTitle("Calibrazione in energia vs2");
    calEn4->SetPoint(4,0.,0.);
    calEn4->SetPointError(4,0.,0.);

    calEn4->GetXaxis()->SetTitle("Amplitude [mV]");
    calEn4->GetYaxis()->SetTitle("Energy [KeV]");
    calEn4->SetMarkerSize(0.8);
    calEn4->SetMarkerStyle(21);
    calEn4->GetXaxis()->SetRangeUser(0.,130.);
    TCanvas* cr1=new TCanvas("cr1","retta di calibrazioen in energia",1650,900);
    cr1->cd();
    calEn4->Draw("AP");

      TF1* retcal4=new TF1("retcal4","[0]*x+[1]",5.,130.);//questa non viene a meno di tagli sul ft a 320 ps
    retcal4->SetLineColor(63);
    retcal4->SetParameter(0,0.06);
    calEn4->Fit(retcal4,"RM+");
    double d=retcal4->GetParameter(0);
    //printf("d vale %f",d);
    for(int ient=0; ient<tree->GetEntriesFast(); ient++){
      tree->GetEvent(ient); 
      for(int k=0; k<4; k++){
      if(k==0) w=5;
    if(k==1) w=6;
    if(k==2) w=9;
    if(k==3) w=10;
        if(amplVec[k]>8. && (t50[k]-t10[k])*2<1500.){
         hAmplChantoten->Fill((amplVec[k]/p[k])*d);
      }
    }
   }

   TCanvas* cE=new TCanvas("cE","Energy",1650,900);// viene solo se viene retcal4 quindi solo se si taglia sa ft
    cE->cd();
    hAmplChantoten->SetLineColor(4);
    hAmplChantoten->Draw();

   



    
    /* gAalleq[k]->SetPoint(0,eA,funcKaeq->GetParameter(2));
    gAalleq[k]->SetPointError(0,0.,funcKaeq->GetParError(2));
    gAalleq[k]->SetPoint(1,eB,funcKbeq->GetParameter(1));
    gAalleq[k]->SetPointError(1,0.,funcKbeq->GetParError(1));
    ComputeLinearParams(gAalleq[k],slopeAlleq[k],constAlleq[k]);*/
  

  //cA->SaveAs("Fit-Kab-AllClu.png");
  
  /*TCanvas* c1 = new TCanvas("c1","Amplitudes CluSiz1",1650,900);
  c1->Divide(2,2);
  for(int k=0; k<4; k++){
    c1->cd(k+1);
    DoFit(hAmplChanCluSiz1[k],funcKa,funcKb);
    hAmplChanCluSiz1[k]->Draw();
    funcKa->DrawClone("same");
    funcKb->DrawClone("same");
    gAclu1[k]->SetPoint(0,eA,funcKa->GetParameter(2));
    gAclu1[k]->SetPointError(0,0.,funcKa->GetParError(2));
    gAclu1[k]->SetPoint(1,eB,funcKb->GetParameter(1));
    gAclu1[k]->SetPointError(1,0.,funcKb->GetParError(1));
    ComputeLinearParams(gAclu1[k],slopeClu1[k],constClu1[k]);
  }*/
  //c1->SaveAs("Fit-Kab-CluSiz1.png");

  /*TCanvas* cF = new TCanvas("cF","Amplitudes Fast CluSiz1",1650,900);
  cF->Divide(2,2);
  for(int k=0; k<4; k++){
    cF->cd(k+1);
    DoFit(hAmplChanFastCluSiz1[k],funcKa,funcKb);
    hAmplChanFastCluSiz1[k]->Draw();
    funcKa->DrawClone("same");
    funcKb->DrawClone("same");
    gAfast[k]->SetPoint(0,eA,funcKa->GetParameter(2));
    gAfast[k]->SetPointError(0,0.,funcKa->GetParError(2));
    gAfast[k]->SetPoint(1,eB,funcKb->GetParameter(1));
    gAfast[k]->SetPointError(1,0.,funcKb->GetParError(1));
    ComputeLinearParams(gAfast[k],slopeFast[k],constFast[k]);
  }*/
  //cF->SaveAs("Fit-Kab-Fast.png");

  
  }