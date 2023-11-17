#include <TFile.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMultiGraph.h>
#include <TKey.h>
#include <TList.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
  
 
 
  
  #include<string>
  #include<TString.h>
  
 
 
  
  #include<TH1F.h>
using namespace std;

TGraph* Smooth(TGraph* g, int nsm=2){
  TGraph* gsm=new TGraph(0);
  int npts=0;
  for(int j=nsm; j<g->GetN()-nsm; j++){
    double x,y,x1,y1;
    g->GetPoint(j,x,y);
    double sum=y;
    for(int k=1; k<=nsm;k++){
      g->GetPoint(j-k,x1,y1);
      sum+=y1;
      g->GetPoint(j+k,x1,y1);
      sum+=y1;
    }
    sum/=(2*nsm+1);
    gsm->SetPoint(npts++,x,sum);
  }
  return gsm;
}

double GetMaxX(TGraph* g){
  double xmax=0.;
  for(int j=0; j<g->GetN(); j++){
    double x,y;
    g->GetPoint(j,x,y);
    if(x>xmax) xmax=x;
  }
  return xmax;
}

double ComputeDerivative(TGraph* g, int j, int npts=5){
  if(npts==3){
    double xm1,ym1,xp1,yp1;
    g->GetPoint(j-1,xm1,ym1);
    g->GetPoint(j+1,xp1,yp1);
    double der=(yp1-ym1)/(xp1-xm1);
    return der;
  }else{
    double xm2,ym2,xm1,ym1,xp1,yp1,xp2,yp2;
    g->GetPoint(j-2,xm2,ym2);
    g->GetPoint(j-1,xm1,ym1);
    g->GetPoint(j+1,xp1,yp1);
    g->GetPoint(j+2,xp2,yp2);
    double der=(ym2-8*ym1+8*yp1-yp2)/(xp2-xp1)/12.;
    return der;
  }
}

TGraph* GetDerivative(TGraph* g, int nsm=4){
  TGraph* gder=new TGraph(0);
  int npts=0;
  for(int j=2; j<g->GetN()-nsm-2; j++){
    double x,y;
    g->GetPoint(j,x,y);
    double der=0;
    double nnn=0;
    for(int k=0; k<=nsm;k++){
      der+=ComputeDerivative(g,j+k);
      nnn+=1.;
    }
    if(nnn>0){
      der/=nnn;
      gder->SetPoint(npts++,x,der);
    }
  }
  return gder;
}

TGraph* CountNextNegativeDer(TGraph* g){
  TGraph* gn=new TGraph(0);
  for(int j=0; j<g->GetN()-1; j++){
    double x,y;
    g->GetPoint(j,x,y);
    int cntneg=0;
    for(Int_t k=j; k<g->GetN()-1; k++){
      double der=ComputeDerivative(g,k);
      if(der>=0) break;
      else cntneg++;
    }
    gn->SetPoint(j,x,cntneg);
  }
  return gn;
}

void GetMeanAndRMSCounts(TGraph* g, double xmin, double xmax, double& mean, double& rms){
  double sum=0,sum2=0,cnts=0;
  for(int j=0; j<g->GetN(); j++){
    double x,c;
    g->GetPoint(j,x,c);
    if(x>xmin && x<xmax){
      cnts+=1.;
      sum+=c;
      sum2+=(c*c);
    }
  }
  if(cnts>0){
    mean=sum/cnts;
    rms=TMath::Sqrt(sum2/cnts-mean*mean);
  }else{
    mean=0;
    rms=0;
  }
  return;
}
double FindOnGraph(TGraph* gcount, double y, double xmin, double xmax, int interpolate, bool backw=kFALSE){
  int jfirst=0;
  int dstep=1;
  if(backw){
    jfirst=gcount->GetN();
    dstep=-1;
  }
  printf("   level=%f. xmin=%f  xmax=%f\n",y,xmin,xmax);
  for(int jstep=0; jstep<gcount->GetN(); jstep++){
    int j=jfirst+dstep*jstep;
    double x,c,xbef,cbef,xaft,caft,xaft2,caft2;
    gcount->GetPoint(j,x,c);
    gcount->GetPoint(j-dstep,xbef,cbef);
    gcount->GetPoint(j+dstep,xaft,caft);
    gcount->GetPoint(j+2*dstep,xaft2,caft2);
    if((dstep==1 && c<y && cbef>y && caft<y) || (dstep==-1 && c>y && cbef<y && caft>y) ){
      if(interpolate==0) return x;
      else{
	double sumx=0,sumx2=0,sumy=0,sumxy=0,npts=0;
	for(int k=j-interpolate; k<=j+interpolate; k++){
	  double xP,yP;
	  gcount->GetPoint(k,xP,yP);
	  sumx+=xP;
	  sumy+=yP;
	  sumxy+=(xP*yP);
	  sumx2+=(xP*xP);
	  npts+=1;
	}
	double m=(npts*sumxy-sumx*sumy)/(npts*sumx2-sumx*sumx);
	double q=(sumy*sumx2-sumx*sumxy)/(npts*sumx2-sumx*sumx);
	double xinterp=(y-q)/m;
	if(xinterp<xmin || xinterp>xmax || TMath::Abs(xinterp-x)>1000.) continue;
	return xinterp;
      }
    }
  }
  return -999.;
}



void FindEdge(TGraph* gcount,TGraph* gnegder, TGraph* gder, double& endplateau, double& edgeleft, double& edgeright){
  
  // first very rough: compute flat levels on the left and on the right and check their difference
  double maxTime=GetMaxX(gcount);
  double levleft,rmsleft,levright,rmsright; 
  printf("maxTime=%f\n",maxTime);
  GetMeanAndRMSCounts(gcount,0.,2000.,levleft,rmsleft);
  GetMeanAndRMSCounts(gcount,maxTime-2000,maxTime,levright,rmsright);
  double y50=0.5*(levleft+levright);
  printf("y50=%f\n",y50);
  double t50fromleft=FindOnGraph(gcount,y50,0.,maxTime,4);
  double t50fromright=FindOnGraph(gcount,y50,0.,maxTime,4,kTRUE);
  double roughsig=levleft-levright;
  printf("roughsig=%f\n",roughsig);
  printf("Rough signal = %f Rough edge position = %f %f\n",roughsig,t50fromleft,t50fromright);
  double minSearchWindow=0;
  double maxSearchWindow=maxTime;
  if(roughsig>0.0005){
    minSearchWindow=TMath::Min(t50fromleft,t50fromright)-6000.;
    if(minSearchWindow<0) minSearchWindow=0;
    maxSearchWindow=TMath::Max(t50fromleft,t50fromright)+6000.;
    if(maxSearchWindow>maxTime) maxSearchWindow=maxTime;
  }
  printf("Search window = %f %f\n",minSearchWindow,maxSearchWindow);
  
  // second step: search for accumulation of adjacent points with negative derivative
  double xmaxn=-1;
  double cmaxn=-1;
  int jmaxn=-1;
  if(gnegder){
    for(int j=0; j<gnegder->GetN(); j++){
      double x,c;
      gnegder->GetPoint(j,x,c);
      if(x<minSearchWindow || x>maxSearchWindow) continue;
      if(c>cmaxn){
	cmaxn=c;
	xmaxn=x;
	jmaxn=j;
      }
      if(c==cmaxn){
	int sum0=0;
	int sum1=0;
	for(int k=1; k<20; k++){
	  double xk,ck;
	  gnegder->GetPoint(jmaxn+k,xk,ck);
	  sum0+=ck;
	  gnegder->GetPoint(j+k,xk,ck);
	  sum1+=ck;
	}
	if(sum1>sum0){
	  cmaxn=c;
	  xmaxn=x;
	  jmaxn=j;
	}
      }
    }
    printf("Maximum adjacent points with negative derivative: t_maxn=%f   n_neg=%f\n",xmaxn,cmaxn);
  }
  
  // third step: search for minimum of derivative and range where derivative differs from 0
  double xminder=-1;
  double dermin=99999.;
  int jminder=-1;
  for(int j=0; j<gder->GetN(); j++){
    double x,d;
    gder->GetPoint(j,x,d);
    if(x<minSearchWindow || x>maxSearchWindow) continue;
    if(d<dermin){
      dermin=d;
      xminder=x;
      jminder=j;
    }
  }
  if(jminder<0){
    endplateau=0;
    edgeleft=0;
    edgeright=0;
    return;
  }
  printf("Minimum of derivative: xminder=%f   dermin=%f\n",xminder,dermin);
  int jleft=-1;
  double dthresh=-1e-7;
  for(int j=jminder; j>0; j--){
    double x,d;
    gder->GetPoint(j,x,d);
    if(d>dthresh){
      jleft=j;
      break;
    }
  }
  int jright=-1;
  for(int j=jminder; j<gder->GetN(); j++){
    double x,d;
    gder->GetPoint(j,x,d);
    if(d>dthresh){
      jright=j;
      break;
    }
  }
  double xleft,xright,dum;
  gder->GetPoint(jleft,xleft,dum);
  gder->GetPoint(jright,xright,dum);
  printf("Region of negative derivative: xleft=%f   xright=%f\n",xleft,xright);
  if(xmaxn>0 && TMath::Abs(xmaxn-xleft)<5000 && TMath::Abs(xmaxn-xright)<5000){
    if(xleft>xmaxn) xleft=xmaxn;
    if(xright<xmaxn) xright=xmaxn;
  }
  printf("Edge range after analysis of derivative: xleft=%f   xright=%f\n",xleft,xright);

  // Fourth step: start from left and seach for N points with couns < baseline-3sigma
  double cmean,crms;
  GetMeanAndRMSCounts(gcount,0.,xleft,cmean,crms);
  printf("Mean before edge = %f rms = %f\n",cmean,crms);
  double thresh=cmean-3*crms;
  int threshbin=TMath::Nint(TMath::Max(10.,cmaxn/3.));
  if(cmaxn<0) threshbin=3;
  double xleft2=-1.;
  for(int j=0; j<gcount->GetN(); j++){
    double x,c;
    gcount->GetPoint(j,x,c);
    int nbelow=0;
    for(Int_t k=j; k<gcount->GetN()-1; k++){
      double x2,c2;
      gcount->GetPoint(k,x2,c2);
      if(c2<thresh) nbelow++;
      else break;
    }
    if(nbelow>threshbin){
      xleft2=x;
      break;
    }
  }
  printf("Left Edge from baseline-N*rms = %f\n",xleft2);
  if(xleft2>0){
    endplateau=TMath::Min(xleft,xleft2);
    edgeleft=TMath::Max(xleft,xleft2);
    edgeright=xright;
    printf("Edge range after all steps: endplateau=%f   edgeleft=%f   edgeright=%f\n",endplateau,edgeleft,edgeright);
    return;
  }
  endplateau=0;
  edgeleft=0;
  edgeright=0;
  return;
}

void ProcessEvent(TGraph* g, double params[20], bool plot){

  TGraph* gs=Smooth(g,10);
  gs->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  gs->SetTitle(g->GetTitle());
  gs->GetYaxis()->SetTitle("Amplitude (smoothened)");
  TGraph* gnegd=CountNextNegativeDer(g);
  gnegd->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  gnegd->SetTitle(g->GetTitle());
  gnegd->GetYaxis()->SetTitle("N. adjacent samplings with negative derivative");
  TGraph* gsnegd=CountNextNegativeDer(gs);
  TGraph* gd=GetDerivative(gs,40);
  gd->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  gd->SetTitle(g->GetTitle());
  gd->GetYaxis()->SetTitle("Amplitude derivative (smoothened)");
  double endpl,edgeleft,edgeright;
  printf("--- Edge finding on fine graph ---\n");
  FindEdge(g,gnegd,gd,endpl,edgeleft,edgeright);
  TF1* fbas=new TF1("fbas","[0]");
  TF1* flow=new TF1("flow","[0]");
  fbas->SetLineColor(2);
  flow->SetLineColor(4);
  double maxTime=GetMaxX(g);
  if(endpl<0.0001) endpl=maxTime;
  
  printf("--- Edge parameters ---\n");
  g->Fit(flow,"Q+","",edgeright,maxTime);
  g->Fit(fbas,"Q+","",0.,endpl);
  double basel=fbas->GetParameter(0);
  double baselsig=fbas->GetParError(0);
  double baselc,baselsigc;
  GetMeanAndRMSCounts(g,0.,endpl,baselc,baselsigc);
  printf("Baseline from fit = %f+-%f    from counts=%f rms=%f\n",basel,baselsig,baselc,baselsigc);
  double sign=flow->GetParameter(0);
  double signsig=flow->GetParError(0);
  double signc,signsigc;
  GetMeanAndRMSCounts(g,edgeright,25000,signc,signsigc);
  printf("Signal level from fit = %f+-%f    from counts=%f rms %f\n",sign,signsig,signc,signsigc);

  bool isSignal=kFALSE;
  if(edgeright>edgeleft && TMath::Abs(sign-basel)>0.001) isSignal=kTRUE;//sign-basel>1
  printf("IsSignal = %d  edgeleft=%f  edgeright=%f\n",isSignal,edgeleft,edgeright);
   double amplitude=-999.;
  double t10=-999.;
  double t50=-999.;
  double t90=-999.;
  double t55=-999.;
  double t60=-999.;
  double t65=-999.;
  double t70=-999.;
  double t75=-999.;
  double t80=-999.;
  double t85=-999.;
  if(isSignal){
   amplitude=TMath::Abs(sign-basel);
    t10=FindOnGraph(g,basel-0.1*amplitude,0.,maxTime,4,kTRUE);
    t50=FindOnGraph(g,basel-0.5*amplitude,0.,maxTime,4);
    t55=FindOnGraph(g,basel-0.55*amplitude,0.,maxTime,4);
    t60=FindOnGraph(g,basel-0.6*amplitude,0.,maxTime,4);
    t65=FindOnGraph(g,basel-0.65*amplitude,0.,maxTime,4);
    t70=FindOnGraph(g,basel-0.7*amplitude,0.,maxTime,4);
    t75=FindOnGraph(g,basel-0.75*amplitude,0.,maxTime,4);
    t80=FindOnGraph(g,basel-0.8*amplitude,0.,maxTime,4);
    t85=FindOnGraph(g,basel-0.85*amplitude,0.,maxTime,4);
    t90=FindOnGraph(g,basel-0.9*amplitude,0.,maxTime,4);
  
  }
  printf("t10 = %f  t50=%f  t55=%f t60=%f t65=%f t70=%f t75=%f t80=%f t85=%f t90=%f\n",t10,t50,t55,t60,t65,t70,t75,t80,t85,t90);
  
  params[0]=basel;
  params[1]=sign;
  params[2]=amplitude;
  params[3]=t10;
  params[4]=t50;
  params[5]=t55;
  params[6]=t60;
  params[7]=t65;
  params[8]=t70;
  params[9]=t75;
  params[10]=t80;
  params[11]=t85;
  params[12]=t90;
  params[13]=t90-t10;
}

void PlotWaveforms1(){
    TFile* inTree = new TFile("AO10P_B6_Vbb_4.0V_09092022_TTree.root");
    TFile* inTree1 = new TFile("AO10P_B6_v.1_Vbb_4V_19092022_TTree.root");
    TFile* inTree2 = new TFile("AO10P_B6_v.1_Vbb_4V_200922_prox_v1_TTree.root");

    TTree* tree=(TTree*)inTree->Get("treeParams");
    TTree* tree1=(TTree*)inTree1->Get("treeParams");
    TTree* tree2=(TTree*)inTree2->Get("treeParams");
    int ev;
    double fallTime;
    double ampl;
     int ev1;
    double fallTime1;
    double ampl1;
     int ev2;
    double fallTime2;
    double ampl2;
    tree->SetBranchAddress("Event",&ev);
    tree->SetBranchAddress("FallTimeCh1",&fallTime);
    tree->SetBranchAddress("SignalAmplCh1",&ampl);

    tree->SetBranchAddress("Event",&ev1);
    tree->SetBranchAddress("FallTimeCh1",&fallTime1);
    tree->SetBranchAddress("SignalAmplCh1",&ampl1);

    tree->SetBranchAddress("Event",&ev2);
    tree->SetBranchAddress("FallTimeCh1",&fallTime2);
    tree->SetBranchAddress("SignalAmplCh1",&ampl2);
    vector<TGraph*> graph_slow;
    vector<TGraph*> graph_fast;
    vector<TGraph*> graph_v1v1;
    vector<TGraph*> graph_v1v2;
    TMultiGraph *ftt=new TMultiGraph("Ftt","Fall Time vs range of time");
    TLegend* legLineb=new TLegend(0.60,0.13,0.80,0.33);
    legLineb->SetFillColor(kWhite);
    legLineb->SetFillStyle(1001);
    TLegendEntry* entVtxb;
    int h=0,d=0;
    for(int ient=0; ient<tree->GetEntriesFast(); ient++){
        tree->GetEvent(ient);
        if(h<30){
        if(fallTime < 1500 && fallTime > 600 ){//100
            TFile* infile = new TFile("AO10P_B6_Vbb_4.0V_09092022.root");
            TGraph* g = (TGraph*)infile->Get(Form("grEv%dChanC1samp25", ient));
            graph_slow.push_back(g);
            h++;
        }
    }
    
    if(d<30){
        if(fallTime < 600 && fallTime > 0){ //&& ampl > 1 ,300
            TFile* infile = new TFile("AO10P_B6_Vbb_4.0V_09092022.root");
            TGraph* g = (TGraph*)infile->Get(Form("grEv%dChanC1samp25", ient));
            graph_fast.push_back(g);
            d++;
        }
        }
    }
    d=0;
    h=0;
    for(int ient=0; ient<tree1->GetEntriesFast(); ient++)
    {
    if(d<30 && ampl1>10)
    {
      TFile* infile = new TFile("AO10P_B6_v.1_Vbb_4V_200922_prox_v1.root");
            TGraph* g = (TGraph*)infile->Get(Form("grEv%dChanC1samp25", ient));
            graph_v1v1.push_back(g);
            d++;
    }
    }
    for(int ient=0; ient<tree2->GetEntriesFast(); ient++)
    {
       if(h<30 && ampl2>10)
    {
      TFile* infile = new TFile("AO10P_B6_v.1_Vbb_4V_19092022.root");
            TGraph* g = (TGraph*)infile->Get(Form("grEv%dChanC1samp25", ient));
            graph_v1v2.push_back(g);
            h++;
    }
    }
    

    cout<<graph_slow.size()<<endl;

     TCanvas* c1 = new TCanvas("c1", "c1", 80, 80, 1500, 1000);
     TMultiGraph* mg_slow = new TMultiGraph("mg_slow", Form("Wafeforms of %d signals with Fall time between 600 ps and 1.5 ns", (int)graph_slow.size()));
     for(int i = 0; i < (int)graph_slow.size(); i++){
         graph_slow.at(i)->SetMarkerStyle(8);
         graph_slow.at(i)->SetLineWidth(2);
         graph_slow.at(i)->SetMarkerColor(15);//i+50
         graph_slow.at(i)->SetLineColor(15);//i+50
         mg_slow->Add(graph_slow.at(i));
     }
     mg_slow->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
     mg_slow->GetYaxis()->SetTitle("Signal [mV]");
     mg_slow->GetXaxis()->SetRangeUser(15000, 30000);
     mg_slow->Draw("AL");

     cout<<graph_fast.size()<<endl;

     TCanvas* c2 = new TCanvas("c2", "c2", 80, 80, 1500, 1000);
     TMultiGraph* mg_fast = new TMultiGraph("mg_fast", Form("Wafeforms of %d signals with Fall time below 600 ps", (int)graph_fast.size()));
     for(int i = 0; i < (int)graph_fast.size(); i++){
        graph_fast.at(i)->SetMarkerStyle(8);
         graph_fast.at(i)->SetLineWidth(2);
         graph_fast.at(i)->SetMarkerColor(1);//i+50
         graph_fast.at(i)->SetLineColor(1);
         mg_fast->Add(graph_fast.at(i));
     }
     mg_fast->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
     mg_fast->GetYaxis()->SetTitle("Signal [mV]");
     mg_fast->GetXaxis()->SetRangeUser(15000, 30000);
     mg_fast->Draw("AL");
 
    vector<TGraph*> smooth_slow;
    vector<TGraph*> smooth_fast;
    double params[14];

    TCanvas* c99=new TCanvas("c99","c99", 80, 80, 1500, 1000);
    
    //gPad->SetGrid();
    TMultiGraph* mg_v1v1 = new TMultiGraph("mg_v1v1", Form("Wafeforms of %d signals  of car V1 prox v1", (int)graph_v1v1.size()));
     for(int i = 0; i < (int)graph_v1v1.size(); i++){
        graph_v1v1.at(i)->SetMarkerStyle(8);
         graph_v1v1.at(i)->SetLineWidth(2);
         graph_v1v1.at(i)->SetMarkerColor(4);//i+50
         graph_v1v1.at(i)->SetLineColor(4);
         mg_v1v1->Add(graph_v1v1.at(i));
     }
     mg_v1v1->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
     mg_v1v1->GetYaxis()->SetTitle("Signal [mV]");
     mg_v1v1->GetXaxis()->SetRangeUser(15000, 30000);
     mg_v1v1->Draw("AL");

       TCanvas* c98=new TCanvas("c98","c98", 80, 80, 1500, 1000);
    
    //gPad->SetGrid();
    TMultiGraph* mg_v1v2 = new TMultiGraph("mg_v1v2", Form("Wafeforms of %d signals  of car V1 prox v2", (int)graph_v1v2.size()));
     for(int i = 0; i < (int)graph_v1v2.size(); i++){
        graph_v1v2.at(i)->SetMarkerStyle(8);
         graph_v1v2.at(i)->SetLineWidth(2);
         graph_v1v2.at(i)->SetMarkerColor(8);//i+50
         graph_v1v2.at(i)->SetLineColor(8);
         mg_v1v2->Add(graph_v1v1.at(i));
     }
     mg_v1v2->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
     mg_v1v2->GetYaxis()->SetTitle("Signal [mV]");
     mg_v1v2->GetXaxis()->SetRangeUser(15000, 30000);
     mg_v1v2->Draw("AL");

     /*TCanvas* c100=new TCanvas("c99","c99", 80, 80, 1500, 1000);
    c100->cd()
    gPad->SetGrid();
    TMultiGraph* mg_v1v11 = new TMultiGraph("mg_v1v1", Form("Wafeforms of %d signals  of car V1 prox v1", (int)graph_v1v1.size()));
     for(int i = 0; i < (int)graph_v1v1.size(); i++){
        graph_v1v1.at(i)->SetMarkerStyle(8);
         graph_v1v1.at(i)->SetLineWidth(2);
         graph_v1v1.at(i)->SetMarkerColor(4);//i+50
         graph_v1v1.at(i)->SetLineColor(4);
         mg_v1v11->Add(graph_v1v1.at(i));
     }
     mg_v1v1->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
     mg_v1v1->GetYaxis()->SetTitle("Signal [mV]");
     mg_v1v1->GetXaxis()->SetRangeUser(15000, 30000);
     mg_v1v1->Draw("AL");

     entVtxb=legLineb->AddEntry(graph_v1v1.at(0),"Prox v1 car v1");
     entVtxb->SetTextColor(4);
       //TCanvas* c98=new TCanvas("c99","c99", 80, 80, 1500, 1000);
    //c98->cd()
    //gPad->SetGrid();
    TMultiGraph* mg_v1v21 = new TMultiGraph("mg_v1v21", Form("Wafeforms of %d signals  of car V1 prox v2/v1 and car v2 prox v2", (int)graph_v1v2.size()));
     for(int i = 0; i < (int)graph_v1v2.size(); i++){
      if(i<30)
      {
        graph_v1v2.at(i)->SetMarkerStyle(8);
         graph_v1v2.at(i)->SetLineWidth(2);
         graph_v1v2.at(i)->SetMarkerColor(8);//i+50
         graph_v1v2.at(i)->SetLineColor(8);
         mg_v1v21->Add(graph_v1v1.at(i));
     }
     mg_v1v21->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
     mg_v1v21->GetYaxis()->SetTitle("Signal [mV]");
     mg_v1v21->GetXaxis()->SetRangeUser(15000, 30000);
     mg_v1v21->Draw("AL");
     }

     entVtxb=legLineb->AddEntry(graph_v1v2.at(0),"Prox v2 car v1");
     entVtxb->SetTextColor(8);

     TMultiGraph* mg_v1v21 = new TMultiGraph("mg_v1v2", Form("Wafeforms of %d signals  of car V2 v2", (int)graph_fast.size()));
     for(int i = 0; i < (int)graph_v1v2.size(); i++){
      if(i<30)
      {
         graph_v1v2.at(i)->SetMarkerStyle(8);
         graph_v1v2.at(i)->SetLineWidth(2);
         graph_v1v2.at(i)->SetMarkerColor(8);//i+50
         graph_v1v2.at(i)->SetLineColor(8);
         mg_v1v21->Add(graph_v1v1.at(i));
     }
     mg_v1v21->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
     mg_v1v21->GetYaxis()->SetTitle("Signal [mV]");
     mg_v1v21->GetXaxis()->SetRangeUser(15000, 30000);
     mg_v1v21->Draw("AL");
}

     entVtxb=legLineb->AddEntry(graph_v1v2.at(0),"Prox v2 car v1");
     entVtxb->SetTextColor(8);

*/

    TCanvas* c5 = new TCanvas("c5", "c5", 80, 80, 1500, 1000);
    c5->cd();
    gPad->SetGrid();
    TLine* lt10_slow;
    TLine* lt90_slow;
    TLine* lt50_slow;
    TLine* llow_slow;
    TLine* lbas_slow;
    TLine* l50_slow;
     printf("************ SLOW WAVEFORMS ****************\n");
    for(int i = 0; i < (int)graph_slow.size(); i++){
        smooth_slow.push_back(Smooth(graph_slow.at(i), 10));
        // printf("Analysing signal %d\n", i);
        if(i == 7){
            ProcessEvent(graph_slow.at(i),params,kTRUE);
            cout<<"l'ampiezza slow vale "<<params[2]<<endl;
            lt10_slow = new TLine(params[3], -20, params[3], 60);
            lt90_slow = new TLine(params[12], -20, params[12], 60);
            lt50_slow = new TLine(params[4], -20, params[4], 60);
            llow_slow = new TLine(params[12], params[1], 25000, params[1]);
            lbas_slow = new TLine(21500, params[0], 25000, params[0]);
            l50_slow = new TLine(21500, params[0]-0.5*params[2], 25000, params[0]-0.5*params[2]);
            lt10_slow->SetLineColor(kOrange+1);
            lt90_slow->SetLineColor(kOrange+1);
            lt50_slow->SetLineColor(kOrange+1);
            llow_slow->SetLineColor(kOrange+1);
            lbas_slow->SetLineColor(kOrange+1);
            l50_slow->SetLineColor(kOrange+1);
            lt10_slow->SetLineWidth(2);
            lt50_slow->SetLineWidth(2);
            lt90_slow->SetLineWidth(2);
            llow_slow->SetLineWidth(2);
            lbas_slow->SetLineWidth(2);
            l50_slow->SetLineWidth(2);
            smooth_slow.at(i)->SetMarkerStyle(8);
            smooth_slow.at(i)->SetLineWidth(2);
            smooth_slow.at(i)->SetMarkerColor(kOrange+4);
            smooth_slow.at(i)->SetLineColor(kOrange+4);
            //smooth_slow.at(i)->Draw("AL");
            graph_slow.at(i)->SetMarkerColor(kOrange+1);
            graph_slow.at(i)->SetLineColor(kOrange+1);
            graph_slow.at(i)->SetLineWidth(2);
            graph_slow.at(i)->GetXaxis()->SetTitle("Time [ps]");
            graph_slow.at(i)->GetYaxis()->SetTitle("Signal [mV]");
            graph_slow.at(i)->Draw("AL");
            lt10_slow->Draw("same");
            lt90_slow->Draw("same");
            lt50_slow->Draw("same");
            llow_slow->Draw("same");
            lbas_slow->Draw("same");
            l50_slow->Draw("same");
            TGraphErrors* ftr2= new TGraphErrors(9);
            ftr2->SetTitle("Fall Time vs range");
            for(int j=4;j<13;j++)
            {
              ftr2->SetPoint(j,5*j+30,(params[j]-params[3])*80/(j*5+20));
              ftr2->SetPointError(j,0.,0.);
            }
          
            ftr2->SetMarkerStyle(8);
            ftr2->SetMarkerColor(94);
            ftr2->SetLineColor(94);
            ftt->Add(ftr2);
            entVtxb=legLineb->AddEntry(ftr2,"Slow");
            entVtxb->SetTextColor(94);
        }
    }
    printf("\n");
    TLine* lt10_fast;
    TLine* lt90_fast;
    TLine* lt50_fast;
    TLine* llow_fast;
    TLine* lbas_fast;
    TLine* l50_fast;
     printf("************ FAST WAVEFORMS ****************\n");
    for(int i = 0; i < (int)graph_fast.size(); i++){
        smooth_fast.push_back(Smooth(graph_fast.at(i), 10));
         printf("Analysing signal %d\n", i);
         //ProcessEvent(graph_fast.at(i),params,kTRUE);
        if(i == 14){
          ProcessEvent(graph_fast.at(i),params,kTRUE);
          cout<<"l'ampiezza fast vale "<<params[2]<<endl;
            lt10_fast= new TLine(params[3], -20, params[3], 60);
            lt90_fast= new TLine(params[12], -20, params[12], 60);
            lt50_fast= new TLine(params[4], -20, params[4], 60);
            llow_fast= new TLine(params[12], params[1], 25000, params[1]);
            lbas_fast= new TLine(21500, params[0], 25000, params[0]);
            l50_fast= new TLine(21500, params[0]-0.5*params[2], 25000, params[0]-0.5*params[2]);
            lt10_fast->SetLineColor(kBlue);
            lt90_fast->SetLineColor(kBlue);
            lt50_fast->SetLineColor(kBlue);
            llow_fast->SetLineColor(kBlue);
            lbas_fast->SetLineColor(kBlue);
            l50_fast->SetLineColor(kBlue);
            lt10_fast->SetLineWidth(2);
            lt90_fast->SetLineWidth(2);
            lt50_fast->SetLineWidth(2);
            llow_fast->SetLineWidth(2);
            lbas_fast->SetLineWidth(2);
            smooth_fast.at(i)->SetMarkerStyle(8);
            smooth_fast.at(i)->SetLineWidth(2);
            smooth_fast.at(i)->SetMarkerColor(kBlue);
            smooth_fast.at(i)->SetLineColor(kBlue);
            //smooth_fast.at(i)->Draw("Lsame");
            graph_fast.at(i)->SetMarkerColor(kBlue);
            graph_fast.at(i)->SetLineColor(kBlue);
            graph_fast.at(i)->SetLineWidth(2);
            graph_fast.at(i)->Draw("Lsame");
            lt10_fast->Draw("same");
            lt90_fast->Draw("same");
            lt50_fast->Draw("same");
            llow_fast->Draw("same");
            lbas_fast->Draw("same");
            l50_fast->Draw("same");
          
          TGraphErrors* ftr= new TGraphErrors(9);
          ftr->SetTitle("Fall Time vs range");
          for(int j=4;j<13;j++)
          {
            ftr->SetPoint(j,5*j+30,(params[j]-params[3])*80/(j*5+20));
            ftr->SetPointError(j,0.,0.);
            cout<<"(params[j]-params[3])*80/(j*5+20)"<<(params[j]-params[3])*80/(j*5+20)<<endl;
          }
          
          ftr->SetMarkerStyle(8);
          ftr->SetMarkerColor(4);
          ftr->SetLineColor(4);
          ftt->Add(ftr);
          entVtxb=legLineb->AddEntry(ftr,"Fast");
          entVtxb->SetTextColor(4);
       
          




            // ProcessEvent(graph_fast.at(i),params,kTRUE);
            // lt10_fast = new TLine(params[3], -20, params[3], 60);
            // lt90_fast = new TLine(params[5], -20, params[5], 60);
            // llow_fast = new TLine(params[3], params[1], 25000, params[1]);
            // lbas_fast = new TLine(21500, params[0], 23000, params[0]);
            // lt10_fast->SetLineColor(kBlue);
            // lt90_fast->SetLineColor(kBlue);
            // llow_fast->SetLineColor(kBlue);
            // lbas_fast->SetLineColor(kBlue);
            // lt10_fast->SetLineWidth(2);
            // lt90_fast->SetLineWidth(2);
            // llow_fast->SetLineWidth(2);
            // lbas_fast->SetLineWidth(2);
            // smooth_fast.at(i)->SetMarkerStyle(8);
            // smooth_fast.at(i)->SetLineWidth(2);
            // smooth_fast.at(i)->SetMarkerColor(kBlue);
            // smooth_fast.at(i)->SetLineColor(kBlue);
            // graph_fast.at(i)->SetMarkerColor(i+50);
            // graph_fast.at(i)->SetLineColor(i+50);
            // graph_fast.at(i)->SetLineWidth(2);
            // smooth_fast.at(i)->Draw("Lsame");
            // graph_fast.at(i)->Draw("Lsame");
            // lt10_fast->Draw("same");
            // lt90_fast->Draw("same");
            // llow_fast->Draw("same");
            // lbas_fast->Draw("same");
        }
    }
     TCanvas* cff=new TCanvas("cff","falltimes",80,80,1500,1000);
          cff->cd();
          gPad->SetGrid();
          ftt->GetXaxis()->SetNdivisions(1020);
          ftt->GetXaxis()->SetTitle("%");
          ftt->GetYaxis()->SetTitle("Fall Time [ps]");
          ftt->GetXaxis()->SetRangeUser(40.,100.);
          ftt->GetYaxis()->SetRangeUser(0.,1000.);
          ftt->Draw("AP");
          legLineb->Draw();

     
          

    TCanvas* c3 = new TCanvas("c3", "c3", 80, 80, 1500, 1000);
    TMultiGraph* mg_smooth_slow = new TMultiGraph("mg_smooth_slow", Form("Smoothed Wafeforms of %d signals with Fall time between 600 ps and 1.5 ns", (int)graph_slow.size()));
    for(int i = 0; i < (int)smooth_slow.size(); i++){
        smooth_slow.at(i)->SetMarkerStyle(8);
        smooth_slow.at(i)->SetLineWidth(2);
        smooth_slow.at(i)->SetMarkerColor(i+50);
        smooth_slow.at(i)->SetLineColor(i+50);
        mg_smooth_slow->Add(smooth_slow.at(i));
    }
    mg_smooth_slow->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
    mg_smooth_slow->GetYaxis()->SetTitle("Signal [mV]");
    mg_smooth_slow->GetXaxis()->SetRangeUser(15000, 30000);
    mg_smooth_slow->Draw("AL");

    TCanvas* c4 = new TCanvas("c4", "c4", 80, 80, 1500, 1000);
    TMultiGraph* mg_smooth_fast = new TMultiGraph("mg_smooth_fast", Form("Smoothed Wafeforms of %d signals with Fall time below 600 ps", (int)graph_fast.size()));
    for(int i = 0; i < (int)smooth_fast.size(); i++){
        smooth_fast.at(i)->SetMarkerStyle(8);
        smooth_fast.at(i)->SetLineWidth(2);
        smooth_fast.at(i)->SetMarkerColor(i+50);
        smooth_fast.at(i)->SetLineColor(i+50);
        mg_smooth_fast->Add(smooth_fast.at(i));
    }
    mg_smooth_fast->GetXaxis()->SetTitle("Time sample (a semple every 25 ps)");
    mg_smooth_fast->GetYaxis()->SetTitle("Signal [mV]");
    mg_smooth_fast->GetXaxis()->SetRangeUser(15000, 30000);
    mg_smooth_fast->Draw("AL");

}