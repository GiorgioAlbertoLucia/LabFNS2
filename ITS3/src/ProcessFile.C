#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>

/*
Script to Prossess 55Fe acquisitions selecting signals greater than 1 mV:
  > SCOPE:
    - Baseline and Amplitude measured as defined by Francesco Prino and calculated on the raw signal
    - Fall Time defined as 2*(t50-t10)
  > ADC:
    - Baseline measured as average on the first 100 samples
    - Amplitude measured as the difference between the Baseline and the 101th sample
    - Fall time set manually to -999.
*/

// function to smoth the signal to find the falling edge
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

// Get max x point on the TGraph
double GetMaxX(TGraph* g){
  double xmax=0.;
  for(int j=0; j<g->GetN(); j++){
    double x,y;
    g->GetPoint(j,x,y);
    if(x>xmax) xmax=x;
  }
  return xmax;
}

// Compute the signal derivative to find the time in which the signal starts
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

// Creates the TGraph of the signal derivative
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

// Create the TGraph of the negative derivative of the signal
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

// Compute average and RMS
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

// Find a specific point on a TGraph
double FindOnGraph(TGraph* gcount, double y, double xmin, double xmax, int interpolate, bool backw=kFALSE){
  int jfirst=0;
  int dstep=1;
  if(backw){
    jfirst=gcount->GetN();
    dstep=-1;
  }
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

// Find edges on TGraph
void FindEdge(TGraph* gcount,TGraph* gnegder, TGraph* gder, double& endplateau, double& edgeleft, double& edgeright){
  
  // first very rough: compute flat levels on the left and on the right and check their difference
  double maxTime=GetMaxX(gcount);
  double levleft,rmsleft,levright,rmsright; 
  // printf("maxTime=%f\n",maxTime);
  GetMeanAndRMSCounts(gcount,0.,2000.,levleft,rmsleft);
  GetMeanAndRMSCounts(gcount,maxTime-2000,maxTime,levright,rmsright);
  double y50=0.5*(levleft+levright);
  // printf("y50=%f\n",y50);
  double t50fromleft=FindOnGraph(gcount,y50,0.,maxTime,4);
  double t50fromright=FindOnGraph(gcount,y50,0.,maxTime,4,kTRUE);
  double roughsig=levleft-levright;
  // printf("roughsig=%f\n",roughsig);
  // printf("Rough signal = %f Rough edge position = %f %f\n",roughsig,t50fromleft,t50fromright);
  double minSearchWindow=0;
  double maxSearchWindow=maxTime;
  if(roughsig>0.0005){
    minSearchWindow=TMath::Min(t50fromleft,t50fromright)-6000.;
    if(minSearchWindow<0) minSearchWindow=0;
    maxSearchWindow=TMath::Max(t50fromleft,t50fromright)+6000.;
    if(maxSearchWindow>maxTime) maxSearchWindow=maxTime;
  }
  // printf("Search window = %f %f\n",minSearchWindow,maxSearchWindow);
  
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
    // printf("Maximum adjacent points with negative derivative: t_maxn=%f   n_neg=%f\n",xmaxn,cmaxn);
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
  // printf("Minimum of derivative: xminder=%f   dermin=%f\n",xminder,dermin);
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
  // printf("Region of negative derivative: xleft=%f   xright=%f\n",xleft,xright);
  if(xmaxn>0 && TMath::Abs(xmaxn-xleft)<5000 && TMath::Abs(xmaxn-xright)<5000){
    if(xleft>xmaxn) xleft=xmaxn;
    if(xright<xmaxn) xright=xmaxn;
  }
  // printf("Edge range after analysis of derivative: xleft=%f   xright=%f\n",xleft,xright);

  // Fourth step: start from left and seach for N points with couns < baseline-3sigma
  double cmean,crms;
  GetMeanAndRMSCounts(gcount,0.,xleft,cmean,crms);
  // printf("Mean before edge = %f rms = %f\n",cmean,crms);
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
  // printf("Left Edge from baseline-N*rms = %f\n",xleft2);
  if(xleft2>0){
    endplateau=TMath::Min(xleft,xleft2);
    edgeleft=TMath::Max(xleft,xleft2);
    edgeright=xright;
    // printf("Edge range after all steps: endplateau=%f   edgeleft=%f   edgeright=%f\n",endplateau,edgeleft,edgeright);
    return;
  }
  endplateau=0;
  edgeleft=0;
  edgeright=0;
  return;
}

// Function to compute the baseline, amplitude and fall time for pixels read by the Scope
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
  // printf("--- Edge finding on fine graph ---\n");
  FindEdge(g,gnegd,gd,endpl,edgeleft,edgeright);
  TF1* fbas=new TF1("fbas","[0]");
  TF1* flow=new TF1("flow","[0]");
  fbas->SetLineColor(2);
  flow->SetLineColor(4);
  double maxTime=GetMaxX(g);
  if(endpl<0.0001) endpl=maxTime;
  
  // printf("--- Edge parameters ---\n");
  g->Fit(flow,"Q+","",edgeright,maxTime);
  g->Fit(fbas,"Q+","",0.,endpl);
  double basel=fbas->GetParameter(0);
  double baselsig=fbas->GetParError(0);
  double baselc,baselsigc;
  GetMeanAndRMSCounts(g,0.,endpl,baselc,baselsigc);
  // printf("Baseline from fit = %f+-%f    from counts=%f rms=%f\n",basel,baselsig,baselc,baselsigc);
  double sign=flow->GetParameter(0);
  double signsig=flow->GetParError(0);
  double signc,signsigc;
  GetMeanAndRMSCounts(g,edgeright,25000,signc,signsigc);
  // printf("Signal level from fit = %f+-%f    from counts=%f rms %f\n",sign,signsig,signc,signsigc);

  bool isSignal=kFALSE;
  if(edgeright>edgeleft && TMath::Abs(sign-basel)>1) isSignal=kTRUE;//sign-basel>1
  // printf("IsSignal = %d  edgeleft=%f  edgeright=%f\n",isSignal,edgeleft,edgeright);
  double amplitude=-999.;
  double t10=-999.;
  double t50=-999.;
  double t90=-999.;

  if(isSignal){
    amplitude=TMath::Abs(sign-basel);
    t10=FindOnGraph(g,basel-0.1*amplitude,0.,maxTime,4,kTRUE);
    t50=FindOnGraph(g,basel-0.5*amplitude,0.,maxTime,4);
    t90=FindOnGraph(g,basel-0.9*amplitude,0.,maxTime,4);
  }
  
  params[0]=basel;
  params[1]=sign;
  params[2]=amplitude;
  params[3]=t10;
  params[4]=t50;
  params[5]=t90;
  params[6]=2*(t50-t10);

  if(!plot){
    delete gs;
    delete gnegd;
    delete gsnegd;
    delete gd;
    delete fbas;
    delete flow;
    return;
  }
  
  TCanvas* c1= new TCanvas("c1","",1600,800);
  c1->Divide(3,2);
  c1->cd(1);
  gPad->SetRightMargin(0.04);
  g->Draw();
  flow->Draw("same");
  fbas->Draw("same");
  TLine* lp=new TLine(endpl,g->GetYaxis()->GetXmin(),endpl,g->GetYaxis()->GetXmax());
  lp->SetLineColor(2);
  lp->SetLineStyle(2);
  lp->Draw();
  TLine* ll=new TLine(edgeleft,g->GetYaxis()->GetXmin(),edgeleft,g->GetYaxis()->GetXmax());
  ll->SetLineColor(kGreen+1);
  ll->SetLineStyle(9);
  ll->Draw();
  TLine* lr=new TLine(edgeright,g->GetYaxis()->GetXmin(),edgeright,g->GetYaxis()->GetXmax());
  lr->SetLineColor(kGreen+1);
  lr->SetLineStyle(9);
  lr->Draw();
  TLine* l10=new TLine(t10,g->GetYaxis()->GetXmin(),t10,g->GetYaxis()->GetXmax());
  l10->SetLineColor(kMagenta+1);
  l10->SetLineStyle(7);
  l10->Draw();
  TLine* l90=new TLine(t90,g->GetYaxis()->GetXmin(),t90,g->GetYaxis()->GetXmax());
  l90->SetLineColor(kMagenta+1);
  l90->SetLineStyle(7);
  l90->Draw();
  TLatex* t1=new TLatex(edgeright+1400,basel,Form("Baseline = %.2f mV\n",basel*1000.));
  t1->SetTextColor(2);
  t1->Draw();
  TLatex* t2=new TLatex(edgeright+1400,basel-0.01*basel,Form("Signal ampl = %.2f mV\n",amplitude*1000.));
  t2->SetTextColor(4);
  t2->Draw();
  TLatex* t3=new TLatex(edgeright+1400,basel-0.02*basel,Form("Fall time = %.1f ps\n",t90-t10));
  t3->SetTextColor(kMagenta+1);
  t3->Draw();
  c1->cd(2);
  gPad->SetRightMargin(0.04);
  gs->Draw();
  c1->cd(4);
  gPad->SetRightMargin(0.04);
  gnegd->Draw();
  c1->cd(5);
  gPad->SetRightMargin(0.04);
  gd->Draw();
  
}

// Function to calculate amplitude, baseline and fall time on pixels read by the ADC 
void ProcessADCEvent(TGraph* g, double params[20]){
  TF1* fbas=new TF1("fbas","[0]");
  fbas->SetLineColor(2);
  double maxTime=250*99;
  
  // printf("--- Edge parameters ---\n");
  g->Fit(fbas,"Q+","",0,maxTime);
  double basel=fbas->GetParameter(0);
  vector<double> amp;

  for(int i = 0; i < 200; i++){
    amp.push_back(g->GetPointY(i));
  }
  double sign = *std::min_element(amp.begin(), amp.end());

  double amplitude=-999.;
  double t10=-999.;
  double t50=-999.;
  double t90=-999.;

  if(basel-sign > 1){
    amplitude = basel-sign;
  }
  else{
    sign = -999;
  }
  
  // printf("bl = %f ul=%f amp=%f\n",basel,sign,amplitude);
  
  params[0]=basel;
  params[1]=sign;
  params[2]=amplitude;
  params[3]=t10;
  params[4]=t50;
  params[5]=t90;
  params[6]=-999;
}

// Main funtion of the file: manage the read and write part
void ProcessFile(TString filnam, int maxEv=999999){

  const int nParsPerChan=7;
  TString varNames[nParsPerChan]={"Baseline","MinLevel","SignalAmpl","t10","t50","t90","FallTime"};  //"t55","t60","t65","t70","t75","t80","t85",,"BaselCoarse","RecovTimeExpo","RecovTimeLin"
  double paramschan[7];
  double params[16*nParsPerChan];
  int ev;
  long timest;

  TString outfilnam=filnam.Data();
  outfilnam.ReplaceAll(".root","_TTree.root");
  TFile* outFile=new TFile(outfilnam.Data(),"recreate");
  TTree* outTree=new TTree("treeParams","tree of parameters");
  outTree->Branch("Event",&ev,"Event/I");
  outTree->Branch("Timestamp",&timest,"Timestamp/l");
  for(int ipx=0; ipx<16; ipx++){
    for(int ivar=0; ivar<nParsPerChan; ivar++){
      outTree->Branch(Form("%sPx%d",varNames[ivar].Data(),ipx),&params[ivar+ipx*nParsPerChan],Form("%sPx%d/D",varNames[ivar].Data(),ipx));
    }
  }

  TFile* f=new TFile(filnam.Data());
  int totEv=0;
  int nkeys=f->GetNkeys();
  TList* lkeys=f->GetListOfKeys();
  int period;
  int lastev=0;
  int px;
  // printf("----- Number of keys in file = %d\n",nkeys);
  TKey* k=(TKey*)lkeys->At(nkeys-1);
  TString cname=k->GetClassName();
  TString oname=k->GetName();
  if(cname=="TGraph"){
    sscanf(oname.Data(),"grEv%dPx%dsamp%d",&ev,&px,&period);
    if(ev>lastev){ 
      lastev=ev;
    }
  }  
  // printf("----- Number of events = %d\n",lastev);
  if(lastev>maxEv) lastev=maxEv;

  int dum;
  char ddum[2];
  for(int iev=0; iev<=lastev; iev++){
    if(iev%100==0) printf("----- Processing Event %d -----\n",iev);
    for(int ipx=0; ipx<16; ipx++){
      // printf("--- Pixel %d ---\n",ipx);
      if(ipx==5 || ipx==6 || ipx==9 || ipx==10){
        TGraph* g=(TGraph*)f->Get(Form("grEv%dPx%dsamp25",iev,ipx));
        const char* grtit=g->GetTitle();
        sscanf(grtit,"Event %d Pixel %2s Time %ld",&dum,ddum,&timest);
        ev=iev;
        ProcessEvent(g,paramschan,kFALSE);
        for(int ivar=0; ivar<nParsPerChan; ivar++) params[ivar+ipx*nParsPerChan]=paramschan[ivar];
        delete g;
      }
      else{
        TGraph* glong=(TGraph*)f->Get(Form("grEv%dPx%dsamp250000",iev,ipx));
        const char* grtit=glong->GetTitle();
        sscanf(grtit,"Event %d Pixel %2s Time %ld",&dum,ddum,&timest);
        ev=iev;
        ProcessADCEvent(glong,paramschan);
        for(int ivar=0; ivar<nParsPerChan; ivar++) params[ivar+ipx*nParsPerChan]=paramschan[ivar];
        delete glong;
      }
    }
    outTree->Fill();
  }
  
  outFile->cd();
  outTree->Write();
  outFile->Close();
  
}
