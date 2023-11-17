 #include <TFile.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TF1.h>
#include <TProfile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <fstream>

const int maxFiles=5;
TH1F* hFallTimeCh1[maxFiles];
TH1F* hFallTimeCh2[maxFiles];
TH1F* hFallTimeCh3[maxFiles];
TH1F* hFallTimeCh4[maxFiles];
TH1F* hFallTime10_50Ch1[maxFiles];
TH1F* hFallTime10_55Ch1[maxFiles];
TH1F* hFallTime10_60Ch1[maxFiles];
TH1F* hFallTime10_65Ch1[maxFiles];
TH1F* hFallTime10_70Ch1[maxFiles];
TH1F* hFallTime10_75Ch1[maxFiles];
TH1F* hFallTime10_80Ch1[maxFiles];
TH1F* hFallTime10_85Ch1[maxFiles];
TH1F* hFallTime10_90Ch1[maxFiles];

TH1F* hFallTime10_50Ch1f[maxFiles];
TH1F* hFallTime10_55Ch1f[maxFiles];
TH1F* hFallTime10_60Ch1f[maxFiles];
TH1F* hFallTime10_65Ch1f[maxFiles];
TH1F* hFallTime10_70Ch1f[maxFiles];
TH1F* hFallTime10_75Ch1f[maxFiles];
TH1F* hFallTime10_80Ch1f[maxFiles];
TH1F* hFallTime10_85Ch1f[maxFiles];

TH1F* hFallTime10_50Ch2[maxFiles];
TH1F* hFallTime10_50Ch3[maxFiles];
TH1F* hFallTime10_50Ch4[maxFiles];

TH1F* hFallTimeClu1Ch1[maxFiles];
TH1F* hFallTimeClu1Ch2[maxFiles];
TH1F* hFallTimeClu1Ch3[maxFiles];
TH1F* hFallTimeClu1Ch4[maxFiles];

TH1F* hFallTimeCluxCh1[maxFiles];//x = >1
TH1F* hFallTimeCluxCh2[maxFiles];
TH1F* hFallTimeCluxCh3[maxFiles];
TH1F* hFallTimeCluxCh4[maxFiles];

TH1F* hAmplCh1[maxFiles];
TH1F* hAmplCh2[maxFiles];
TH1F* hAmplCh3[maxFiles];
TH1F* hAmplCh4[maxFiles];
TH1F* hAmplClu1Ch1[maxFiles];
TH1F* hAmplClu1Ch2[maxFiles];
TH1F* hAmplClu1Ch3[maxFiles];
TH1F* hAmplClu1Ch4[maxFiles];
TH1F* hBlCh1[maxFiles];
TH1F* hBlCh2[maxFiles];
TH1F* hBlCh3[maxFiles];
TH1F* hBlCh4[maxFiles];
TH1F* hBlClu1Ch1[maxFiles];
TH1F* hBlClu1Ch2[maxFiles];
TH1F* hBlClu1Ch3[maxFiles];
TH1F* hBlClu1Ch4[maxFiles];
TH1F* hFttot[maxFiles];
TH1F* hAmptot[maxFiles];
TH2F* hFtAmpCh1[maxFiles];
TH2F* hFtAmpCh2[maxFiles];
TH2F* hFtAmpCh3[maxFiles];
TH2F* hFtAmpCh4[maxFiles];
TH2F* hFtAmpClu1Ch1[maxFiles];
TH2F* hFtAmpClu1Ch2[maxFiles];
TH2F* hFtAmpClu1Ch3[maxFiles];
TH2F* hFtAmpClu1Ch4[maxFiles];
TH2F* hFtAmpCluxCh1[maxFiles];
TH2F* hFtAmpCluxCh2[maxFiles];
TH2F* hFtAmpCluxCh3[maxFiles];
TH2F* hFtAmpCluxCh4[maxFiles];
TH2F* hFtAmpClutot[maxFiles];

void FormatHistos(TObjArray* arrHisto, int theColor){
  int nHist=arrHisto->GetEntries();
  for(int jh=0; jh<nHist; jh++){
    TString cname=((TObject*)arrHisto->At(jh))->ClassName();
    if(cname.Contains("TH1")){
      TH1* h=(TH1*)arrHisto->At(jh);
      h->SetLineWidth(3);
      h->SetLineColor(theColor);
      h->GetYaxis()->SetTitleOffset(1.25);
    }else if(cname.Contains("TH2")){
      TH2* h=(TH2*)arrHisto->At(jh);
      h->SetLineColor(theColor);
    }    
  }
}

void WriteHistos(TObjArray* arrHisto, TString filnam){
  TFile* outf=new TFile(filnam.Data(),"recreate");
  int nHist=arrHisto->GetEntries();
  for(int jh=0; jh<nHist; jh++){
    TString cname=((TObject*)arrHisto->At(jh))->ClassName();
    if(cname.Contains("TH1")){
      TH1* h=(TH1*)arrHisto->At(jh);
      h->Write();
    }else if(cname.Contains("TH2")){
      TH2* h=(TH2*)arrHisto->At(jh);
      h->Write();
    }
  }
  outf->Close();
  delete outf;
}

void NormalizeHistos(TObjArray* arrHisto){
 
  int nHist=arrHisto->GetEntries();
  printf("nHist=%d\n",nHist);
  for(int jh=0; jh<nHist; jh++){
    TString cname=((TObject*)arrHisto->At(jh))->ClassName();
    if(cname.Contains("TH1")){
      TH1* h=(TH1*)arrHisto->At(jh);
      TString hname=h->GetName();
      if(hname.Contains("CluSiz") || hname.Contains("CluTyp") || hname.Contains("hMaxSigPix")) continue;
      double tot=h->Integral();
      printf("tot=%f\n",tot);
      if(tot>0){
	h->Scale(1./tot);
	h->GetYaxis()->SetTitle("Entries (a.u.)");
  
      }
    }
  }
}
// void SetHistoMaximum(TH1* hcur, TH1* href, double scal=1.05){
//   if(hcur->GetMaximum() > href->GetMaximum()) href->SetMaximum(hcur->GetMaximum()*scal);
// }

// void SetProfMaximum(TProfile* pcur, TProfile* pref, double scal1, double scal2){
//   if(pcur->GetMaximum() > pref->GetMaximum()*scal1) pref->SetMaximum(pcur->GetMaximum()*scal2);
// }
// void SetMaxima(int nfils){

//   for(int j=0; j<nfils; j++){
//     SetHistoMaximum(hFallTimeCh1[j],hFallTimeCh1[0]);
//     SetHistoMaximum(hFallTimeCh2[j],hFallTimeCh2[0]);
//     SetHistoMaximum(hFallTimeCh3[j],hFallTimeCh3[0]);
//     SetHistoMaximum(hFallTimeCh4[j],hFallTimeCh4[0]);
//     SetHistoMaximum(hAmplCh1[j],hAmplCh1[0]);
//     SetHistoMaximum(hAmplCh2[j],hAmplCh2[0]);
//     SetHistoMaximum(hAmplCh3[j],hAmplCh3[0]);
//     SetHistoMaximum(hAmplCh4[j],hAmplCh4[0]);
//     }

// }
void FillHistosFromTree(TFile* f, int jfil, int iTrigChan,string dfile="DerivativeCalTot.txt" ){
  printf("Fill histos with iTrigChan = %d\n",iTrigChan);
   double gain[3][4];
  double j5,j6,j9,j10;
  int aa=0;
    ifstream inputd(dfile.c_str());
    if(!inputd){
        cout<<"il file "<<dfile<<" non esiste"<<endl;
    }
  while(inputd>>j5>>j6>>j9>>j10)
  {

  gain[aa][0]=j5;
  gain[aa][1]=j6;
  gain[aa][2]=j9;
  gain[aa][3]=j10;
  aa++;

  }
cout<<"aa vale"<<aa<<endl;
  TTree* tree=(TTree*)f->Get("treeParams");
  int ev;
    for(int yy=0;yy<3;yy++)
    {
      for(int gg=0;gg<4;gg++)
      {
cout<<gain[yy][gg]<<" ";
      }
      cout<<endl;
    }
  
  double t10Vec[4];
  double t50Vec[4];
  double t55Vec[4];
  double t60Vec[4];
  double t65Vec[4];
  double t70Vec[4];
  double t75Vec[4];
  double t80Vec[4];
  double t85Vec[4];
  double t90Vec[4];
  double amplVec[4];
  double baselVec[4];
  double recoTimeVec[4];
  double fallTimeVec[4];
  tree->SetBranchAddress("Event",&ev);
  for(int k=0; k<4; k++){
    tree->SetBranchAddress(Form("t10Ch%d",k+1),&t10Vec[k]);
    tree->SetBranchAddress(Form("t50Ch%d",k+1),&t50Vec[k]);
    tree->SetBranchAddress(Form("t55Ch%d",k+1),&t55Vec[k]);
    tree->SetBranchAddress(Form("t60Ch%d",k+1),&t60Vec[k]);
    tree->SetBranchAddress(Form("t65Ch%d",k+1),&t65Vec[k]);
    tree->SetBranchAddress(Form("t70Ch%d",k+1),&t70Vec[k]);
    tree->SetBranchAddress(Form("t75Ch%d",k+1),&t75Vec[k]);
    tree->SetBranchAddress(Form("t80Ch%d",k+1),&t80Vec[k]);
    tree->SetBranchAddress(Form("t85Ch%d",k+1),&t85Vec[k]);
    tree->SetBranchAddress(Form("t90Ch%d",k+1),&t90Vec[k]);
    tree->SetBranchAddress(Form("FallTimeCh%d",k+1),&fallTimeVec[k]);
    tree->SetBranchAddress(Form("SignalAmplCh%d",k+1),&amplVec[k]);
    tree->SetBranchAddress(Form("BaselineCh%d",k+1),&baselVec[k]);
  }

  for(int ient=0; ient<tree->GetEntriesFast(); ient++){
    tree->GetEvent(ient);
    //if(amplVec[0]>=0 && fallTimeVec[0]>=0){
      if(amplVec[0]>=8. && fallTimeVec[0]>=0){ //&& fallTimeVec[0]<=600.){
	    hFallTimeCh1[jfil]->Fill(fallTimeVec[0]/1000.);
      hFttot[jfil]->Fill((t50Vec[0]/1000-t10Vec[0]/1000)*2);
      hFallTime10_50Ch1[jfil]->Fill((t50Vec[0]/1000-t10Vec[0]/1000)*2);
      hFallTime10_55Ch1[jfil]->Fill((t55Vec[0]/1000-t10Vec[0]/1000)*80/45);
      hFallTime10_60Ch1[jfil]->Fill((t60Vec[0]/1000-t10Vec[0]/1000)*80/50);
      hFallTime10_65Ch1[jfil]->Fill((t65Vec[0]/1000-t10Vec[0]/1000)*80/55);
      hFallTime10_70Ch1[jfil]->Fill((t70Vec[0]/1000-t10Vec[0]/1000)*80/60);
      hFallTime10_75Ch1[jfil]->Fill((t75Vec[0]/1000-t10Vec[0]/1000)*80/65);
      hFallTime10_80Ch1[jfil]->Fill((t80Vec[0]/1000-t10Vec[0]/1000)*80/70);
      hFallTime10_85Ch1[jfil]->Fill((t85Vec[0]/1000-t10Vec[0]/1000)*80/75);
     
      if(fallTimeVec[0]<600.){
        hFallTime10_50Ch1f[jfil]->Fill((t50Vec[0]/1000-t10Vec[0]/1000)*2);
        hFallTime10_55Ch1f[jfil]->Fill((t55Vec[0]/1000-t10Vec[0]/1000)*80/45);
        hFallTime10_60Ch1f[jfil]->Fill((t60Vec[0]/1000-t10Vec[0]/1000)*80/50);
        hFallTime10_65Ch1f[jfil]->Fill((t65Vec[0]/1000-t10Vec[0]/1000)*80/55);
        hFallTime10_70Ch1f[jfil]->Fill((t70Vec[0]/1000-t10Vec[0]/1000)*80/60);
        hFallTime10_75Ch1f[jfil]->Fill((t75Vec[0]/1000-t10Vec[0]/1000)*80/65);
        hFallTime10_80Ch1f[jfil]->Fill((t80Vec[0]/1000-t10Vec[0]/1000)*80/70);
        hFallTime10_85Ch1f[jfil]->Fill((t85Vec[0]/1000-t10Vec[0]/1000)*80/75);
      }
      
      hAmptot[jfil]->Fill(amplVec[0]/gain[jfil][0]);
	    hAmplCh1[jfil]->Fill(amplVec[0]);
      hBlCh1[jfil]->Fill(baselVec[0]);
      hFtAmpCh1[jfil]->Fill(amplVec[0],(t50Vec[0]/1000-t10Vec[0]/1000)*2);
      hFtAmpClutot[jfil]->Fill(amplVec[0],(t50Vec[0]/1000-t10Vec[0]/1000)*2);
     // passati a Ft 10-50
      //if(fallTimeVec[0]>0.85) printf("l'evento %d ha Fall Time alto\n",ient);
    }
    if(amplVec[1]>=8. && fallTimeVec[1]>=0){ // && fallTimeVec[1]<=600.){
      hAmptot[jfil]->Fill(amplVec[1]/gain[jfil][1]);
	    hFallTimeCh2[jfil]->Fill(fallTimeVec[1]/1000.);
      hFttot[jfil]->Fill((t50Vec[1]/1000-t10Vec[1]/1000)*2);
      hFallTime10_50Ch2[jfil]->Fill((t50Vec[1]/1000-t10Vec[1]/1000)*2);
	    hAmplCh2[jfil]->Fill(amplVec[1]);
      hBlCh2[jfil]->Fill(baselVec[1]);
      hFtAmpCh2[jfil]->Fill(amplVec[1],(t50Vec[1]/1000-t10Vec[1]/1000)*2);
       hFtAmpClutot[jfil]->Fill(amplVec[1],(t50Vec[1]/1000-t10Vec[1]/1000)*2);
    }     
    if(amplVec[2]>=8. && fallTimeVec[2]>=0) {// {&& fallTimeVec[2]<=600.){
      hAmptot[jfil]->Fill(amplVec[2]/gain[jfil][2]);
	    hFallTimeCh3[jfil]->Fill(fallTimeVec[2]/1000.);
      hFallTime10_50Ch3[jfil]->Fill((t50Vec[2]/1000-t10Vec[2]/1000)*2);
      hFttot[jfil]->Fill((t50Vec[2]/1000-t10Vec[2]/1000)*2);
	    hAmplCh3[jfil]->Fill(amplVec[2]);
      hBlCh3[jfil]->Fill(baselVec[2]);
      hFtAmpCh3[jfil]->Fill(amplVec[2],(t50Vec[2]/1000-t10Vec[2]/1000)*2);
      hFtAmpClutot[jfil]->Fill(amplVec[2],(t50Vec[1]/1000-t10Vec[2]/1000)*2);
    }     
    if(amplVec[3]>=8. && fallTimeVec[3]>=0) { //&& fallTimeVec[3]<=600.){
      hAmptot[jfil]->Fill(amplVec[3]/gain[jfil][3]);
	    hFallTimeCh4[jfil]->Fill(fallTimeVec[3]/1000.);
      hFttot[jfil]->Fill((t50Vec[3]/1000-t10Vec[3]/1000)*2);
      hFallTime10_50Ch4[jfil]->Fill((t50Vec[3]/1000-t10Vec[3]/1000)*2);
	    hAmplCh4[jfil]->Fill(amplVec[3]);
      hBlCh4[jfil]->Fill(baselVec[3]);
      hFtAmpCh4[jfil]->Fill(amplVec[3],(t50Vec[3]/1000-t10Vec[3]/1000)*2);
      hFtAmpClutot[jfil]->Fill(amplVec[2],(t50Vec[3]/1000-t10Vec[3]/1000)*2);
    }
    int clusiz=0;
    int clutyp=0;
    double totampl=0;
    int maxsigpix=-1;
    double maxsig=-999.;
    for(int k=0; k<4; k++){
      if(amplVec[k]>0){
  	clutyp+=1<<k;
  	clusiz+=1;
  	totampl+=amplVec[k];
	
      }
    }
    if(clusiz==1)
    {
    if(amplVec[0]>=8. && fallTimeVec[0]>=0){ //&& fallTimeVec[0]<=600.){
	    hFallTimeClu1Ch1[jfil]->Fill(fallTimeVec[0]/1000.);
	    hAmplClu1Ch1[jfil]->Fill(amplVec[0]);
      hBlClu1Ch1[jfil]->Fill(baselVec[0]);
      hFtAmpClu1Ch1[jfil]->Fill(amplVec[0],fallTimeVec[0]/1000.);
    }
    if(amplVec[1]>=8. && fallTimeVec[1]>=0){ // && fallTimeVec[1]<=600.){
	    hFallTimeClu1Ch2[jfil]->Fill(fallTimeVec[1]/1000.);
	    hAmplClu1Ch2[jfil]->Fill(amplVec[1]);
      hBlClu1Ch2[jfil]->Fill(baselVec[1]);
      hFtAmpClu1Ch2[jfil]->Fill(amplVec[1],fallTimeVec[1]/1000.);
    }     
    if(amplVec[2]>=8. && fallTimeVec[2]>=0) {// {&& fallTimeVec[2]<=600.){
	    hFallTimeClu1Ch3[jfil]->Fill(fallTimeVec[2]/1000.);
	    hAmplClu1Ch3[jfil]->Fill(amplVec[2]);
      hBlClu1Ch3[jfil]->Fill(baselVec[2]);
      hFtAmpClu1Ch3[jfil]->Fill(amplVec[2],fallTimeVec[2]/1000.);
    }     
    if(amplVec[3]>=8. && fallTimeVec[3]>=0) { //&& fallTimeVec[3]<=600.){
	    hFallTimeClu1Ch4[jfil]->Fill(fallTimeVec[3]/1000.);
	    hAmplClu1Ch4[jfil]->Fill(amplVec[3]);
      hBlClu1Ch4[jfil]->Fill(baselVec[3]);
      hFtAmpClu1Ch4[jfil]->Fill(amplVec[3],fallTimeVec[3]/1000.);
    }
    }

    if(clusiz>1)
    {
    if(amplVec[0]>=8. && fallTimeVec[0]>=0){ //&& fallTimeVec[0]<=600.){
	    hFallTimeCluxCh1[jfil]->Fill(fallTimeVec[0]/1000.);
	    
    }
    if(amplVec[1]>=8. && fallTimeVec[1]>=0){ // && fallTimeVec[1]<=600.){
	    hFallTimeCluxCh2[jfil]->Fill(fallTimeVec[1]/1000.);
	    
    }     
    if(amplVec[2]>=8. && fallTimeVec[2]>=0) {// {&& fallTimeVec[2]<=600.){
	    hFallTimeCluxCh3[jfil]->Fill(fallTimeVec[2]/1000.);
	    
    }     
    if(amplVec[3]>=8. && fallTimeVec[3]>=0) { //&& fallTimeVec[3]<=600.){
	    hFallTimeCluxCh4[jfil]->Fill(fallTimeVec[3]/1000.);
	    
    }
    }
  }
  
}

void PlotFromTree1(string dfile="DerivativeCalTot.txt",TString configFile="configuration.txt", bool normalizeToArea=kTRUE){

  int nFiles=0;
  TString fileNames[maxFiles];//maxfile è 5, dichiarato prima
  TString trigChan[maxFiles];
  int jTrigChan[maxFiles];
  int cols[maxFiles]={kMagenta+1,1};
  TString legEntry[maxFiles];
 

    
  FILE* cFile=fopen(configFile.Data(),"r");
  char suff[100];
  char line[200];
  fgets(line,200,cFile);
  sscanf(line,"%d %s",&nFiles,suff);
  printf("nFiles vale=%d \n",nFiles);
  if(nFiles>maxFiles){
    printf("ERROR: maximum number of files is %d\n",maxFiles);
    return;
  }
  TString suffix=suff;
  int readFiles=0;
  
  for(int jf=0; jf<nFiles; jf++){
    fgets(line,200,cFile);
    TString theLine(line);
    TObjArray* arrEnt=theLine.Tokenize(";");
    int nEnt=arrEnt->GetEntries();
    if(nEnt!=4){
      printf("ERROR: expect filename ; trigchan ; color ; legendtext\n");
      break;
    }
    for(int k=0; k<nEnt; k++){
      TObjString* str=(TObjString*)arrEnt->At(k);
      TString theStr=str->GetString();
      theStr.ReplaceAll("\n","");
      if(k==0) fileNames[jf]=theStr.Data();
      if(k==1) trigChan[jf]=theStr.Data();
      else if(k==2) cols[jf]=theStr.Atoi();
      else if(k==3) legEntry[jf]=theStr.Data();
      jTrigChan[jf]=0;
      if(trigChan[jf].Contains("J5")) jTrigChan[jf]=1;
      else if(trigChan[jf].Contains("J10")) jTrigChan[jf]=3;
      else if(trigChan[jf].Contains("OR")) jTrigChan[jf]=0;
    }
    //    arrEnt->Delete();
    delete arrEnt;
    readFiles++;
    if(feof(cFile)) break;
  }
  fclose(cFile);
  if(readFiles!=nFiles){
    printf("ERROR: mismatch between number of expected files (%d) and number of read lines (%d)\n",nFiles,readFiles);
    return;
  }
  printf("Number of files to be analyzed = %d suffix for plots = %s\n",nFiles,suffix.Data());
  for(int jf=0; jf<nFiles; jf++){
    printf("File %d = %s  trigger channel = %s(%d)  Color = %d  legend Entry = %s\n",jf, fileNames[jf].Data(),trigChan[jf].Data(),jTrigChan[jf],cols[jf],legEntry[jf].Data());
  }
  if(nFiles==0) return;
 
  TObjArray* arrHisto = new TObjArray();
  double cnt04[4];
  double cntall[4];
  double cntFT04[4];
  double cntFTall[4];
  double cntMaxFT04[4];
  double cntMax[4];
  for(int j=0; j<nFiles; j++){//fa istogrammi
    hFallTimeClu1Ch1[j]=new TH1F(Form("hFallTimeClu1Ch1_%d",j)," Fall Time_CluSize1 ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTimeClu1Ch2[j]=new TH1F(Form("hFallTimeClu1Ch2_%d",j)," Fall Time_CluSize1 ; Fall Time Ch2 (ns) ; Entries",100,0.,10.);
    hFallTimeClu1Ch3[j]=new TH1F(Form("hFallTimeClu1Ch3_%d",j)," Fall Time_CluSize1 ; Fall Time Ch3 (ns) ; Entries",100,0.,10.);
    hFallTimeClu1Ch4[j]=new TH1F(Form("hFallTimeClu1Ch4_%d",j)," Fall Time_CluSize1 ; Fall Time Ch4 (ns) ; Entries",100,0.,10.);
    hFttot[j]=new TH1F(Form("hFttot_%d",j)," Fall Time ; Fall Time  (ns) ; Entries",200,0.,10.);
    hFallTimeCluxCh1[j]=new TH1F(Form("hFallTimeCluxCh1_%d",j)," Fall Time_CluSizex ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTimeCluxCh2[j]=new TH1F(Form("hFallTimeCluxCh2_%d",j)," Fall Time_CluSizex ; Fall Time Ch2 (ns) ; Entries",100,0.,10.);
    hFallTimeCluxCh3[j]=new TH1F(Form("hFallTimeCluxCh3_%d",j)," Fall Time_CluSizex ; Fall Time Ch3 (ns) ; Entries",100,0.,10.);
    hFallTimeCluxCh4[j]=new TH1F(Form("hFallTimeCluxCh4_%d",j)," Fall Time_CluSizex ; Fall Time Ch4 (ns) ; Entries",100,0.,10.);
    hAmptot[j]=new TH1F(Form("hAmptot_%d",j)," Amplitude all chan equalized ;Amplitude (mV); Entries",150,0.,150.);

    hAmplClu1Ch1[j]=new TH1F(Form("hAmplClu1Ch1_%d",j)," ; Signal Amplitude CluSize1 Ch1 (mV) ; Entries",100.,0.,100.);
    hAmplClu1Ch2[j]=new TH1F(Form("hAmplClu1Ch2_%d",j)," ; Signal Amplitude CluSize1 Ch2 (mV) ; Entries",100.,0.,100.);
    hAmplClu1Ch3[j]=new TH1F(Form("hAmplClu1Ch3_%d",j)," ; Signal Amplitude CluSize1 Ch3 (mV) ; Entries",100.,0.,100.);
    hAmplClu1Ch4[j]=new TH1F(Form("hAmplClu1Ch4_%d",j)," ; Signal Amplitude CluSize1 Ch4 (mV) ; Entries",100.,0.,100.);
    hBlClu1Ch1[j]=new TH1F(Form("hBlClu1Ch1_%d",j)," ; Signal Baseline Ch1 (mV) CluSize1 ; Entries",100,30.,70.);
    hBlClu1Ch2[j]=new TH1F(Form("hBlClu1Ch2_%d",j)," ; Signal baseline Ch2 (mV) CluSize1 ; Entries",100,30.,70.);
    hBlClu1Ch3[j]=new TH1F(Form("hBlClu1Ch3_%d",j)," ; Signal Baseline Ch3 (mV) CluSize1 ; Entries",100,30.,70.);
    hBlClu1Ch4[j]=new TH1F(Form("hBlClu1Ch4_%d",j)," ; Signal Baseline Ch4 (mV) CluSize1 ; Entries",100,30.,70.);
    hFtAmpClu1Ch1[j]=new TH2F(Form("hFtAmpClu1Ch1_%d",j)," ; Signal Amplitude CluSize1  (mV) ; Fall Time  (ns) ; Entries",100.,0.,100.,100,0.,10.);
    hFtAmpClu1Ch2[j]=new TH2F(Form("hFtAmpClu1Ch2_%d",j)," ; Signal Amplitude CluSize1  (mV) ; Fall Time  (ns) ; Entries",100.,0.,100.,100,0.,10.);
    hFtAmpClu1Ch3[j]=new TH2F(Form("hFtAmpClu1Ch3_%d",j)," ; Signal Amplitude CluSize1  (mV) ; Fall Time  (ns) ; Entries",100.,0.,100.,100,0.,10.);
    hFtAmpClu1Ch4[j]=new TH2F(Form("hFtAmpClu1Ch4_%d",j)," ; Signal Amplitude CluSize1  (mV) ; Fall Time  (ns) ; Entries",100.,0.,100.,100,0.,10.);

    hFallTimeCh1[j]=new TH1F(Form("hFallTimeCh1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTimeCh2[j]=new TH1F(Form("hFallTimeCh2_%d",j)," All clusters ; Fall Time Ch2 (ns) ; Entries",100,0.,10.);
    hFallTimeCh3[j]=new TH1F(Form("hFallTimeCh3_%d",j)," All clusters ; Fall Time Ch3 (ns) ; Entries",100,0.,10.);
    hFallTimeCh4[j]=new TH1F(Form("hFallTimeCh4_%d",j)," All clusters ; Fall Time Ch4 (ns) ; Entries",100,0.,10.);
    hFallTime10_50Ch1[j]=new TH1F(Form("hFallTime (50-10)*2_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTime10_50Ch2[j]=new TH1F(Form("hFallTime (50-10)*2_Ch2__%d",j)," All clusters ; Fall Time Ch2 (ns) ; Entries",100,0.,10.);
    hFallTime10_50Ch3[j]=new TH1F(Form("hFallTime (50-10)*2_Ch3__%d",j)," All clusters ; Fall Time Ch3 (ns) ; Entries",100,0.,10.);
    hFallTime10_50Ch4[j]=new TH1F(Form("hFallTime (50-10)*2_Ch4_%d",j)," All clusters ; Fall Time Ch4 (ns) ; Entries",100,0.,10.);

    hFallTime10_50Ch1f[j]=new TH1F(Form("hFallTime 50-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);
    hFallTime10_55Ch1[j]=new TH1F(Form("hFallTime 55-10_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTime10_60Ch1[j]=new TH1F(Form("hFallTime 60-10_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTime10_65Ch1[j]=new TH1F(Form("hFallTime 65-10_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTime10_70Ch1[j]=new TH1F(Form("hFallTime 70-10_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTime10_75Ch1[j]=new TH1F(Form("hFallTime 75-10_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTime10_80Ch1[j]=new TH1F(Form("hFallTime 80-10_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTime10_85Ch1[j]=new TH1F(Form("hFallTime 85-10_Ch1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);

    hFallTime10_55Ch1f[j]=new TH1F(Form("hFallTime 55-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);
    hFallTime10_60Ch1f[j]=new TH1F(Form("hFallTime 60-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);
    hFallTime10_65Ch1f[j]=new TH1F(Form("hFallTime 65-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);
    hFallTime10_70Ch1f[j]=new TH1F(Form("hFallTime 70-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);
    hFallTime10_75Ch1f[j]=new TH1F(Form("hFallTime 75-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);
    hFallTime10_80Ch1f[j]=new TH1F(Form("hFallTime 80-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);
    hFallTime10_85Ch1f[j]=new TH1F(Form("hFallTime 85-10_Ch1f_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",25,0.,0.6);

    hAmplCh1[j]=new TH1F(Form("hAmplCh1_%d",j)," ; Signal Amplitude Ch1 (mV) ; Entries",100.,0.,100.);
    hAmplCh2[j]=new TH1F(Form("hAmplCh2_%d",j)," ; Signal Amplitude Ch2 (mV) ; Entries",100.,0.,100.);
    hAmplCh3[j]=new TH1F(Form("hAmplCh3_%d",j)," ; Signal Amplitude Ch3 (mV) ; Entries",100.,0.,100.);
    hAmplCh4[j]=new TH1F(Form("hAmplCh4_%d",j)," ; Signal Amplitude Ch4 (mV) ; Entries",100.,0.,100.);
    hBlCh1[j]=new TH1F(Form("hBlCh1_%d",j)," ; Signal Baseline Ch1 (mV) ; Entries",100,30.,70.);
    hBlCh2[j]=new TH1F(Form("hBlCh2_%d",j)," ; Signal baseline Ch2 (mV) ; Entries",100,30.,70.);
    hBlCh3[j]=new TH1F(Form("hBlCh3_%d",j)," ; Signal Baseline Ch3 (mV) ; Entries",100,30.,70.);
    hBlCh4[j]=new TH1F(Form("hBlCh4_%d",j)," ; Signal Baseline Ch4 (mV) ; Entries",100,30.,70.);
    hFtAmpCh1[j]=new TH2F(Form("hFtAmpCh1_%d",j)," ; Signal Amplitude  (mV) ; Fall Time  (ns) ; Entries",100,0.,100.,100,0.,10.);
    hFtAmpCh2[j]=new TH2F(Form("hFtAmpCh2_%d",j)," ; Signal Amplitude  (mV) ; Fall Time  (ns) ; Entries",100,0.,100.,100,0.,10.);
    hFtAmpCh3[j]=new TH2F(Form("hFtAmpCh3_%d",j)," ; Signal Amplitude  (mV) ; Fall Time  (ns) ; Entries",100,0.,100.,100,0.,10.);
    hFtAmpCh4[j]=new TH2F(Form("hFtAmpCh4_%d",j)," ; Signal Amplitude  (mV) ; Fall Time  (ns) ; Entries",100,0.,100.,100,0.,10.);
    hFtAmpClutot[j]=new TH2F(Form("hFtAmpAllCh_%d",j)," ; Signal Amplitude  (mV) ; Fall Time  (ns) ; Entries",100,0.,100.,100,0.,10.);
  
  

    int indexh=0;
    arrHisto->AddAtAndExpand(hFallTime10_50Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_55Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_60Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_65Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_70Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_75Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_80Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_85Ch1[j],indexh++);

    arrHisto->AddAtAndExpand(hFallTime10_50Ch1f[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_55Ch1f[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_60Ch1f[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_65Ch1f[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_70Ch1f[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_75Ch1f[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_80Ch1f[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_85Ch1f[j],indexh++);

    arrHisto->AddAtAndExpand(hFallTime10_50Ch2[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_50Ch3[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTime10_50Ch4[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh4[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh4[j],indexh++);
    arrHisto->AddAtAndExpand(hBlCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hBlCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hBlCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hBlCh4[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpCh4[j],indexh++);


    arrHisto->AddAtAndExpand(hFallTimeClu1Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeClu1Ch2[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeClu1Ch3[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeClu1Ch4[j],indexh++);

     arrHisto->AddAtAndExpand(hFallTimeCluxCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCluxCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCluxCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCluxCh4[j],indexh++);

    arrHisto->AddAtAndExpand(hAmplClu1Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplClu1Ch2[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplClu1Ch3[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplClu1Ch4[j],indexh++);
    arrHisto->AddAtAndExpand(hBlClu1Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hBlClu1Ch2[j],indexh++);
    arrHisto->AddAtAndExpand(hBlClu1Ch3[j],indexh++);
    arrHisto->AddAtAndExpand(hBlClu1Ch4[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpClu1Ch1[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpClu1Ch2[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpClu1Ch3[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpClu1Ch4[j],indexh++);
    arrHisto->AddAtAndExpand(hFtAmpClutot[j],indexh++);
    arrHisto->AddAtAndExpand(hFttot[j],indexh++);
    arrHisto->AddAtAndExpand(hAmptot[j],indexh++);
   

    printf("indexh=%d\n",indexh);
    int nHist1=arrHisto->GetEntries();
  printf("nHist1=%d\n",nHist1);

    
    TFile* f=new TFile(fileNames[j].Data());//ti mette tutti i grafici in un Tfile
    FillHistosFromTree(f,j,jTrigChan[j]-1,dfile);
    TString outHisFilNam=fileNames[j];
    outHisFilNam.ReplaceAll("_TTree.root","_Histos.root");
    WriteHistos(arrHisto,outHisFilNam);
    fileNames[j].ReplaceAll("_TTree.root","");
    FormatHistos(arrHisto,cols[j]);
    //arrHisto->Clear();
    hFtAmpCh1[j]->SetStats(0);
    hFtAmpCh2[j]->SetStats(0);
    hFtAmpCh3[j]->SetStats(0);
    hFtAmpCh4[j]->SetStats(0);
    hFtAmpClutot[j]->SetStats(0);
hFtAmpClutot[j]->GetZaxis()->SetRangeUser(0,600.);
    
    hFtAmpClu1Ch1[j]->SetStats(0);
    hFtAmpClu1Ch2[j]->SetStats(0);
    hFtAmpClu1Ch3[j]->SetStats(0);
    hFtAmpClu1Ch4[j]->SetStats(0);
    
      printf("PRIMA \n");
    if(normalizeToArea) NormalizeHistos(arrHisto);
    
    printf("DOPO \n");
  }
  // SetMaxima(nFiles);

  TCanvas* cA = new TCanvas("cA","Amplitudes",1650,900);
  cA->Divide(2,2);
  cA->cd(1);
  TLegend* leg= new TLegend(0.12,0.6,0.5,0.89);
  leg->SetMargin(0.1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh1[j]->Draw("histo");
    else hAmplCh1[j]->Draw("histo,sames");
    leg->AddEntry(hAmplCh1[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  cA->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh2[j]->Draw("histo");
    else hAmplCh2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh2[j]->GetLineColor());
    gPad->Modified();
  }
  cA->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh3[j]->Draw("histo");
    else hAmplCh3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh3[j]->GetLineColor());
    gPad->Modified();
  } 
  cA->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh4[j]->Draw("histo");
    else hAmplCh4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh4[j]->GetLineColor());
    gPad->Modified();
  }
  //cA->SaveAs(Form("Amplitudes_%s_cut600ps.png",suffix.Data()));

  TCanvas* cAc = new TCanvas("cAc","Amplitudes CluSize_1",1650,900);
  cAc->Divide(2,2);
  cAc->cd(1);
  
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplClu1Ch1[j]->Draw("histo");
    else hAmplClu1Ch1[j]->Draw("histo,sames");
    leg->AddEntry(hAmplClu1Ch1[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplClu1Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplClu1Ch1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  cAc->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplClu1Ch2[j]->Draw("histo");
    else hAmplClu1Ch2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplClu1Ch2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplClu1Ch2[j]->GetLineColor());
    gPad->Modified();
  }
    cAc->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplClu1Ch3[j]->Draw("histo");
    else hAmplClu1Ch3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplClu1Ch3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplClu1Ch3[j]->GetLineColor());
    gPad->Modified();
  }
    cAc->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplClu1Ch4[j]->Draw("histo");
    else hAmplClu1Ch4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplClu1Ch4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplClu1Ch4[j]->GetLineColor());
    gPad->Modified();
  }
  cAc->SaveAs(Form("Amplitudes_%s_ClSi_1.png",suffix.Data()));

  TCanvas* cAx = new TCanvas("cAx","Fall Time CluSize>1",1650,900);
  cAx->Divide(2,2);
  cAx->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCluxCh1[j]->Draw("histo");
    else hFallTimeCluxCh1[j]->Draw("histo,sames");
    leg->AddEntry(hFallTimeCluxCh1[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCluxCh1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCluxCh1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();

    cAx->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCluxCh2[j]->Draw("histo");
    else hFallTimeCluxCh2[j]->Draw("histo,sames");
    leg->AddEntry(hFallTimeCluxCh2[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCluxCh2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCluxCh2[j]->GetLineColor());
    gPad->Modified();
  }

   cAx->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCluxCh3[j]->Draw("histo");
    else hFallTimeCluxCh3[j]->Draw("histo,sames");
    leg->AddEntry(hFallTimeCluxCh3[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCluxCh3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCluxCh3[j]->GetLineColor());
    gPad->Modified();
  }

   cAx->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCluxCh4[j]->Draw("histo");
    else hFallTimeCluxCh4[j]->Draw("histo,sames");
    leg->AddEntry(hFallTimeCluxCh4[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCluxCh4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCluxCh4[j]->GetLineColor());
    gPad->Modified();
  }

  TCanvas* cF = new TCanvas("cF","FallTimes",1650,900);
  cF->Divide(2,2);
  cF->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh1[j]->Draw("histo");
    else hFallTimeCh1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  cF->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh2[j]->Draw("histo");
    else hFallTimeCh2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh2[j]->GetLineColor());
    gPad->Modified();
  }
  cF->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh3[j]->Draw("histo");
    else hFallTimeCh3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh3[j]->GetLineColor());
    gPad->Modified();
  } 
  cF->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh4[j]->Draw("histo");
    else hFallTimeCh4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh4[j]->GetLineColor());
    gPad->Modified();
  }
  
  cF->SaveAs(Form("FallTime_%s_Allclu.png",suffix.Data()));

  TCanvas* cFc = new TCanvas("cFc","FallTimes clusize1",1650,900);
  cFc->Divide(2,2);
  cFc->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeClu1Ch1[j]->Draw("histo");
    else hFallTimeClu1Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeClu1Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeClu1Ch1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
   cFc->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeClu1Ch2[j]->Draw("histo");
    else hFallTimeClu1Ch2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeClu1Ch2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeClu1Ch2[j]->GetLineColor());
    gPad->Modified();
  }
  TCanvas* cF2t = new TCanvas("cF2t","FallTimes",1650,900);
 
  cF2t->cd();
  for(int j=0; j<nFiles; j++){
    if(j==0) hFttot[j]->Draw("histo");
    else hFttot[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFttot[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFttot[j]->GetLineColor());
    gPad->Modified();
  }
  TCanvas* cF2 = new TCanvas("cF2","FallTimes (50-10)*2",1650,900);
  cF2->Divide(2,2);
  cF2->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_50Ch1[j]->Draw("histo");
    else hFallTime10_50Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_50Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_50Ch1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
   cF2->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_50Ch2[j]->Draw("histo");
    else hFallTime10_50Ch2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_50Ch2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_50Ch2[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
    cF2->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_50Ch3[j]->Draw("histo");
    else hFallTime10_50Ch3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_50Ch3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_50Ch3[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
    cF2->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_50Ch4[j]->Draw("histo");
    else hFallTime10_50Ch4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_50Ch4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_50Ch4[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  cFc->SaveAs(Form("FallTime10_50_%s.png",suffix.Data()));
    cFc->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeClu1Ch3[j]->Draw("histo");
    else hFallTimeClu1Ch3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeClu1Ch3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeClu1Ch3[j]->GetLineColor());
    gPad->Modified();
  }
      cFc->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeClu1Ch4[j]->Draw("histo");
    else hFallTimeClu1Ch4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeClu1Ch4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeClu1Ch4[j]->GetLineColor());
    gPad->Modified();
  }
  cFc->SaveAs(Form("FallTime_%s_clusize1.png",suffix.Data()));

  TCanvas* cB = new TCanvas("cB","Baselines",1650,900);
  cB->Divide(2,2);
  cB->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hBlCh1[j]->Draw("histo");
    else hBlCh1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hBlCh1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hBlCh1[j]->GetLineColor());
    gPad->Modified();
  }

  cB->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hBlCh2[j]->Draw("histo");
    else hBlCh2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hBlCh2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hBlCh2[j]->GetLineColor());
    gPad->Modified();
  }

  cB->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hBlCh3[j]->Draw("histo");
    else hBlCh3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hBlCh3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hBlCh3[j]->GetLineColor());
    gPad->Modified();
  }

   cB->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hBlCh4[j]->Draw("histo");
    else hBlCh4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hBlCh4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hBlCh4[j]->GetLineColor());
    gPad->Modified();
  }
  cB->SaveAs(Form("Baseline_%s_allClu_bin_ps.png",suffix.Data()));

  TCanvas* cfa = new TCanvas("cfa","Histos Fall Time vs Amplitude",1650,900);
  cfa->Divide(2,2);
  cfa->cd(1);
  gPad->SetGrid();//Logz();
  hFtAmpCh1[0]->Draw("colz");
  for(int j=1; j<nFiles; j++) hFtAmpCh1[j]->Draw("same,colz");
  cfa->cd(2);
  gPad->SetGrid();//Logz();
  hFtAmpCh2[0]->Draw("colz");
  for(int j=1; j<nFiles; j++) hFtAmpCh2[j]->Draw("same,colz");
  cfa->cd(3);
  gPad->SetGrid();//Logz();
  hFtAmpCh3[0]->Draw("colz");
  for(int j=1; j<nFiles; j++) hFtAmpCh3[j]->Draw("same,colz");
  cfa->cd(4);
  gPad->SetGrid();//Logz();
  hFtAmpCh4[0]->Draw("colz");
  for(int j=1; j<nFiles; j++) hFtAmpCh4[j]->Draw("same,colz");
  leg->Draw();
  cfa->SaveAs(Form("FtAmpl_%s_AllClu_bin_ps.png",suffix.Data()));

   TCanvas* cfa2 = new TCanvas("cfa2","Histos Fall Time vs Amplitude",1650,900);
  cfa2->cd();
  
  gPad->SetGrid();//Logz();
  hFtAmpClutot[0]->Draw("colz");
  for(int j=1; j<nFiles; j++) hFtAmpClutot[j]->Draw("same,colz");

   TCanvas* caa2 = new TCanvas("caa2","Amplitue equalized",1650,900);
  caa2->cd();
  
  gPad->SetGrid();//Logz();
  
  for(int j=0; j<nFiles; j++) 
  {
    if(j==0) hAmptot[j]->Draw("histo");
    else hAmptot[j]->Draw("histo,sames");
  }
  leg->Draw();

  TCanvas* cFt = new TCanvas("cFt","FallTimes differente range",1650,900);
  cFt->Divide(4,2);
  cFt->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_50Ch1[j]->Draw("histo");
    else hFallTime10_50Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_50Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_50Ch1[j]->GetLineColor());
    gPad->Modified();
  }

  
  cFt->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_55Ch1[j]->Draw("histo");
    else hFallTime10_55Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_55Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_55Ch1[j]->GetLineColor());
    gPad->Modified();
  }

  
  
    cFt->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_60Ch1[j]->Draw("histo");
    else hFallTime10_60Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_60Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_60Ch1[j]->GetLineColor());
    gPad->Modified();
  }

   
  
    cFt->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_65Ch1[j]->Draw("histo");
    else hFallTime10_65Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_65Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_65Ch1[j]->GetLineColor());
    gPad->Modified();
  }

   
  
    cFt->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_70Ch1[j]->Draw("histo");
    else hFallTime10_70Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_70Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_70Ch1[j]->GetLineColor());
    gPad->Modified();
  }

   
  
    cFt->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_75Ch1[j]->Draw("histo");
    else hFallTime10_75Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_75Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_75Ch1[j]->GetLineColor());
    gPad->Modified();
  }

  
      cFt->cd(7);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_80Ch1[j]->Draw("histo");
    else hFallTime10_80Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_80Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_80Ch1[j]->GetLineColor());
    gPad->Modified();
  }

   

      cFt->cd(8);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_85Ch1[j]->Draw("histo");
    else hFallTime10_85Ch1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_85Ch1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_85Ch1[j]->GetLineColor());
    gPad->Modified();
  }

  
  //cFt->SaveAs(Form("FallTime_different_range_%s_.png",suffix.Data()));
  TCanvas* cFtf = new TCanvas("cFtf","FallTimes differente range fast",1650,900);
  cFtf->Divide(4,2);
  cFtf->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_50Ch1f[j]->Draw("histo");
    else hFallTime10_50Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_50Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_50Ch1f[j]->GetLineColor());
    gPad->Modified();
  }
  cFtf->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_55Ch1f[j]->Draw("histo");
    else hFallTime10_55Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_55Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_55Ch1f[j]->GetLineColor());
    gPad->Modified();
  }
  cFtf->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_60Ch1f[j]->Draw("histo");
    else hFallTime10_60Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_60Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_60Ch1f[j]->GetLineColor());
    gPad->Modified();
  }
    cFtf->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_65Ch1f[j]->Draw("histo");
    else hFallTime10_65Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_65Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_65Ch1f[j]->GetLineColor());
    gPad->Modified();
  }
      cFtf->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_70Ch1f[j]->Draw("histo");
    else hFallTime10_70Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_70Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_70Ch1f[j]->GetLineColor());
    gPad->Modified();
  }
         cFtf->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_75Ch1f[j]->Draw("histo");
    else hFallTime10_75Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_75Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_75Ch1f[j]->GetLineColor());
    gPad->Modified();
  }
  

          cFtf->cd(7);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_80Ch1f[j]->Draw("histo");
    else hFallTime10_80Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_80Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_80Ch1f[j]->GetLineColor());
    gPad->Modified();
  
  }
         cFtf->cd(8);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime10_85Ch1f[j]->Draw("histo");
    else hFallTime10_85Ch1f[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime10_85Ch1f[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime10_85Ch1f[j]->GetLineColor());
    gPad->Modified();
  }
  
  
  





  
}