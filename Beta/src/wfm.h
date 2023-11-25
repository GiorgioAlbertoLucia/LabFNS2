//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 20 23:00:11 2023 by ROOT version 6.26/10
// from TTree wfm/
// found on file: Beta/data/output/BetaMeas_Lab2_clean.root
//////////////////////////////////////////////////////////

#ifndef wfm_h
#define wfm_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class wfm {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           nt2;
   Double_t        t2[1002];   //[nt2]
   Int_t           nt3;
   Double_t        t3[1002];   //[nt3]
   Int_t           nw2;
   Double_t        w2[1002];   //[nw2]
   Int_t           nw3;
   Double_t        w3[1002];   //[nw3]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_nt2;   //!
   TBranch        *b_t2;   //!
   TBranch        *b_nt3;   //!
   TBranch        *b_t3;   //!
   TBranch        *b_nw2;   //!
   TBranch        *b_w2;   //!
   TBranch        *b_nw3;   //!
   TBranch        *b_w3;   //!

   wfm(TTree *tree=0);
   virtual ~wfm();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   virtual void    BuildTree(const char * filename);
};

#endif

#ifdef wfm_cxx
wfm::wfm(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Beta/data/output/BetaMeas_Lab2_clean.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Beta/data/output/BetaMeas_Lab2_clean.root");
      }
      f->GetObject("wfm",tree);

   }
   Init(tree);
}

wfm::~wfm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t wfm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t wfm::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void wfm::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nt2", &nt2, &b_nt2);
   fChain->SetBranchAddress("t2", t2, &b_t2);
   fChain->SetBranchAddress("nt3", &nt3, &b_nt3);
   fChain->SetBranchAddress("t3", t3, &b_t3);
   fChain->SetBranchAddress("nw2", &nw2, &b_nw2);
   fChain->SetBranchAddress("w2", w2, &b_w2);
   fChain->SetBranchAddress("nw3", &nw3, &b_nw3);
   fChain->SetBranchAddress("w3", w3, &b_w3);
   Notify();
}

Bool_t wfm::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void wfm::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef wfm_cxx
