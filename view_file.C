#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH2F.h>
#include <TTree.h>
#include <TApplication.h>
#include <TKey.h>

using namespace std;

void loopdir(TDirectoryFile* f1, string indent = "") {
    TIter next(f1->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        cout << indent << key->GetName() << endl;
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TDirectoryFile")) continue;
        loopdir(dynamic_cast<TDirectoryFile*>(key->ReadObj()), indent + key->GetName() + "/");
    }
}
