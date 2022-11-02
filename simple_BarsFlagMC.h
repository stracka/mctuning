#ifndef simple_BarsFlagMC_h
#define simple_BarsFlagMC_h

#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <bitset>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class simple_BarsFlagMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        xo;
   Double_t        yo;
   Double_t        zo;
   Double_t        pxo;
   Double_t        pyo;
   Double_t        pzo;
   Int_t           event;
   Double_t        t_event;
   Bool_t          pbar;
   Bool_t          mc;
   Int_t           nTPCHits;
   Int_t           nTPCMCHits;
   Int_t           nDigi;
   Int_t           nBars;
   Int_t           nCompleteBars;
   Int_t           nBarEnds;
   vector<int>     *bars_id;
   vector<int>     *bars_ntrks;
   vector<int>     *bars_flag;
   vector<double>  *bars_edep;
   vector<double>  *bars_path;
   vector<double>  *bars_z;
   vector<double>  *bars_t;
   vector<double>  *bars_phi;
   vector<double>  *bars_atop;
   vector<double>  *bars_abot;
   vector<double>  *bars_ttop;
   vector<double>  *bars_tbot;
   vector<double>  *bars_ttop_twc;
   vector<double>  *bars_tbot_twc;
   vector<double>  *pairs_tof;
   vector<double>  *pairs_dphi;
   vector<double>  *pairs_dzeta;
   vector<double>  *pairs_dist;
   Double_t        tof_min;
   Double_t        tof_max;
   Double_t        tof_mean;
   Double_t        tof_std;
   Double_t        dphi_min;
   Double_t        dphi_max;
   Double_t        dphi_mean;
   Double_t        dphi_std;
   Double_t        dzeta_min;
   Double_t        dzeta_max;
   Double_t        dzeta_mean;
   Double_t        dzeta_std;
   Double_t        dist_min;
   Double_t        dist_max;
   Double_t        dist_mean;
   Double_t        dist_std;

   // List of branches
   TBranch        *b_xo;   //!
   TBranch        *b_yo;   //!
   TBranch        *b_zo;   //!
   TBranch        *b_pxo;   //!
   TBranch        *b_pyo;   //!
   TBranch        *b_pzo;   //!
   TBranch        *b_event;   //!
   TBranch        *b_t_event;   //!
   TBranch        *b_pbar;   //!
   TBranch        *b_mc;   //!
   TBranch        *b_nTPCHits;   //!
   TBranch        *b_nTPCMCHits;   //!
   TBranch        *b_nDigi;   //!
   TBranch        *b_nBars;   //!
   TBranch        *b_nCompleteBars;   //!
   TBranch        *b_nBarEnds;   //!
   TBranch        *b_bars_id;   //!
   TBranch        *b_bars_ntrks;   //!
   TBranch        *b_bars_flag;   //!
   TBranch        *b_bars_edep;   //!
   TBranch        *b_bars_path;   //!
   TBranch        *b_bars_z;   //!
   TBranch        *b_bars_t;   //!
   TBranch        *b_bars_phi;   //!
   TBranch        *b_bars_atop;   //!
   TBranch        *b_bars_abot;   //!
   TBranch        *b_bars_ttop;   //!
   TBranch        *b_bars_tbot;   //!
   TBranch        *b_bars_ttop_twc;   //!
   TBranch        *b_bars_tbot_twc;   //!
   TBranch        *b_pairs_tof;   //!
   TBranch        *b_pairs_dphi;   //!
   TBranch        *b_pairs_dzeta;   //!
   TBranch        *b_pairs_dist;   //!
   TBranch        *b_TOF_MIN;   //!
   TBranch        *b_TOF_MAX;   //!
   TBranch        *b_TOF_MEAN;   //!
   TBranch        *b_TOF_STD;   //!
   TBranch        *b_DPHI_MIN;   //!
   TBranch        *b_DPHI_MAX;   //!
   TBranch        *b_DPHI_MEAN;   //!
   TBranch        *b_DPHI_STD;   //!
   TBranch        *b_DZETA_MIN;   //!
   TBranch        *b_DZETA_MAX;   //!
   TBranch        *b_DZETA_MEAN;   //!
   TBranch        *b_DZETA_STD;   //!
   TBranch        *b_DIST_MIN;   //!
   TBranch        *b_DIST_MAX;   //!
   TBranch        *b_DIST_MEAN;   //!
   TBranch        *b_DIST_STD;   //!

   simple_BarsFlagMC(TString, TTree *tree=0);
   virtual ~simple_BarsFlagMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef simple_BarsFlagMC_cxx
simple_BarsFlagMC::simple_BarsFlagMC(TString filename, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename);
      }
      f->GetObject("tDataBV",tree);

   }
   Init(tree);
}

simple_BarsFlagMC::~simple_BarsFlagMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t simple_BarsFlagMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t simple_BarsFlagMC::LoadTree(Long64_t entry)
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

void simple_BarsFlagMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   bars_id = 0;
   bars_ntrks = 0;
   bars_flag = 0;
   bars_edep = 0;
   bars_path = 0;
   bars_z = 0;
   bars_t = 0;
   bars_phi = 0;
   bars_atop = 0;
   bars_abot = 0;
   bars_ttop = 0;
   bars_tbot = 0;
   bars_ttop_twc = 0;
   bars_tbot_twc = 0;
   pairs_tof = 0;
   pairs_dphi = 0;
   pairs_dzeta = 0;
   pairs_dist = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("xo", &xo, &b_xo);
   fChain->SetBranchAddress("yo", &yo, &b_yo);
   fChain->SetBranchAddress("zo", &zo, &b_zo);
   fChain->SetBranchAddress("pxo", &pxo, &b_pxo);
   fChain->SetBranchAddress("pyo", &pyo, &b_pyo);
   fChain->SetBranchAddress("pzo", &pzo, &b_pzo);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("t_event", &t_event, &b_t_event);
   fChain->SetBranchAddress("pbar", &pbar, &b_pbar);
   fChain->SetBranchAddress("mc", &mc, &b_mc);
   fChain->SetBranchAddress("nTPCHits", &nTPCHits, &b_nTPCHits);
   fChain->SetBranchAddress("nTPCMCHits", &nTPCMCHits, &b_nTPCMCHits);
   fChain->SetBranchAddress("nDigi", &nDigi, &b_nDigi);
   fChain->SetBranchAddress("nBars", &nBars, &b_nBars);
   fChain->SetBranchAddress("nCompleteBars", &nCompleteBars, &b_nCompleteBars);
   fChain->SetBranchAddress("nBarEnds", &nBarEnds, &b_nBarEnds);
   fChain->SetBranchAddress("bars_id", &bars_id, &b_bars_id);
   fChain->SetBranchAddress("bars_ntrks", &bars_ntrks, &b_bars_ntrks);
   fChain->SetBranchAddress("bars_flag", &bars_flag, &b_bars_flag);
   fChain->SetBranchAddress("bars_edep", &bars_edep, &b_bars_edep);
   fChain->SetBranchAddress("bars_path", &bars_path, &b_bars_path);
   fChain->SetBranchAddress("bars_z", &bars_z, &b_bars_z);
   fChain->SetBranchAddress("bars_t", &bars_t, &b_bars_t);
   fChain->SetBranchAddress("bars_phi", &bars_phi, &b_bars_phi);
   fChain->SetBranchAddress("bars_atop", &bars_atop, &b_bars_atop);
   fChain->SetBranchAddress("bars_abot", &bars_abot, &b_bars_abot);
   fChain->SetBranchAddress("bars_ttop", &bars_ttop, &b_bars_ttop);
   fChain->SetBranchAddress("bars_tbot", &bars_tbot, &b_bars_tbot);
   fChain->SetBranchAddress("bars_ttop_twc", &bars_ttop_twc, &b_bars_ttop_twc);
   fChain->SetBranchAddress("bars_tbot_twc", &bars_tbot_twc, &b_bars_tbot_twc);
   fChain->SetBranchAddress("pairs_tof", &pairs_tof, &b_pairs_tof);
   fChain->SetBranchAddress("pairs_dphi", &pairs_dphi, &b_pairs_dphi);
   fChain->SetBranchAddress("pairs_dzeta", &pairs_dzeta, &b_pairs_dzeta);
   fChain->SetBranchAddress("pairs_dist", &pairs_dist, &b_pairs_dist);
   fChain->SetBranchAddress("tof_min", &tof_min, &b_TOF_MIN);
   fChain->SetBranchAddress("tof_max", &tof_max, &b_TOF_MAX);
   fChain->SetBranchAddress("tof_mean", &tof_mean, &b_TOF_MEAN);
   fChain->SetBranchAddress("tof_std", &tof_std, &b_TOF_STD);
   fChain->SetBranchAddress("dphi_min", &dphi_min, &b_DPHI_MIN);
   fChain->SetBranchAddress("dphi_max", &dphi_max, &b_DPHI_MAX);
   fChain->SetBranchAddress("dphi_mean", &dphi_mean, &b_DPHI_MEAN);
   fChain->SetBranchAddress("dphi_std", &dphi_std, &b_DPHI_STD);
   fChain->SetBranchAddress("dzeta_min", &dzeta_min, &b_DZETA_MIN);
   fChain->SetBranchAddress("dzeta_max", &dzeta_max, &b_DZETA_MAX);
   fChain->SetBranchAddress("dzeta_mean", &dzeta_mean, &b_DZETA_MEAN);
   fChain->SetBranchAddress("dzeta_std", &dzeta_std, &b_DZETA_STD);
   fChain->SetBranchAddress("dist_min", &dist_min, &b_DIST_MIN);
   fChain->SetBranchAddress("dist_max", &dist_max, &b_DIST_MAX);
   fChain->SetBranchAddress("dist_mean", &dist_mean, &b_DIST_MEAN);
   fChain->SetBranchAddress("dist_std", &dist_std, &b_DIST_STD);
   Notify();
}

Bool_t simple_BarsFlagMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void simple_BarsFlagMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t simple_BarsFlagMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef simple_BarsFlagMC_cxx
