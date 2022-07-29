#ifndef simple_h
#define simple_h

#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <bitset>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// // Header file for the classes stored in the TTree if any.
// #include "TClonesArray.h"
// #include "Riostream.h"

class simple {
public :
   TTree           *treeMCBV;
   std::string     filename_core;
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t           pbar_file_flag;
   Int_t           mc_file_flag;
   Int_t           event_selection;
   Long64_t        n_tot_events, n_selected;
   ///< parameters for histos
   Double_t        t_event_max;
   Double_t        A_max;
   // Declaration of leaf types
   Int_t           event;
   Double_t        t_event;
   Bool_t          pbar; ///< 0 = cosmic, 1 = pbar annihilation
   Bool_t          mc;
   ///< MC generation info (primary particle information: pos, mom)
   Double_t xo, yo, zo, pxo, pyo, pzo;
   ///< TPC overall variables
   Int_t nTPCHits; ///< Number of TPC hits
   Int_t nTPCMCHits; ///< Number of TPC MC hits
   ///< MLU Trigger logic (4 bars OR)
   std::vector<uint16_t> mlu_table; // veto MLU table
   Bool_t BVMLU[16];
   ///< Barrel Veto (BV) variables
   Int_t           nDigi;
   Int_t           nBars;
   Int_t           nBarEnds;
   vector<Int_t>     *bars_id;
   vector<Int_t>     *bars_ntrks;
   vector<Double_t>   *bars_edep;
   vector<Double_t>   *bars_path;
   vector<Double_t>   *bars_z;
   vector<Double_t>   *bars_t;
   vector<Double_t>   *bars_phi;
   vector<Double_t>   *bars_atop;
   vector<Double_t>   *bars_abot;
   vector<Double_t>   *bars_ttop;
   vector<Double_t>   *bars_tbot;
   vector<Double_t>   *pairs_tof;
   vector<Double_t>   *pairs_dphi;
   vector<Double_t>   *pairs_dzeta;
   vector<Double_t>   *pairs_dist;
   Double_t         tof_min;
   Double_t         tof_max;
   Double_t         tof_mean;
   Double_t         tof_std;
   Double_t         dphi_min;
   Double_t         dphi_max;
   Double_t         dphi_mean;
   Double_t         dphi_std;
   Double_t         dzeta_min;
   Double_t         dzeta_max;
   Double_t         dzeta_mean;
   Double_t         dzeta_std;
   Double_t         dist_min;
   Double_t         dist_max;
   Double_t         dist_mean;
   Double_t         dist_std;
   // List of branches
   TBranch        *b_xo;   //!
   TBranch        *b_yo;   //!
   TBranch        *b_zo;   //!
   TBranch        *b_pxo;  //!
   TBranch        *b_pyo;  //!
   TBranch        *b_pzo;  //!
   TBranch        *b_event;   //!
   TBranch        *b_t_event;   //!
   TBranch        *b_pbar;   //!
   TBranch        *b_mc;   //!
   TBranch        *b_nTPCHits; //!
   TBranch        *b_nTPCMCHits; //!
   TBranch        *b_nDigi;   //!
   TBranch        *b_nBars;   //!
   TBranch        *b_nBarEnds;   //!
   TBranch        *b_bars_id;   //!
   TBranch        *b_bars_ntrks;   //!
   TBranch        *b_bars_edep;   //!
   TBranch        *b_bars_path;   //!
   TBranch        *b_bars_z;   //!
   TBranch        *b_bars_t;   //!
   TBranch        *b_bars_phi;   //!
   TBranch        *b_bars_atop;   //!
   TBranch        *b_bars_abot;   //!
   TBranch        *b_bars_ttop;   //!
   TBranch        *b_bars_tbot;   //!
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


   ///< single bar histos
  TH1F *hnDigi[2], *hnBars[2], *hEdep[2], *hEvTotEne[2], *hBarID[2];
  TH1F *ht[2], *hZ[2];
  TH1F *hATop[2], *hABot[2], *hlnA[2], *hsqrtA[2], *hdt[2]; ///< histos for Atop, Abot, tTop, tBot;
  ///< pair of bar histos
  TH2F *hTOFvsDIST[2];
  TH1F *hTOF[2], *hDIST[2], *hDz[2], *hDPhi[2];
  TH1F *hTOFMIN[2], *hTOFMAX[2], *hTOFMEAN[2], *hTOFSTD[2];
  TH1F *hDPHIMIN[2], *hDPHIMAX[2], *hDPHIMEAN[2], *hDPHISTD[2];
  TH1F *hDZETAMIN[2], *hDZETAMAX[2], *hDZETAMEAN[2], *hDZETASTD[2];
  TH1F *hDISTMIN[2], *hDISTMAX[2], *hDISTMEAN[2], *hDISTSTD[2];
  
   simple(std::string, Int_t);
   virtual ~simple();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadMCBVTREE(Long64_t entry);
   virtual void     InitMCBVTREE(TTree *tree);
   virtual void     ShowBV();

   virtual void   SetMLUTable(string);
   virtual void   CreateHistos();
   virtual void   FillHistos();
   virtual Bool_t EventSelection();
   virtual void   ShowHistos();
   uint16_t       bool_array_to_uint16 (bool ar[], size_t ar_size);
   bool           mlu(uint16_t);
   bool           findInVector(const std::vector<uint16_t>  &, const uint16_t  &);

};

#endif

#ifdef simple_cxx

simple::simple(std::string file_name, Int_t selection) : treeMCBV(0) 
{


   event_selection = selection;

   pbar_file_flag = -1;
   if(file_name.find("pbar") != std::string::npos) pbar_file_flag = 1;
   if(file_name.find("cosmic") != std::string::npos) pbar_file_flag = 0;
   mc_file_flag = -1;
   if(file_name.find("MC") != std::string::npos) mc_file_flag = 1;
   if(file_name.find("DATA") != std::string::npos) mc_file_flag = 0;
   // used to generate this class and read the Tree.

   TTree *tree = 0;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file_name.c_str());
   if (!f || !f->IsOpen()) {
      f = new TFile(file_name.c_str());
   }

   f->GetObject("tDataBV",tree);
   InitMCBVTREE(tree);

   ///< Set the filename_core
   filename_core=file_name;
   Int_t dot, slash, len;
   slash = filename_core.find_first_of('/');
   len = filename_core.size();
   filename_core = filename_core.substr(slash+1,len);
   dot = filename_core.find_last_of('.');
   filename_core = filename_core.substr(0,dot);

}

simple::~simple()
{
   if (!treeMCBV) return;
   delete treeMCBV->GetCurrentFile();
}

Int_t simple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!treeMCBV) return 0;
   return treeMCBV->GetEntry(entry);
}

Long64_t simple::LoadMCBVTREE(Long64_t entry)
{
// Set the environment to read one entry
   if (!treeMCBV) return -5;
   Long64_t centry = treeMCBV->LoadTree(entry);
   if (centry < 0) return centry;
   if (treeMCBV->GetTreeNumber() != fCurrent) {
      fCurrent = treeMCBV->GetTreeNumber();
   }
   return centry;
}

void simple::InitMCBVTREE(TTree *tree)
{
   // Set object pointer
   bars_id = 0;
   bars_ntrks = 0;
   bars_edep = 0;
   bars_path = 0;
   bars_z = 0;
   bars_t = 0;
   bars_phi = 0;
   bars_atop = 0;
   bars_abot = 0;
   bars_ttop = 0;
   bars_tbot = 0;
   pairs_tof = 0;
   pairs_dphi = 0;
   pairs_dzeta = 0;
   pairs_dist = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   treeMCBV = tree;
   fCurrent = -1;
   treeMCBV->SetMakeClass(1);
   if(mc_file_flag) {
   ///< Generation information on primary particle (pos, mom)
      treeMCBV->SetBranchAddress("xo",    &xo,   &b_xo);
      treeMCBV->SetBranchAddress("yo",    &yo,   &b_yo);
      treeMCBV->SetBranchAddress("zo",    &zo,   &b_zo);
      treeMCBV->SetBranchAddress("pxo",   &pxo,  &b_pxo);
      treeMCBV->SetBranchAddress("pyo",   &pyo,  &b_pyo);
      treeMCBV->SetBranchAddress("pzo",   &pzo,  &b_pzo);
   }
   treeMCBV->SetBranchAddress("event", &event, &b_event);
   treeMCBV->SetBranchAddress("t_event", &t_event, &b_t_event);
   treeMCBV->SetBranchAddress("pbar", &pbar, &b_pbar);
   treeMCBV->SetBranchAddress("mc", &mc, &b_mc);
   if(mc_file_flag) {
      ///< TPC info
      treeMCBV->SetBranchAddress("nTPCHits",  &nTPCHits,   &b_nTPCHits);
      treeMCBV->SetBranchAddress("nTPCMCHits",&nTPCMCHits, &b_nTPCMCHits);
   }
   treeMCBV->SetBranchAddress("nDigi", &nDigi, &b_nDigi);
   treeMCBV->SetBranchAddress("nBars", &nBars, &b_nBars);
   treeMCBV->SetBranchAddress("nBarEnds", &nBarEnds, &b_nBarEnds);
   treeMCBV->SetBranchAddress("bars_id", &bars_id, &b_bars_id);
   treeMCBV->SetBranchAddress("bars_ntrks", &bars_ntrks, &b_bars_ntrks);
   treeMCBV->SetBranchAddress("bars_edep", &bars_edep, &b_bars_edep);
   treeMCBV->SetBranchAddress("bars_path", &bars_path, &b_bars_path);
   treeMCBV->SetBranchAddress("bars_z", &bars_z, &b_bars_z);
   treeMCBV->SetBranchAddress("bars_t", &bars_t, &b_bars_t);
   treeMCBV->SetBranchAddress("bars_phi", &bars_phi, &b_bars_phi);
   treeMCBV->SetBranchAddress("bars_atop", &bars_atop, &b_bars_atop);
   treeMCBV->SetBranchAddress("bars_abot", &bars_abot, &b_bars_abot);
   treeMCBV->SetBranchAddress("bars_ttop", &bars_ttop, &b_bars_ttop);
   treeMCBV->SetBranchAddress("bars_tbot", &bars_tbot, &b_bars_tbot);
   treeMCBV->SetBranchAddress("pairs_tof", &pairs_tof, &b_pairs_tof);
   treeMCBV->SetBranchAddress("pairs_dphi", &pairs_dphi, &b_pairs_dphi);
   treeMCBV->SetBranchAddress("pairs_dzeta", &pairs_dzeta, &b_pairs_dzeta);
   treeMCBV->SetBranchAddress("pairs_dist", &pairs_dist, &b_pairs_dist);
   treeMCBV->SetBranchAddress("tof_min", &tof_min, &b_TOF_MIN);
   treeMCBV->SetBranchAddress("tof_max", &tof_max, &b_TOF_MAX);
   treeMCBV->SetBranchAddress("tof_mean", &tof_mean, &b_TOF_MEAN);
   treeMCBV->SetBranchAddress("tof_std", &tof_std, &b_TOF_STD);
   treeMCBV->SetBranchAddress("dphi_min", &dphi_min, &b_DPHI_MIN);
   treeMCBV->SetBranchAddress("dphi_max", &dphi_max, &b_DPHI_MAX);
   treeMCBV->SetBranchAddress("dphi_mean", &dphi_mean, &b_DPHI_MEAN);
   treeMCBV->SetBranchAddress("dphi_std", &dphi_std, &b_DPHI_STD);
   treeMCBV->SetBranchAddress("dzeta_min", &dzeta_min, &b_DZETA_MIN);
   treeMCBV->SetBranchAddress("dzeta_max", &dzeta_max, &b_DZETA_MAX);
   treeMCBV->SetBranchAddress("dzeta_mean", &dzeta_mean, &b_DZETA_MEAN);
   treeMCBV->SetBranchAddress("dzeta_std", &dzeta_std, &b_DZETA_STD);
   treeMCBV->SetBranchAddress("dist_min", &dist_min, &b_DIST_MIN);
   treeMCBV->SetBranchAddress("dist_max", &dist_max, &b_DIST_MAX);
   treeMCBV->SetBranchAddress("dist_mean", &dist_mean, &b_DIST_MEAN);
   treeMCBV->SetBranchAddress("dist_std", &dist_std, &b_DIST_STD);

}

uint16_t simple::bool_array_to_uint16 (bool ar[], size_t ar_size) {
    uint16_t ret {};

    for (size_t i = 0; i < ar_size; ++i) {
        uint16_t s {ar[i]};
        s <<= i;
        ret |= s;
    }
    return ret;
}

bool simple::mlu(uint16_t data) {
   //    return findInVector<uint16_t>(comptable, data);
   return !findInVector(mlu_table, data);   
}

bool simple::findInVector(const std::vector<uint16_t>  & vecOfElements, const uint16_t  & element)
{
  bool result;
  // Find given element in vector
  auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
  if (it != vecOfElements.end())
    result = true;
  else
    result = false;
  return result;
}

void simple::SetMLUTable(string mlufile) {
    std::vector<uint16_t> comptable;   // trigger table 
    // fill trigger table (complementary of veto table)
    std::ifstream rfile(mlufile.c_str());
    string pattern; 
    while (rfile >> pattern){
      comptable.push_back((uint16_t)(strtol(pattern.c_str(),NULL,16)));
    }
    // fill veto table
    for (int i = 0; i<65536; i++){
      uint16_t pat = (uint16_t)i; 
      if (!findInVector(comptable, pat))
	    mlu_table.push_back(pat);
    }
    std::cout << mlu_table.size() << " " << comptable.size() << std::endl;
}

#endif // #ifdef simple_cxx
