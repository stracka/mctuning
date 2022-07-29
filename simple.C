#define simple_cxx
#include "simple.h"
#include <TH2.h>
#include <TRandom.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>

///< CONSTANTS
static const Int_t   NBARS=64;
static const Double_t C        = 299792458*1.e-9; //< speed of light (m/ns)
static const Double_t V_EJ200    = C/1.93; //< Accordint to Gareth (private communication)
static const Double_t s_to_ns = 1.e9;

///< FLAGS for CANVAS to be shown
static const bool bSummary       = false;
static const bool bTopBot        = true;
static const bool bPairs         = false;
static const bool bMinMaxMeanSTD = false;
static const bool bTOFDIST       = false;
///< FLAG to write the canvas output to a file
static const bool c_save         = false;
///< ############################################################
///< HISTOS AND NUMBERS FOR THE BV BARS
///< ############################################################
void simple::ShowBV() {

   if (treeMCBV == 0) return;
   n_tot_events = treeMCBV->GetEntriesFast();
   n_selected = 0;

   if(event_selection==3) SetMLUTable("./mlu_pbars_bottom_trap.dat");
   CreateHistos();

   for (Long64_t ievt=0; ievt<n_tot_events;ievt++) { //loop on Events
      LoadMCBVTREE(ievt);
      treeMCBV->GetEntry(ievt);
      FillHistos();
   }//end loop on events
   
   if(n_tot_events>0) {
     // PrepareHistos();
     // ShowHistos(); 
   }
   
//________________________________________________________________________

}

void simple::ShowHistos() {

   ///< Defining canvas
   TCanvas *cTOF, *cDPHI, *cDZETA, *cDIST;
   TCanvas *cSummary, *cTOFDIST, *cPairs, *cTopBot;

   if(bTopBot) {cTopBot = new TCanvas("cTopBot","cTopBot", 1600, 1200); cTopBot->Divide(3,2);}

   std::string  stat_opt = "imr";     // i = integral, m = mean, r = rms, etc. 
   gStyle->SetOptStat(stat_opt.c_str());
   gStyle->SetOptFit(111111111);
   std::string draw_opt = "H,SAME";
   std::ostringstream s_c_out;
   for(unsigned int i=0; i<2; i++) {
      if(hnBars[i]->GetEntries()<=0) continue;


      if(bTopBot) {
         cTopBot->cd(1);      hATop[i]->Draw(draw_opt.c_str());
         cTopBot->cd(2);      hABot[i]->Draw(draw_opt.c_str());
         cTopBot->cd(3);      hsqrtA[i]->Draw(draw_opt.c_str());
         cTopBot->cd(4);      hlnA[i]->Draw(draw_opt.c_str());
         cTopBot->cd(6);      hnBars[i]->Draw(draw_opt.c_str());
      }

   }
   if(c_save) {
     s_c_out.clear(); s_c_out.str("");
     s_c_out << "root_output_files/" << filename_core << "_cTopBot.root";
     TFile* fout = new TFile(s_c_out.str().c_str(),"RECREATE");
     fout->cd();
     if(bTopBot) {
       for(unsigned int i=0; i<2; i++) {
	 hATop[i]->Write(); 
	 hABot[i]->Write(); 
	 hsqrtA[i]->Write(); 
	 hlnA[i]->Write(); 
	 hnBars[i]->Write(); 		 
       }      
     }
     fout->Close(); 

   }

   std::cout << "Total number of events:  " << n_tot_events << std::endl;
   if(n_tot_events!=0) {
      std::cout << "Total number of selected events: " << n_selected;
      std::cout << "(" << (int)(n_selected*100.)/n_tot_events << "\%)" << std::endl;
   }
}

void simple::FillHistos() {
   Bool_t fill = EventSelection();
   if(!fill) return;
   n_selected++;

   hnBars[pbar]->Fill(nBars);

   for(unsigned int i = 0; i<bars_id->size(); i++) {
      hATop[pbar]->Fill(bars_atop->at(i));
      hABot[pbar]->Fill(bars_abot->at(i));
      hlnA[pbar]->Fill(0.5*TMath::Log(bars_atop->at(i)/bars_abot->at(i)));
      hsqrtA[pbar]->Fill(TMath::Sqrt(bars_atop->at(i)*bars_abot->at(i)));
   }

}

Bool_t simple::EventSelection() {
   ///< To simulate the trigger [they work only on MC]
   ///< 0 -> no event selection
   ///< 1 => BV grand OR (any two ends) [nBarEnds>=2] [is it working?]
   ///< 2 => TPC grand OR [is it working?]
   ///< 3 => MLU [is it working?]
   ///< For other event selection 
   ///< 11 => nBars >=1

   ///< No selection
   if(event_selection == 0) return true;

   ///< event_selection < 10 => TRIGGER EMULATION => ONLY FOR MC
   if(event_selection<10) {
      if(!mc_file_flag) return true; 
      ///< BV grand OR
      if(event_selection==1) { ///< BV GrandOR?
         if(nBarEnds>=2) {return true;} else {return false;}
      }
      if(event_selection==2) { ///< TPC GrandOR?
         if(nTPCHits>200) {return true;} else {return false;}
      }
      if(event_selection==3) { ///< MLU2+
         for(unsigned int i=0; i<16; i++) {BVMLU[i] = false;}
         for(unsigned int i = 0; i<bars_id->size(); i++) {
            unsigned int pad = (int)(bars_id->at(i)/4);
            BVMLU[pad] = 1;
         }
         uint16_t necklace = bool_array_to_uint16(BVMLU,16);
         return mlu(necklace);
      }
   }
   /*
   ///< From now on => event_selection > 10
   if(event_selection==11) { ///< 11
      if(bars_id->size()>=1) {return true;} else {return false;}
   }
   */
   std::cout << "WARNING: No recognized event_selection flag (" << event_selection << ") => REJECTING THE EVENT" << std::endl;
   return false;

}


void simple::CreateHistos() {
   ostringstream name, title;
   Int_t colors[2]; colors[0] = 8; colors[1] = 38; ///< Color code for cosmics and pbars [38/46]
   Int_t stile[2]; stile[0] = 3001; stile[1] = 3008; ///< Fill style code for cosmics and pbars

   for(unsigned int i=0; i<2; i++) {

      name.clear(); name.str(""); title.clear(); title.str("");
      name << "hnBars" << i; title << "Number of bars per event "; if(i==0) {title << "(cosmic)";} else {title << "(pbar)";}
      hnBars[i]      = new TH1F(name.str().c_str(),title.str().c_str(), 15, -0.5, 14.5);
      hnBars[i]->SetFillColor(colors[i]); hnBars[i]->SetFillStyle(stile[i]); 
      hnBars[i]->Sumw2(); 
      
      name.clear(); name.str(""); title.clear(); title.str("");
      name << "hATop" << i; title << "A_{top} "; if(i==0) {title << "(cosmic)";} else {title << "(pbar)";}
      hATop[i]          = new TH1F(name.str().c_str(),title.str().c_str(), 100, 0., 2.);
      hATop[i]->SetXTitle("a.u.");
      hATop[i]->SetFillColor(colors[i]); hATop[i]->SetFillStyle(stile[i]);
      hATop[i]->Sumw2(); 
      
      name.clear(); name.str(""); title.clear(); title.str("");
      name << "hABot" << i; title << "A_{bot} "; if(i==0) {title << "(cosmic)";} else {title << "(pbar)";}
      hABot[i]          = new TH1F(name.str().c_str(),title.str().c_str(), 100, 0., 2.);
      hABot[i]->SetXTitle("a.u.");
      hABot[i]->SetFillColor(colors[i]); hABot[i]->SetFillStyle(stile[i]);
      hABot[i]->Sumw2(); 
      
      
      name.clear(); name.str(""); title.clear(); title.str("");
      name << "hsqrtA" << i; title << "#sqrt{A_{top} #times A_{bot}} #propto #sqrt{e^{-L/#lambda}} "; if(i==0) {title << "(cosmic)";} else {title << "(pbar)";}
      hsqrtA[i]          = new TH1F(name.str().c_str(),title.str().c_str(), 100, 0., 2.);
      hsqrtA[i]->SetXTitle("a.u.");
      hsqrtA[i]->SetFillColor(colors[i]); hsqrtA[i]->SetFillStyle(stile[i]);
      hsqrtA[i]->Sumw2(); 

      name.clear(); name.str(""); title.clear(); title.str("");
      name << "hlnA" << i; title << "#frac{1}{2} log(A_{top}/A_{bot}) = #frac{z_{hit}}{#lambda} "; if(i==0) {title << "(cosmic)";} else {title << "(pbar)";}
      hlnA[i]          = new TH1F(name.str().c_str(),title.str().c_str(), 100, -2., +2.);
      hlnA[i]->SetXTitle("a.u.");
      hlnA[i]->SetFillColor(colors[i]); hlnA[i]->SetFillStyle(stile[i]);

   }

}

