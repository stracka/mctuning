#include "simple.C"
#include "TLegend.h"

#define NBARS 64

void DrawDistrs(int  filenumber = 0)
{
  string filename;

  if(filenumber == 0) filename = "DATA_BV_TREE_cosmics_run_output07947_old_thrs800.root";
  if(filenumber == 1) filename = "DATA_BV_TREE_cosmics_run_output07947_thrs800.root";
  if(filenumber == 2) filename = "DATA_BV_TREE_cosmics_run_output08768_thrs800.root";
  if(filenumber == 3) filename = "DATA_BV_TREE_cosmics_run_output08926_thrs800.root";
  if(filenumber == 4) filename = "DATA_BV_TREE_cosmics_run_output08926_std_sw_thrs800.root";
  if(filenumber == 5) filename = "DATA_BV_TREE_cosmics_run_output08926_low_sw_thrs300.root";
  if(filenumber == 6) filename = "DATA_BV_TREE_cosmics_run_output08925_low_sw_thrs300.root";
  if(filenumber == 7) filename = "DATA_BV_TREE_cosmics_runs_8925_8926_low_sw_thrs300.root";
  if(filenumber == 8) filename = "DATA_BV_TREE_cosmics_run_output08925_sw_thrs500.root";
  if(filenumber == 9) filename = "DATA_BV_TREE_cosmics_run_output08926_sw_thrs500.root";
  if(filenumber == 10) filename = "DATA_BV_TREE_cosmics_runs_8925_8926_sw_thrs500.root";

  cout<<endl;
  cout<<"Filename: "<<filename<<endl;
  cout<<endl;

  TH1F *hatop[NBARS];
  TH1F *habot[NBARS];
  TH1F *hlog[NBARS];
  TH1F *httop_tbot[NBARS];
  TH2F *hlog_vs_ttop_tbot[NBARS];
  TH2F *hratio_vs_ttop_tbot[NBARS];
  
  for (int i = 0; i < NBARS; i++)
  {
    hatop[i] = new TH1F(Form("hatop%d", i), Form("hatop%d", i), 500, 0, 1);
    hatop[i]->SetTitle("Amplitude top ADCs");
    hatop[i]->SetTitleSize(0.07262164);
    hatop[i]->GetXaxis()->SetTitle("A_{top} (V)");
    hatop[i]->GetXaxis()->SetTitleSize(0.045);
    hatop[i]->SetLineColor(kRed);
    hatop[i]->SetLineWidth(2);

    habot[i] = new TH1F(Form("habot%d", i), Form("habot%d", i), 500, 0, 1);
    habot[i]->SetTitle("Amplitude bot ADCs");
    habot[i]->SetTitleSize(0.07262164);
    habot[i]->GetXaxis()->SetTitle("A_{bot} (V)");
    habot[i]->GetXaxis()->SetTitleSize(0.045);
    habot[i]->SetLineColor(kRed);
    habot[i]->SetLineWidth(2);

    hlog[i] = new TH1F(Form("hlog%d", i), Form("hlog%d", i), 400, -4, 4);
    hlog[i]->SetTitle("log of bot/top Amplitudes");
    hlog[i]->SetTitleSize(0.07262164);
    hlog[i]->GetXaxis()->SetTitle("log(A_{bot}/A_{top})");
    hlog[i]->GetXaxis()->SetTitleSize(0.045);
    hlog[i]->SetLineColor(kRed);
    hlog[i]->SetLineWidth(2);

    httop_tbot[i] = new TH1F(Form("httop_tbot%d", i), Form("httop_tbot%d", i), 500, -3e-8, 3e-8);
    httop_tbot[i]->SetTitle("time TDCs top - TDCs bot");
    httop_tbot[i]->SetTitleSize(0.07262164);
    httop_tbot[i]->GetXaxis()->SetTitle("t_{top}-t_{bot} (s)");
    httop_tbot[i]->GetXaxis()->SetTitleSize(0.045);
    httop_tbot[i]->SetLineColor(kRed);
    httop_tbot[i]->SetLineWidth(2);
  
    hlog_vs_ttop_tbot[i] = new TH2F(Form("hlog_vs_ttop_tbot%d", i), Form("hlog_vs_ttop_tbot%d", i), 100, -3e-8, 3e-8, 100, -4, 4);
    hlog_vs_ttop_tbot[i]->SetTitle("log ratio vs time");
    hlog_vs_ttop_tbot[i]->SetTitleSize(0.07262164);
    hlog_vs_ttop_tbot[i]->GetXaxis()->SetTitle("t_{top}-t_{bot} (s)");
    hlog_vs_ttop_tbot[i]->GetXaxis()->SetTitleSize(0.045);
    hlog_vs_ttop_tbot[i]->GetYaxis()->SetTitle("log(A_{bot}/A_{top})");
    hlog_vs_ttop_tbot[i]->GetYaxis()->SetTitleSize(0.045);

    hratio_vs_ttop_tbot[i] = new TH2F(Form("hratio_vs_ttop_tbot%d", i), Form("hratio_vs_ttop_tbot%d", i), 100, -3e-8, 3e-8, 100, 0, 25);
    hratio_vs_ttop_tbot[i]->SetTitle("ratio vs time");
    hratio_vs_ttop_tbot[i]->SetTitleSize(0.07262164);
    hratio_vs_ttop_tbot[i]->GetXaxis()->SetTitle("t_{top}-t_{bot} (s)");
    hratio_vs_ttop_tbot[i]->GetXaxis()->SetTitleSize(0.045);
    hratio_vs_ttop_tbot[i]->GetYaxis()->SetTitle("A_{bot}/A_{top}");
    hratio_vs_ttop_tbot[i]->GetYaxis()->SetTitleSize(0.045);
    
    
  }

  TH1F *hn;
  hn = new TH1F("hn", "hn", 12, 1, 13);
  hn->SetTitle("Event multiplicity");
  hn->SetTitleSize(0.07262164);
  hn->GetXaxis()->SetTitle("#Bars hit");
  hn->GetXaxis()->SetTitleSize(0.045);
  hn->SetLineColor(kRed);
  hn->SetLineWidth(2);

  simple test(filename.c_str(), 0);

  for (int i = 0; i < test.treeMCBV->GetEntries(); i++)
  {
    test.GetEntry(i);

    std::vector<double> cal_atop;
    std::vector<double> cal_abot;
    std::vector<double> cal_ttop_tbot;

    std::vector<int> cal_id;

    for (int j = 0; j < test.bars_id->size(); j++)
    {
      double atop = test.bars_atop->at(j);
      double abot = test.bars_abot->at(j);

      double ttop = test.bars_ttop->at(j);
      double tbot = test.bars_tbot->at(j);
      double ttop_tbot = ttop - tbot;

      cal_atop.push_back(atop);
      cal_abot.push_back(abot);
      cal_ttop_tbot.push_back(ttop_tbot);

      cal_id.push_back(test.bars_id->at(j));

      hatop[test.bars_id->at(j)]->Fill(cal_atop.at(j));
      habot[test.bars_id->at(j)]->Fill(cal_abot.at(j));
      hlog[test.bars_id->at(j)]->Fill(TMath::Log(cal_abot.at(j) / cal_atop.at(j)));

      httop_tbot[test.bars_id->at(j)]->Fill(cal_ttop_tbot.at(j));

      hlog_vs_ttop_tbot[test.bars_id->at(j)]->Fill(cal_ttop_tbot.at(j),TMath::Log(cal_abot.at(j) / cal_atop.at(j)));
      hratio_vs_ttop_tbot[test.bars_id->at(j)]->Fill(cal_ttop_tbot.at(j), cal_abot.at(j) / cal_atop.at(j));
    }

    hn->Fill(cal_id.size());
  }

  // gStyle->SetOptStat(0);

  TCanvas *c[NBARS];

  for (int i = 0; i < NBARS; i++)
  {
    c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 1800, 1200);
    c[i]->Divide(3, 2);

    c[i]->cd(1);
    hatop[i]->DrawNormalized();

    c[i]->cd(2);
    habot[i]->DrawNormalized();
    
    c[i]->cd(3);
    hn->DrawNormalized();
    
    c[i]->cd(4);
    hlog[i]->DrawNormalized();
    
    c[i]->cd(5);
    httop_tbot[i]->DrawNormalized();

    c[i]->cd(6);
    hlog_vs_ttop_tbot[i]->Draw("colz");    
    //c[i].SaveAs(Form("status_canvas%d.png",i));
  }

  return;

}
