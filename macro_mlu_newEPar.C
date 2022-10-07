#include "simple.C"
#include "TLegend.h"

double macro_mlu_newEPar(double mlu_scale = 0.5,
                         double top_gain = 0.1944,
                         double qE_top = 0.2, 
                         double kE_top = 0.07, 
                         double bot_gain = 0.1944,
                         double qE_bot = 0.2,
                         double kE_bot = 0.07, 
                         double gain_c = 0.1,
                         double lambda = 1.3,
                         double c_Ej = 299792458 / 1.93,
                         double t_smearing = 500e-12,
                         int filenumber = 0) // introduci anche q nella risoluzione in energia
{
  TString filename;

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

  int pointIndex = filename.Index(".");
  int thrsindex = filename.Index("thrs");
  TString Sthreshold = (TString)(filename(thrsindex+4,pointIndex-(thrsindex+4)));
  double threshold = Sthreshold.Atof();

  double thres = threshold/TMath::Power(2.0,14);
  double zmin = -1.3;
  double zmax = 1.3;

  // double gain_top[64];
  // double gain_bot[64];

  // for (int i = 0; i < 64; i++)
  // {
  //   gain_top[i] = top_gain + gRandom->Gaus(0., top_gain * top_spread);
  //   gain_bot[i] = bot_gain + gRandom->Gaus(0., bot_gain * bot_spread);
  // }

  // arXiv: 1905.06032 SiPM + Scint non linearity
  // https://ieeexplore.ieee.org/document/5874113 resolution
  // deltaE/E = a E^-b
  // <Q> = gain/c * ( 1 - exp(-cE))

  double err_scale = 0.2;

  string gname[2] = {"MC_BV_TREE_cosmics_ecomug_hsphere_1M_B1T_set999.root", (string)filename};
  bool isMC[2] = {true, false};

  TH1F *hatop[2];
  hatop[0] = new TH1F("hatop0", "hatop0", 50, 0, 1);
  hatop[1] = new TH1F("hatop1", "hatop1", 50, 0, 1);
  hatop[1]->SetTitle("Amplitude top ADCs");
  hatop[1]->SetTitleSize(0.07262164);
  hatop[1]->GetXaxis()->SetTitle("A_{top} (V)");
  hatop[1]->GetXaxis()->SetTitleSize(0.045);
  hatop[1]->SetLineColor(kRed);
  hatop[0]->SetLineWidth(2);
  hatop[1]->SetLineWidth(2);

  TH1F *habot[2];
  habot[0] = new TH1F("habot0", "habot0", 50, 0, 1);
  habot[1] = new TH1F("habot1", "habot1", 50, 0, 1);
  habot[1]->SetTitle("Amplitude bot ADCs");
  habot[1]->SetTitleSize(0.07262164);
  habot[1]->GetXaxis()->SetTitle("A_{bot} (V)");
  habot[1]->GetXaxis()->SetTitleSize(0.045);
  habot[1]->SetLineColor(kRed);
  habot[0]->SetLineWidth(2);
  habot[1]->SetLineWidth(2);

  TH1F *httop_tbot[2];
  httop_tbot[0] = new TH1F("httop_tbot0", "httop_tbot0", 50, -3e-8, 3e-8);
  httop_tbot[1] = new TH1F("httop_tbot1", "httop_tbot1", 50, -3e-8, 3e-8);
  httop_tbot[1]->SetTitle("time TDCs bot - TDCs top");
  httop_tbot[1]->SetTitleSize(0.07262164);
  httop_tbot[1]->GetXaxis()->SetTitle("t_{bot}-t_{top} (s)");
  httop_tbot[1]->GetXaxis()->SetTitleSize(0.045);
  httop_tbot[1]->SetLineColor(kRed);
  httop_tbot[0]->SetLineWidth(2);
  httop_tbot[1]->SetLineWidth(2);


  TH1F *hn[2];
  hn[0] = new TH1F("hn0", "hn0", 12, 1, 13);
  hn[1] = new TH1F("hn1", "hn1", 12, 1, 13);
  hn[1]->SetTitle("Event multiplicity");
  hn[1]->SetTitleSize(0.07262164);
  hn[1]->GetXaxis()->SetTitle("#Bars hit");
  hn[1]->GetXaxis()->SetTitleSize(0.045);
  hn[1]->SetLineColor(kRed);
  hn[0]->SetLineWidth(2);
  hn[1]->SetLineWidth(2);


  TH1F *hlog[2];
  hlog[0] = new TH1F("hlog0", "hlog0", 40, -4, 4);
  hlog[1] = new TH1F("hlog1", "hlog1", 40, -4, 4);
  hlog[1]->SetTitle("log of top/bot Amplitudes");
  hlog[1]->SetTitleSize(0.07262164);
  hlog[1]->GetXaxis()->SetTitle("log(A_{top}/A_{bot})");
  hlog[1]->GetXaxis()->SetTitleSize(0.045);
  hlog[1]->SetLineColor(kRed);
  hlog[0]->SetLineWidth(2);
  hlog[1]->SetLineWidth(2);

  TH2F *hlog_vs_ttop_tbot[2];
  hlog_vs_ttop_tbot[0] = new TH2F("h0log_vs_ttop_tbot", "h0log_vs_ttop_tbot", 100, -3e-8, 3e-8, 100, -4, 4);
  hlog_vs_ttop_tbot[1] = new TH2F("h1log_vs_ttop_tbot", "h1log_vs_ttop_tbot", 100, -3e-8, 3e-8, 100, -4, 4);
  hlog_vs_ttop_tbot[1]->SetLineColor(kRed);

  TH2F *hratio_vs_ttop_tbot[2];
  hratio_vs_ttop_tbot[0] = new TH2F("h0ratio_vs_ttop_tbot", "h0ratio_vs_ttop_tbot", 100, -3e-8, 3e-8, 100, 0, 25);
  hratio_vs_ttop_tbot[1] = new TH2F("h1ratio_vs_ttop_tbot", "h1ratio_vs_ttop_tbot", 100, -3e-8, 3e-8, 100, 0, 25);
  hratio_vs_ttop_tbot[1]->SetLineColor(kRed);


  for (int num = 0; num < 2; num++)
  {

    simple test(gname[num].c_str(), 0);

    if (isMC[num])
      test.SetMLUTable("./mlu_pbars_bottom_trap.dat");

    for (int i = 0; i < test.treeMCBV->GetEntries(); i++)
    {
      test.GetEntry(i);

      std::vector<double> cal_atop;
      std::vector<double> cal_abot;
      std::vector<double> cal_ttop_tbot;

      std::vector<int> cal_id;
      std::vector<int> mlu_id;

      for (int j = 0; j < test.bars_id->size(); j++)
      {

        if (isMC[num])
        {
          if(test.bars_id->at(j) == 33 || test.bars_id->at(j) == 31 || test.bars_id->at(j) == 30 || test.bars_id->at(j) == 29 || test.bars_id->at(j) == 28 )
            continue;
          
          double edep = test.bars_edep->at(j);
          double z = test.bars_z->at(j);

          double atop = edep * TMath::Exp(-(zmax - z)/lambda);
          double abot = edep * TMath::Exp(-(z - zmin)/lambda);

          double thit = test.bars_t->at(j);
          double ttop = thit + (zmax - z) / c_Ej + gRandom->Gaus(0, t_smearing); //preparare tuning barre saparate
          double tbot = thit + (z - zmin) / c_Ej + gRandom->Gaus(0, t_smearing);
          double ttop_tbot = ttop-tbot;

          // calibrate amplitudes
          // double Eres_top = deltaE_a * TMath::Exp(-deltaE_b*TMath::Log(atop));
          // double Eres_bot = deltaE_a * TMath::Exp(-deltaE_b*TMath::Log(abot));

          double Eres_top = TMath::Sqrt(qE_top + kE_top / atop);
          double Eres_bot = TMath::Sqrt(qE_bot + kE_bot / abot);

          // include non-linearity and resolution
          atop = top_gain / gain_c * (1 - TMath::Exp(-gain_c * atop)) * gRandom->Gaus(1., Eres_top);
          abot = bot_gain / gain_c * (1 - TMath::Exp(-gain_c * abot)) * gRandom->Gaus(1., Eres_bot);

          if (atop > thres * mlu_scale && abot > thres * mlu_scale)
          {
            mlu_id.push_back(test.bars_id->at(j));
          }
          if (atop > thres && abot > thres)
          {
            cal_atop.push_back(atop);
            cal_abot.push_back(abot);
            cal_ttop_tbot.push_back(ttop_tbot);
            
            cal_id.push_back(test.bars_id->at(j));
          }
        }
        else
        {
          if(test.bars_id->at(j) == 33 || test.bars_id->at(j) == 31 || test.bars_id->at(j) == 30 || test.bars_id->at(j) == 29 || test.bars_id->at(j) == 28 )
            continue;
          double atop = test.bars_atop->at(j);
          double abot = test.bars_abot->at(j);

          double ttop = test.bars_ttop->at(j);
          double tbot = test.bars_tbot->at(j);          
          double ttop_tbot = ttop-tbot;

          cal_atop.push_back(atop);
          cal_abot.push_back(abot);
          cal_ttop_tbot.push_back(ttop_tbot);
          
          cal_id.push_back(test.bars_id->at(j));
        }
      }

      if (isMC[num])
      {
        // apply mlu
        Bool_t BVMLU[16];
        for (unsigned int k = 0; k < 16; k++)
        {
          BVMLU[k] = false;
        }

        for (unsigned int k = 0; k < mlu_id.size(); k++)
        {
          unsigned int pad = (int)(mlu_id.at(k) / 4);
          BVMLU[pad] = 1;
        }
        uint16_t pattern = test.bool_array_to_uint16(BVMLU, 16);
        if (!test.mlu(pattern))
          continue;
      }

      for (int j = 0; j < cal_id.size(); j++)
      {
        hatop[num]->Fill(cal_atop.at(j));
        habot[num]->Fill(cal_abot.at(j));
        
        hlog[num]->Fill(-TMath::Log(cal_abot.at(j) / cal_atop.at(j)));

        hlog_vs_ttop_tbot[num]->Fill(-cal_ttop_tbot.at(j), -TMath::Log(cal_abot.at(j) / cal_atop.at(j)));
        hratio_vs_ttop_tbot[num]->Fill( -cal_ttop_tbot.at(j), cal_atop.at(j)/cal_abot.at(j));
    
        httop_tbot[num]->Fill(-cal_ttop_tbot.at(j));
      }

      hn[num]->Fill(cal_id.size());
    }
  }

  TLegend *legend = new TLegend(0.6745438,0.7575799,0.9372471,0.9300563,NULL,"brNDC");//(0.1,0.7,0.48,0.9);
  legend->SetTextSize(0.04668534);
  legend->AddEntry(hatop[1],"Data","l");
  legend->AddEntry(hatop[0],"Monte Carlo","l");
  
  gStyle->SetOptStat(0);

  TCanvas c("c","c",1800,1200);
  c.Divide(3, 2);

  c.cd(1);
  hatop[1]->DrawNormalized();
  hatop[0]->DrawNormalized("same");
  legend->Draw();

  c.cd(2);
  habot[1]->DrawNormalized();
  habot[0]->DrawNormalized("same");
  legend->Draw();

  c.cd(3);
  hn[1]->DrawNormalized();
  hn[0]->DrawNormalized("same");
  legend->Draw();

  c.cd(4);
  hlog[1]->DrawNormalized();
  hlog[0]->DrawNormalized("same");
  legend->Draw();

  c.cd(5);
  httop_tbot[1]->DrawNormalized();
  httop_tbot[0]->DrawNormalized("same");
  legend->Draw();

  c.SaveAs("status.png");

  TCanvas c2d("c2d","c2d",1000,1000);
  c2d.Divide(2,2);
  c2d.cd(1);
  hlog_vs_ttop_tbot[1]->Draw("colz");
  c2d.cd(2);
  hlog_vs_ttop_tbot[0]->Draw("colz");
  c2d.cd(3);
  hratio_vs_ttop_tbot[1]->Draw("colz");
  c2d.cd(4);
  hratio_vs_ttop_tbot[0]->Draw("colz");

  c2d.SaveAs("status2d.png");



  TH1F *href;
  TH1F *htest;

  href = hn[1];
  htest = hn[0];
  href->Scale(1. / href->GetSumOfWeights());
  htest->Scale(1. / htest->GetSumOfWeights());
  double err = err_scale / href->GetNbinsX();
  double pseudochi2 = 0;
  for (int i = 0; i < href->GetNbinsX() + 2; i++)
  {
    double res = href->GetBinContent(i) - htest->GetBinContent(i);
    pseudochi2 += (res * res / err / err) / href->GetNbinsX();
  }
  //  cout << pseudochi2/href->GetNbinsX() << endl;

  href = hlog[1];
  htest = hlog[0];
  href->Scale(1. / href->GetSumOfWeights());
  htest->Scale(1. / htest->GetSumOfWeights());
  err = err_scale / href->GetNbinsX();
  //  pseudochi2=0;
  for (int i = 0; i < href->GetNbinsX() + 2; i++)
  {
    double res = href->GetBinContent(i) - htest->GetBinContent(i);
    pseudochi2 += (res * res / err / err) / href->GetNbinsX();
  }
  //  cout << pseudochi2/ href->GetNbinsX()<< endl;

  href = hatop[1];
  htest = hatop[0];
  href->Scale(1. / href->GetSumOfWeights());
  htest->Scale(1. / htest->GetSumOfWeights());
  err = err_scale / href->GetNbinsX();
  //  pseudochi2=0;
  for (int i = 0; i < href->GetNbinsX() + 2; i++)
  {
    double res = href->GetBinContent(i) - htest->GetBinContent(i);
    pseudochi2 += (res * res / err / err) / href->GetNbinsX();
  }
  //  cout << pseudochi2/href->GetNbinsX() << endl;

  href = habot[1];
  htest = habot[0];
  href->Scale(1. / href->GetSumOfWeights());
  htest->Scale(1. / htest->GetSumOfWeights());
  err = err_scale / href->GetNbinsX();
  //  pseudochi2=0;
  for (int i = 0; i < href->GetNbinsX() + 2; i++)
  {
    double res = href->GetBinContent(i) - htest->GetBinContent(i);
    pseudochi2 += (res * res / err / err) / href->GetNbinsX();
  }
  //  cout << pseudochi2/href->GetNbinsX() << endl;

  href = httop_tbot[1];
  htest = httop_tbot[0];
  href->Scale(1. / href->GetSumOfWeights());
  htest->Scale(1. / htest->GetSumOfWeights());
  err = err_scale / href->GetNbinsX();
  //  pseudochi2=0;
  for (int i = 0; i < href->GetNbinsX() + 2; i++)
  {
    double res = href->GetBinContent(i) - htest->GetBinContent(i);
    pseudochi2 += (res * res / err / err) / href->GetNbinsX();
  }
  //  cout << pseudochi2/href->GetNbinsX() << endl;

  return pseudochi2;
}
