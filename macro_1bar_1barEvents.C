#include "simple.C"
#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <fstream>
using namespace std;

#define NBARS 64

double macro_1bar_1barEvents(int barnumber=0, int filenumber = 0)
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

  //gStyle->SetOptStat(0);

  double mlu_scale =0.5;
  double c_Ej = 158559135.35984075;
  double t_smearing = 800e-12;

  double top_gain[NBARS];
  //double qE_top[NBARS];
  //double deltaE_top[NBARS];
  double bot_gain[NBARS];
  double qE[NBARS];
  double kE[NBARS];
  double gain_c[NBARS];
  double lambda[NBARS];

  //reading configuration files for parameters
  for(int i=0; i<NBARS; i++)
  {
    ifstream iFile(Form("ConfigFiles_1barEvents/File_%d.txt",i), ios::in);
    if (!iFile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    int cont=0;
    double p;
    while (iFile >> p) 
    {
      if (cont == 0)      top_gain[i] = p;
      // else if (cont == 1) qE_top[i] = p;
      // else if (cont == 2) deltaE_top[i] = p;
      else if (cont == 1) bot_gain[i] = p;
      else if (cont == 2) qE[i] = p;
      else if (cont == 3) kE[i] = p;
      else if (cont == 4) gain_c[i] = p;
      else if (cont == 5) lambda[i] = p;  

      cont++;
    }
    iFile.close(); 
  }


  // if (ParN == 0)  top_gain[BarN] = ParVal;
  // else if (ParN == 1) top_spread[BarN] = ParVal;
  // else if (ParN == 2) bot_gain[BarN] = ParVal;
  // else if (ParN == 3) bot_spread[BarN] = ParVal;
  // else if (ParN == 4) deltaE_a[BarN] = ParVal;
  // else if (ParN == 5) gain_c[BarN] = ParVal;
  // else if (ParN == 6) lambda[BarN] = ParVal;
  // else {cout<<"Error ParN out of range!!!"<<endl; return 0.0;}
  

  // for(int i=0; i<NBARS; i++)
  // {
  //   cout<<i<<")_______"<<endl;
  //   cout<< top_gain[i] <<endl;
  //   cout<< top_spread[i] <<endl;
  //   cout<< bot_gain[i] <<endl;
  //   cout<< bot_spread[i] <<endl;
  //   cout<< deltaE_a[i] <<endl;
  //   cout<< gain_c[i] <<endl;
  //   cout<< lambda[i] <<endl;
  // }

  const double thres = threshold/TMath::Power(2,14);//0.0366;
  const double zmin = -1.3;
  const double zmax = 1.3;

  // double gain_top[NBARS];
  // double gain_bot[NBARS];

  // for (int i = 0; i < NBARS; i++)
  // {
  //   gain_top[i] = top_gain[i] + gRandom->Gaus(0., top_gain[i] * top_spread[i]);
  //   gain_bot[i] = bot_gain[i] + gRandom->Gaus(0., bot_gain[i] * bot_spread[i]);
  // }

  // arXiv: 1905.06032 SiPM + Scint non linearity
  // https://ieeexplore.ieee.org/document/5874113 resolution
  // deltaE/E = a E^-b
  // <Q> = gain/c * ( 1 - exp(-cE))

  double err_scale = 0.2;

  string gname[2] = {"MC_BV_TREE_cosmics_ecomug_hsphere_1M_B1T_set999.root", (string)filename};
  bool isMC[2] = {true, false};

  TH1F *hatop[2][NBARS]; // first index 0->MC, 1->Data; second index bar number
  TH1F *habot[2][NBARS];
  TH1F *hlog[2][NBARS];

  for (int i = 0; i < NBARS; i++)
  {
    hatop[1][i] = new TH1F(Form("hatop1_bar%d", i), Form("hatop1_bar%d", i), 50, 0, 1);
    hatop[1][i]->SetLineColor(kRed);
    habot[1][i] = new TH1F(Form("habot1_bar%d", i), Form("habot1_bar%d", i), 50, 0, 1);
    habot[1][i]->SetLineColor(kRed);
    hlog[1][i] = new TH1F(Form("hlog1_bar%d", i), Form("hlog1_bar%d", i), 50, -4, 4);
    hlog[1][i]->SetLineColor(kRed);

    hatop[0][i] = new TH1F(Form("hatop0_bar%d", i), Form("hatop0_bar%d", i), 50, 0, 1);
    habot[0][i] = new TH1F(Form("habot0_bar%d", i), Form("habot0_bar%d", i), 50, 0, 1);
    hlog[0][i] = new TH1F(Form("hlog0_bar%d", i), Form("hlog0_bar%d", i), 50, -4, 4);
  }

  TH1F *httop_tbot[2][NBARS];
  for (int i = 0; i < NBARS; i++)
  {
    httop_tbot[0][i] = new TH1F(Form("h0ttop_tbot_bar%d", i), Form("h0ttop_tbot_bar%d", i), 50, -3e-8, 3e-8);
    httop_tbot[1][i] = new TH1F(Form("h1ttop_tbot_bar%d", i), Form("h1ttop_tbot_bar%d", i), 50, -3e-8, 3e-8);
    httop_tbot[1][i]->SetLineColor(kRed);
  }

  TH2F *hlog_vs_ttop_tbot[2][NBARS];
  for (int i = 0; i < NBARS; i++)
  {
    hlog_vs_ttop_tbot[0][i] = new TH2F(Form("h0log_vs_ttop_tbot%d", i), Form("h0log_vs_ttop_tbot%d", i), 100, -3e-8, 3e-8, 100, -4, 4);
    hlog_vs_ttop_tbot[1][i] = new TH2F(Form("h1log_vs_ttop_tbot%d", i), Form("h1log_vs_ttop_tbot%d", i), 100, -3e-8, 3e-8, 100, -4, 4);
    hlog_vs_ttop_tbot[1][i]->SetLineColor(kRed);
  }

  TH2F *hratio_vs_ttop_tbot[2][NBARS];
  for (int i = 0; i < NBARS; i++)
  {
    hratio_vs_ttop_tbot[0][i] = new TH2F(Form("h0ratio_vs_ttop_tbot%d", i), Form("h0ratio_vs_ttop_tbot%d", i), 100, -3e-8, 3e-8, 100, 0, 25);
    hratio_vs_ttop_tbot[1][i] = new TH2F(Form("h1ratio_vs_ttop_tbot%d", i), Form("h1ratio_vs_ttop_tbot%d", i), 100, -3e-8, 3e-8, 100, 0, 25);
    hratio_vs_ttop_tbot[1][i]->SetLineColor(kRed);
  }

  TH1F *hn[2];
  hn[0] = new TH1F("hn0", "hn0", 12, 1, 13);
  hn[1] = new TH1F("hn1", "hn1", 12, 1, 13);
  hn[1]->SetLineColor(kRed);

  for (int num = 0; num < 2; num++)
  {

    simple test(gname[num], 0);

    if (isMC[num])
      test.SetMLUTable("./mlu_pbars_bottom_trap.dat");

    Int_t nRep;
    if (isMC[num])
      nRep = 8;
    else
      nRep = 1;

    double rnd_rot = gRandom->Rndm();
    int rot = (int)(rnd_rot * 64);
    Int_t Nentries = test.treeMCBV->GetEntries();
    for (int i = 0; i < Nentries * nRep; i++) //loop on events
    {

      test.GetEntry(i % Nentries);
      // ruota indici barre colpite di i*64/8
      //  oppure casulamente con i in [1,8]

      if (i % Nentries == 0)
      {
        rnd_rot = gRandom->Rndm();
        rot = (int)(rnd_rot * 64);
      }

      std::vector<int> cal_id; //(cal_id+i*64/nRep)%64
      std::vector<int> mlu_id;

      std::vector<double> cal_atop[NBARS];
      std::vector<double> cal_abot[NBARS];
      std::vector<double> cal_ttop_tbot[NBARS];

      std::vector<int> bars_id_rot;

      for (int j = 0; j < test.bars_id->size(); j++) // create a new rotated vector of bars
      {
        if (!isMC[num])
          bars_id_rot.push_back(test.bars_id->at(j));
        else
        {
          int new_id = (test.bars_id->at(j) + rot) % 64;
          bars_id_rot.push_back(new_id);
        }
      }
      
      if(bars_id_rot.size() != 1) continue; //-> selection of just 1 bar events
      //non si può fare in modo corretto perchè la molteplicità dell'evento non è una questione solo geometrica, 
      //ma dipende anche  dalla calibrazione delle altre barre. 
      //Un mu può colpire più di una barra ma produrre un segnale che non è sopra soglia. 
      //Mentre nel MC senza calibrazione prendiamo tutto
      
      for (int j = 0; j < bars_id_rot.size(); j++)
      {
        if (bars_id_rot.at(j) >= NBARS)
          continue;

        if(bars_id_rot.at(j)!=barnumber) continue;

        if (isMC[num])
        {
          double edep = test.bars_edep->at(j);
          double z = test.bars_z->at(j);

          double atop = edep * TMath::Exp(-(zmax - z) / lambda[bars_id_rot.at(j)]);
          double abot = edep * TMath::Exp(-(z - zmin) / lambda[bars_id_rot.at(j)]);

          double thit = test.bars_t->at(j);
          double ttop = thit + (zmax - z) / c_Ej + gRandom->Gaus(0, t_smearing);
          double tbot = thit + (z - zmin) / c_Ej + gRandom->Gaus(0, t_smearing);
          double ttop_tbot = ttop - tbot;

          // calibrate amplitudes
          // double Eres_top = deltaE_a * TMath::Exp(-deltaE_b*TMath::Log(atop));
          // double Eres_bot = deltaE_a * TMath::Exp(-deltaE_b*TMath::Log(abot));

          double Eres_top = TMath::Sqrt(qE[bars_id_rot.at(j)] + kE[bars_id_rot.at(j)] / atop);
          double Eres_bot = TMath::Sqrt(qE[bars_id_rot.at(j)] + kE[bars_id_rot.at(j)] / abot);

          // include non-linearity and resolution
          atop = top_gain[bars_id_rot.at(j)] / gain_c[bars_id_rot.at(j)] * (1 - TMath::Exp(-gain_c[bars_id_rot.at(j)] * atop)) * gRandom->Gaus(1., Eres_top);
          abot = bot_gain[bars_id_rot.at(j)] / gain_c[bars_id_rot.at(j)] * (1 - TMath::Exp(-gain_c[bars_id_rot.at(j)] * abot)) * gRandom->Gaus(1., Eres_bot);

          if (atop > thres * mlu_scale && abot > thres * mlu_scale)
          {
            mlu_id.push_back(bars_id_rot.at(j));
          }
          if (atop > thres && abot > thres)
          {
            cal_atop[bars_id_rot.at(j)].push_back(atop);
            cal_abot[bars_id_rot.at(j)].push_back(abot);

            cal_ttop_tbot[bars_id_rot.at(j)].push_back(ttop_tbot);

            cal_id.push_back(bars_id_rot.at(j));
          }
        }
        else
        {
          double atop = test.bars_atop->at(j);
          double abot = test.bars_abot->at(j);

          double ttop = test.bars_ttop->at(j);
          double tbot = test.bars_tbot->at(j);
          double ttop_tbot = ttop - tbot;

          cal_atop[bars_id_rot.at(j)].push_back(atop);
          cal_abot[bars_id_rot.at(j)].push_back(abot);
          cal_ttop_tbot[bars_id_rot.at(j)].push_back(ttop_tbot);

          cal_id.push_back(bars_id_rot.at(j));
        }
      }

      // if (isMC[num]) // non funziona
      // {
      //   // apply mlu
      //   Bool_t BVMLU[16];
      //   for (unsigned int k = 0; k < 16; k++)
      //   {
      //     BVMLU[k] = false;
      //   }

      //   for (unsigned int k = 0; k < mlu_id.size(); k++)
      //   {
      //     unsigned int pad = (int)(mlu_id.at(k) / 4);
      //     BVMLU[pad] = 1;
      //   }
      //   uint16_t pattern = test.bool_array_to_uint16(BVMLU, 16);
      //   if (!test.mlu(pattern))
      //     continue;
      // }

      for (int k = 0; k < NBARS; k++)
      {
        for (int j = 0; j < cal_ttop_tbot[k].size(); j++)
        {
          httop_tbot[num][k]->Fill(-cal_ttop_tbot[k].at(j));
          hatop[num][k]->Fill(cal_atop[k].at(j));
          habot[num][k]->Fill(cal_abot[k].at(j));

          hlog[num][k]->Fill(-TMath::Log(cal_abot[k].at(j) / cal_atop[k].at(j)));

          hlog_vs_ttop_tbot[num][k]->Fill(-cal_ttop_tbot[k].at(j), -TMath::Log(cal_abot[k].at(j) / cal_atop[k].at(j)));
          hratio_vs_ttop_tbot[num][k]->Fill( -cal_ttop_tbot[k].at(j),  cal_atop[k].at(j)/cal_abot[k].at(j));
        }
      }

      hn[num]->Fill(cal_id.size());
    }
  }

  TCanvas c[NBARS];
  TCanvas c2d[NBARS];
  for (int k = 0; k < NBARS; k++)
  {
    if(k!=barnumber) continue;
    
    c[k].Divide(3, 2);

    c[k].cd(1);
    hatop[1][k]->DrawNormalized();
    hatop[0][k]->DrawNormalized("same");

    c[k].cd(2);
    habot[1][k]->DrawNormalized();
    habot[0][k]->DrawNormalized("same");

    c[k].cd(3);
    hn[1]->DrawNormalized();
    hn[0]->DrawNormalized("same");

    c[k].cd(4);
    hlog[1][k]->DrawNormalized();
    hlog[0][k]->DrawNormalized("same");

    c[k].cd(5);
    httop_tbot[1][k]->DrawNormalized();
    httop_tbot[0][k]->DrawNormalized("same");

    c[k].SaveAs(Form("Plots_Multbars_1barEvents/status_bar%d.png", k));
    
    c2d[k].Divide(2,2);
    c2d[k].cd(1);
    hlog_vs_ttop_tbot[1][k]->Draw("colz");
    c2d[k].cd(2);
    hlog_vs_ttop_tbot[0][k]->Draw("colz");
    c2d[k].cd(3);
    hratio_vs_ttop_tbot[1][k]->Draw("colz");
    c2d[k].cd(4);
    hratio_vs_ttop_tbot[0][k]->Draw("colz");

    c2d[k].SaveAs(Form("Plots_Multbars_1barEvents/status2d_bar%d.png", k));

  }

  TH1F *href;
  TH1F *htest;

  double err;
  double pseudochi2 = 0;

  // href = hn[1];
  // htest = hn[0];
  // href->Scale(1. / href->GetSumOfWeights());
  // htest->Scale(1. / htest->GetSumOfWeights());
  // err = err_scale / href->GetNbinsX();
  // for (int i = 0; i < href->GetNbinsX() + 2; i++)
  // {
  //   double res = href->GetBinContent(i) - htest->GetBinContent(i);
  //   pseudochi2 += (res * res / err / err) / href->GetNbinsX();
  // }
  // //  cout << pseudochi2/href->GetNbinsX() << endl;

  for (int k = 0; k < NBARS; k++)
  {
    if(k!=barnumber) continue;

    href = hlog[1][k];
    htest = hlog[0][k];
    href->Scale(1. / href->GetSumOfWeights());
    htest->Scale(1. / htest->GetSumOfWeights());
    err = err_scale / href->GetNbinsX();
    //  pseudochi2=0;
    for (int i = 0; i < href->GetNbinsX() + 2; i++)
    {
      double res = href->GetBinContent(i) - htest->GetBinContent(i);
      pseudochi2 += (res * res / err / err) / href->GetNbinsX();
    }
    // cout << pseudochi2/ href->GetNbinsX()<< endl;

    href = hatop[1][k];
    htest = hatop[0][k];
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

    href = habot[1][k];
    htest = habot[0][k];
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

    // href = httop_tbot[1][k];
    // htest = httop_tbot[0][k];
    // href->Scale(1. / href->GetSumOfWeights());
    // htest->Scale(1. / htest->GetSumOfWeights());
    // err = err_scale / href->GetNbinsX();
    // //  pseudochi2=0;
    // for (int i = 0; i < href->GetNbinsX() + 2; i++)
    // {
    //   double res = href->GetBinContent(i) - htest->GetBinContent(i);
    //   pseudochi2 += (res * res / err / err) / href->GetNbinsX();
    // }
    //  cout << pseudochi2/href->GetNbinsX() << endl;
  }
  return pseudochi2;
}
