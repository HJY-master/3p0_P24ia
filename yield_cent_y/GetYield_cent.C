#include "TH1F.h"
#include "TH2F.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAttPad.h"
#include "stdio.h"
#include "TLegend.h"
#include "TFile.h"
#include "TCanvas.h"

#include "/Users/junyihan/workspace/myfunc/func.h"
#include "/Users/junyihan/workspace/myfunc/drawfunc.h"

void GetYield_cent()
{
  double yscale = 1.2;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadTopMargin(0.05);
  
  TFile *pol_file   = TFile::Open("../outsig.root");
  TFile *polBg_file = TFile::Open("../outbkg.root");
  //============================================================================
  TH1F *ax;
  
  const int nCent = 7;
  const int nMode = 10;
    
    float nsigma = 3;
  
  TString centName[nCent] = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "5-40%", "10-40%"};
  TString yName[nMode];
    TString ptName[nMode];
    for(int imode = 0; imode < nMode; imode++)
    {
        yName[imode] = Form("%.1f < y < %.1f", -1.0+0.1*imode, -1.0+0.1*(imode+1));
        ptName[imode] = "0.2 < p_{T} < 3.0 GeV/c";
    }
//  TString ptName[nMode] = {"0.4 < p_{T} < 0.8 GeV/c", "0.2 < p_{T} < 2.0 GeV/c", "0.2 < p_{T} < 1.0 GeV/c", "0.2 < p_{T} < 2.6 GeV/c"};
  double yMin[nMode];
  double yMax[nMode];
    double ptMin[nMode];
    double ptMax[nMode];
    for(int imode = 0; imode < nMode; imode++)
    {
        yMin[imode] = -1.0+0.1*imode;
        yMax[imode] = -1.0+0.1*(imode+1);
        ptMin[imode] = 0.2;
        ptMax[imode] = 3.0;
    }
//  double ptMin[nMode] = {0.4, 0.2, 0.2, 0.2};
//  double ptMax[nMode] = {0.8, 2.0, 1.0, 2.6};

  double nXiYield[nMode][nCent] = {0.};
  double nXiYieldErr[nMode][nCent] = {0.};
  //============================================================================
  // get input
  TH3F *h3InvMassXiYvsPt_tot[nCent-2];
  TH3F *h3InvMassXiYvsPt_bg[nCent-2];

  for(int icent=0; icent<nCent-2; icent++){
    h3InvMassXiYvsPt_tot[icent] = (TH3F*)pol_file->Get(Form("hMassLdRapidityvsPt_cent%d", icent));
    h3InvMassXiYvsPt_bg[icent] = (TH3F*)polBg_file->Get(Form("hMassLdRapidityvsPt_cent%d", icent));
  }

  TH1F *invMassVsCent_tot[nMode][nCent];
  TH1F *invMassVsCent_bg[nMode][nCent];
  
  int ptBinIdx_Low, ptBinIdx_High;
  int yBinIdx_Low, yBinIdx_High;

  for(int icent = 0; icent < nCent-2; icent++)
  {
    for(int imode = 0; imode < nMode; imode++)
    {
      auto ybinl = ((TH1F*)h3InvMassXiYvsPt_tot[icent]->ProjectionY())->FindBin(yMin[imode]+0.00001);
      auto ybinr = ((TH1F*)h3InvMassXiYvsPt_tot[icent]->ProjectionY())->FindBin(yMax[imode]-0.00001);
      auto pbinl = ((TH1F*)h3InvMassXiYvsPt_tot[icent]->ProjectionZ())->FindBin(ptMin[imode]+0.00001);
      auto pbinr = ((TH1F*)h3InvMassXiYvsPt_tot[icent]->ProjectionZ())->FindBin(ptMax[imode]-0.00001);

      invMassVsCent_tot[imode][icent] = (TH1F*)h3InvMassXiYvsPt_tot[icent]->ProjectionX(Form("tot_mode%d_cent%d", imode, icent), ybinl, ybinr, pbinl, pbinr);
      invMassVsCent_bg[imode][icent] = (TH1F*)h3InvMassXiYvsPt_bg[icent]->ProjectionX(Form("bg_mode%d_cent%d", imode, icent), ybinl, ybinr, pbinl, pbinr);
    
      invMassVsCent_tot[imode][icent]->Rebin(2);
      invMassVsCent_bg[imode][icent]->Rebin(2);
    }
  }


  // TH1F *invMassVsCent_tot_merge[nMode];
  // TH1F *invMassVsCent_bg_merge[nMode];
    
    // 5-40%
    for(int imode=0; imode<nMode; imode++){
        invMassVsCent_tot[imode][nCent-2] = (TH1F*)invMassVsCent_tot[imode][1]->Clone(Form("sig_mode%d_cent%d", imode, nCent-2));
        invMassVsCent_bg[imode][nCent-2] = (TH1F*)invMassVsCent_bg[imode][1]->Clone(Form("bkg_mode%d_cent%d", imode, nCent-2));
        for(int icent=2; icent<nCent-2; icent++){
          invMassVsCent_tot[imode][nCent-2]->Add(invMassVsCent_tot[imode][icent]);
          invMassVsCent_bg[imode][nCent-2]->Add(invMassVsCent_bg[imode][icent]);
        }
    }
    
    // 10-40%
    for(int imode=0; imode<nMode; imode++){
        invMassVsCent_tot[imode][nCent-1] = (TH1F*)invMassVsCent_tot[imode][2]->Clone(Form("sig_mode%d_cent%d", imode, nCent-1));
        invMassVsCent_bg[imode][nCent-1] = (TH1F*)invMassVsCent_bg[imode][2]->Clone(Form("bkg_mode%d_cent%d", imode, nCent-1));
        for(int icent=3; icent<nCent-2; icent++){
          invMassVsCent_tot[imode][nCent-1]->Add(invMassVsCent_tot[imode][icent]);
          invMassVsCent_bg[imode][nCent-1]->Add(invMassVsCent_bg[imode][icent]);
        }
    }

  //============================================================================
  
  TF1 *f1_cent;
  TF1 *f2_cent;
  TF1 *f1_centFit;
  TF1 *f2_centFit;
  
  TF1 *f1_centRotation;
  TF1 *f2_centRotation;
  TF1 *f2_centSig;
  
  double IMcut_3Sigma[nMode][nCent][2];
  double IMcut_SigmaErr[nMode][nCent];
  double IMcut_mean[nMode][nCent];
  double IMcut_meanErr[nMode][nCent];
  
  TH1F *invMassVsCent_sig[nMode][nCent];
  TF1  *fInvMassVsCent_tot[nMode][nCent];
  TH1F *invMassVsCentReal_sig[nMode][nCent];
  
  TF1 *fInvMassVsCentReal_bg[nMode][nCent];
  TF1 *fInvMassVsCentReal_sig[nMode][nCent];
  
  TH1F *invMassVsCentFit_sig[nMode][nCent];
  TF1 *fInvMassVsCentFit_bg[nMode][nCent];
  
  TH1F *invMassVsCentFitReal_sig[nMode][nCent];
  TF1 *fInvMassVsCentFitReal_bg[nMode][nCent];

  double lrange = 1.09, rrange = 1.15;
  ofstream fout[nCent];
//  fout.open("IMcut.txt");
    for(int icent=0; icent<nCent; icent++){
        fout[icent].open(Form("IMcut_cent%s.txt", centName[icent].Data()));
    }
  for(int imode=0; imode<nMode; imode++){
    for(int icent=0; icent<nCent; icent++){
      double nInclLambda_l = invMassVsCent_tot[imode][icent]->Integral(invMassVsCent_tot[imode][icent]->FindBin(1.09), invMassVsCent_tot[imode][icent]->FindBin(1.105));
      double nBgLambda_l   = invMassVsCent_bg[imode][icent]->Integral(invMassVsCent_bg[imode][icent]->FindBin(1.09), invMassVsCent_bg[imode][icent]->FindBin(1.105));
      double nInclLambda_r = invMassVsCent_tot[imode][icent]->Integral(invMassVsCent_tot[imode][icent]->FindBin(1.13), invMassVsCent_tot[imode][icent]->FindBin(1.15));
      double nBgLambda_r   = invMassVsCent_bg[imode][icent]->Integral(invMassVsCent_bg[imode][icent]->FindBin(1.13), invMassVsCent_bg[imode][icent]->FindBin(1.15));
      double scaler = 0;
      if((nBgLambda_l + nBgLambda_r) != 0 ) scaler = (nInclLambda_l + nInclLambda_r) / (nBgLambda_l + nBgLambda_r);
      if(icent == nCent-1) cout << "scaler: " << scaler << endl;
      invMassVsCent_bg[imode][icent]->Scale(scaler);
      
      //Rotation method
      invMassVsCent_sig[imode][icent] = SubtractBG( invMassVsCent_tot[imode][icent], invMassVsCent_bg[imode][icent] );
      
      double	MeanMin	= 1.1112;
      double	MeanMax	= 1.1202;
      double	WidthMin= 0.0010;
      double	WidthMax= 0.0025;

      f1_centRotation = new TF1("f1_centRotation", Ploy_Student_t, lrange, rrange, 6);
      f1_centRotation->SetParNames("nu", "m0", "Width", "norm", "p0", "p1");
      f1_centRotation->SetParameter(0, 1);
      f1_centRotation->SetParameter(1, 1.11600e+00);
      f1_centRotation->SetParLimits(1, MeanMin, MeanMax);
      f1_centRotation->SetParameter(2, 1.50000e-03);
      f1_centRotation->SetParLimits(2, WidthMin, WidthMax); 
      f1_centRotation->SetParLimits(3, 0, 1e8);
      f1_centRotation->SetParameter(4, 1);
      f1_centRotation->SetParameter(5, 1);// from zuowen
      
      invMassVsCent_sig[imode][icent]->Fit("f1_centRotation", "RNQM", "", lrange, rrange);
      f1_centRotation->SetParameters( f1_centRotation->GetParameter("nu"),
                                      f1_centRotation->GetParameter("m0"),
                                      f1_centRotation->GetParameter("Width"),
                                      f1_centRotation->GetParameter("norm"),
                                      f1_centRotation->GetParameter("p0"),
                                      f1_centRotation->GetParameter("p1")
                                      );
      invMassVsCent_sig[imode][icent]->Fit("f1_centRotation", "RNQM", "", lrange, rrange);
      
      double *parR;
      const double *parErrR;
      parR = f1_centRotation->GetParameters();
      parErrR = f1_centRotation->GetParErrors();
      
      fInvMassVsCent_tot[imode][icent] = (TF1*)f1_centRotation;
      f2_centRotation = new TF1("f2_centRotation", background, lrange, rrange, 2); //fit the residual background by linear function
      f2_centRotation->SetLineColor(kBlue);
      f2_centRotation->SetLineStyle(2);
      f2_centRotation->SetParameters(&parR[4]);
      f2_centRotation->SetParErrors(&parErrR[4]);
      fInvMassVsCentReal_bg[imode][icent] = (TF1*) f2_centRotation;

      f2_centSig = new TF1("f2_centSig", Student_t, lrange, rrange, 3);
      f2_centSig->SetLineColor(kRed);
      f2_centSig->SetLineStyle(1);
      f2_centSig->SetParameters(&parR[0]);
      f2_centSig->SetParErrors(&parErrR[0]);
      fInvMassVsCentReal_sig[imode][icent] = (TF1*) f2_centSig;
      
      invMassVsCentReal_sig[imode][icent] = SubtractBGFn(invMassVsCent_sig[imode][icent], f2_centRotation, lrange, rrange);
      
      IMcut_3Sigma[imode][icent][0] = f1_centRotation->GetParameter("m0") - nsigma * f1_centRotation->GetParameter("Width");
      IMcut_3Sigma[imode][icent][1] = f1_centRotation->GetParameter("m0") + nsigma * f1_centRotation->GetParameter("Width");

      IMcut_SigmaErr[imode][icent] = parErrR[2];

      IMcut_mean[imode][icent] = f1_centRotation->GetParameter("m0");
      IMcut_meanErr[imode][icent] = f1_centRotation->GetParError(1);

      if(icent == 10) cout << "left boundary: " << IMcut_3Sigma[imode][icent][0] << " and right boundary: " << IMcut_3Sigma[imode][icent][1] << endl;
      if(icent == nCent-1) cout << "5-40: 2sigma left boundary: " << f1_centRotation->GetParameter("m0") - 2.0 * f1_centRotation->GetParameter("Width") << " and right boundary: " << f1_centRotation->GetParameter("m0") + 2.0 * f1_centRotation->GetParameter("Width") << endl;
//      if(icent == 10) fout << IMcut_3Sigma[imode][icent][0] << " " << IMcut_3Sigma[imode][icent][1] << endl;
        fout[icent] << IMcut_3Sigma[imode][icent][0] << " " << IMcut_3Sigma[imode][icent][1] << endl;
    }
  }
  
  //============================================================================
  TH1F *hCountsCent[nMode];
  for(int imode=0; imode<nMode; imode++){
    hCountsCent[imode] = new TH1F(Form("hCountsCent_mode%d", imode), "", nCent, 0, nCent);
  }

  TCanvas *ca_invMass = new TCanvas("ca_invMass", "invMass of Lambda", 1500, 900);
  TPad *p0_invMass = new TPad("p0_invMass", "invMass of Lambda", 0.05, 0.07, 0.47, 0.95);
  TPad *p1_invMass = new TPad("p1_invMass", "invMass of Lambda", 0.52, 0.07, 0.97, 0.95);
  TLatex *text_invm = new TLatex(.5, .2, "");
  draw_21pad(ca_invMass, p0_invMass, p1_invMass, text_invm, "M_{inv} [GeV/c^{2}]", "Counts");
  p0_invMass->SetMargin(0.1, 0.03, 0.045, 0.03);
  p1_invMass->SetMargin(0.1, 0.03, 0.045, 0.03);

  TFile *invmass_ncent[nMode];
  double labeloff = 1.092;
  
  for(int imode=0; imode<nMode; imode++){
    invmass_ncent[imode] = new TFile(Form("invmass_mode%d.root", imode), "RECREATE");
    invmass_ncent[imode]->cd();
    for(int icent=0; icent<nCent; icent++){
      int lowerMassBin_Xi = invMassVsCent_tot[imode][icent]->FindBin(IMcut_3Sigma[imode][icent][0]) + 1;
      int upperMassBin_Xi = invMassVsCent_tot[imode][icent]->FindBin(IMcut_3Sigma[imode][icent][1]);
      float nTotal            = invMassVsCent_tot[imode][icent]->Integral(lowerMassBin_Xi, upperMassBin_Xi);
      
      double err = 0.;
      double nSignal = invMassVsCentReal_sig[imode][icent]->IntegralAndError(lowerMassBin_Xi, upperMassBin_Xi, err );
      
      float nBackground  = nTotal - nSignal;
      //========================================
      double chi2 = fInvMassVsCent_tot[imode][icent]->GetChisquare();
      double ndf = fInvMassVsCent_tot[imode][icent]->GetNDF();

      nXiYield[imode][icent] = nSignal;
      nXiYieldErr[imode][icent] = err;
      if(nSignal>0){
        hCountsCent[imode]->SetBinContent(icent+1, nSignal);
        hCountsCent[imode]->SetBinError(icent+1, err);
      }
      //========================================

      p0_invMass->cd();
      float maxCounts       = invMassVsCent_tot[imode][icent]->GetBinContent(invMassVsCent_tot[imode][icent]->GetMaximumBin());
      float maxCountsErr    = invMassVsCent_tot[imode][icent]->GetBinError(invMassVsCent_tot[imode][icent]->GetMaximumBin());

      ax = gPad->DrawFrame( lrange, -0.3*(maxCounts + maxCountsErr), rrange, yscale*(maxCounts + maxCountsErr) ); // no EPcorr

      SetAxis(ax, 1.3, 1.3);
      ax->SetLabelSize(22,"X");
      ax->SetLabelSize(22,"Y");
      ax->Draw();
      
      invMassVsCent_tot[imode][icent]->Draw("ESAME");
      invMassVsCent_tot[imode][icent]->SetMarkerStyle(20);
      invMassVsCent_tot[imode][icent]->SetMarkerColor(kRed);
      
      invMassVsCent_bg[imode][icent]->Draw("ESAME");
      invMassVsCent_bg[imode][icent]->SetMarkerStyle(24);
      invMassVsCent_bg[imode][icent]->SetMarkerColor(kBlue);
      
      TLine *lineLFit = new TLine(IMcut_3Sigma[imode][icent][0], 0., IMcut_3Sigma[imode][icent][0], yscale*(maxCounts + maxCountsErr));
      lineLFit->SetLineColor(kRed);
      lineLFit->SetLineStyle(2);
      lineLFit->SetLineWidth(2);
      lineLFit->Draw("SAME");
      
      TLine *lineRFit = new TLine(IMcut_3Sigma[imode][icent][1], 0., IMcut_3Sigma[imode][icent][1], yscale*(maxCounts + maxCountsErr));
      lineRFit->SetLineColor(kRed);
      lineRFit->SetLineStyle(2);
      lineRFit->SetLineWidth(2);
      lineRFit->Draw("SAME");
      
      maxCounts = yscale*(maxCounts + maxCountsErr);

      TLatex *tex_tip = new TLatex(1.108, 0.9*maxCounts, Form("%.1f#sigma invmass region", nsigma));
      tex_tip->SetTextFont(42);
      tex_tip->SetTextSize(0.03);
      tex_tip->Draw("same");
      
      TLatex *tex_comb = new TLatex(1.126, 0.5*maxCounts, "(rotate #pi^{-} 180#circ 1 time)");
      tex_comb->SetTextFont(42);
      tex_comb->SetTextSize(0.032);
      tex_comb->Draw("same");

      TLatex *tex_decay = new TLatex(labeloff, 0.85*maxCounts, "#Lambda #Rightarrow p + #pi^{-}");
      tex_decay->SetTextFont(42);
      tex_decay->SetTextSize(0.04);
      tex_decay->Draw("same");

      TLatex *tex_sideband = new TLatex(labeloff, -0.13*maxCounts, "Normalize range: [1.09, 1.105] & [1.13, 1.15]");
      tex_sideband->SetTextFont(42);
      tex_sideband->SetTextSize(0.032);
      tex_sideband->Draw("same");

      TLegend *legT = new TLegend(0.62, 0.65, 0.9, 0.92);
      legT->SetFillColor(10);
      // legT->SetFillStyle(0);
      legT->SetLineStyle(3004);
      legT->SetLineColor(10);
      legT->SetLineWidth(0.);
      legT->SetTextFont(42);
      legT->SetTextSize(0.04);
      legT->SetHeader(Form("#Lambda %s Au+Au", centName[icent].Data()));
      legT->AddEntry(invMassVsCent_tot[imode][icent], "total", "pe");
      legT->AddEntry(invMassVsCent_bg[imode][icent], "comb. BG", "pe");
      legT->Draw();

      p1_invMass->cd();
      maxCounts       = invMassVsCent_sig[imode][icent]->GetBinContent(invMassVsCent_sig[imode][icent]->GetMaximumBin());
      maxCountsErr    = invMassVsCent_sig[imode][icent]->GetBinError(invMassVsCent_sig[imode][icent]->GetMaximumBin());

      yscale = 1.6;
      ax = gPad->DrawFrame( lrange, -0.4*(maxCounts + maxCountsErr), rrange, yscale*(maxCounts + maxCountsErr) ); // no EPcorr
      SetAxis(ax, 1.3, 1.3);
      ax->SetLabelSize(22,"X");
      ax->SetLabelSize(22,"Y");
      ax->Draw();

      invMassVsCent_sig[imode][icent]->SetLineColor(kBlack);
      invMassVsCent_sig[imode][icent]->SetMarkerColor(kBlack);
      invMassVsCent_sig[imode][icent]->SetMarkerStyle(20);
      invMassVsCent_sig[imode][icent]->Draw("ESAME");

      fInvMassVsCent_tot[imode][icent]->SetMarkerColor(kBlack);
      fInvMassVsCent_tot[imode][icent]->SetLineColor(kBlack);
      fInvMassVsCent_tot[imode][icent]->Draw("lSAME");

      invMassVsCentReal_sig[imode][icent]->Draw("SAME HIST");
      invMassVsCentReal_sig[imode][icent]->SetLineColor(kRed);
      invMassVsCentReal_sig[imode][icent]->SetFillColor(kRed);
      invMassVsCentReal_sig[imode][icent]->SetFillStyle(3004);

      fInvMassVsCentReal_bg[imode][icent]->Draw("SAME");

      TLine *lineLFit_r = new TLine(IMcut_3Sigma[imode][icent][0], 0., IMcut_3Sigma[imode][icent][0], yscale*(maxCounts + maxCountsErr));
      lineLFit_r->SetLineColor(kRed);
      lineLFit_r->SetLineStyle(2);
      lineLFit_r->SetLineWidth(2);
      lineLFit_r->Draw("SAME");
      
      TLine *lineRFit_r = new TLine(IMcut_3Sigma[imode][icent][1], 0., IMcut_3Sigma[imode][icent][1], yscale*(maxCounts + maxCountsErr));
      lineRFit_r->SetLineColor(kRed);
      lineRFit_r->SetLineStyle(2);
      lineRFit_r->SetLineWidth(2);
      lineRFit_r->Draw("SAME");

      TBox *bbox;
      // if(imode == 0) bbox = new TBox(labeloff+0.035, 0.35*maxCounts, labeloff+0.047, 0.60*maxCounts);
      // else bbox = new TBox(labeloff+0.035, 0.35*maxCounts, labeloff+0.057, 0.60*maxCounts);
      bbox = new TBox(labeloff+0.035, 0.35*maxCounts, labeloff+0.057, 0.60*maxCounts);
      bbox->SetLineColor(2);
      bbox->SetLineWidth(2);
      bbox->SetFillStyle(0);
      bbox->Draw("same");

      TLatex *tex_yield = new TLatex(labeloff, 1.5*maxCounts, Form("S=%.2e, B=%.2e", nSignal, nBackground));
      tex_yield->SetTextFont(42);
      tex_yield->SetTextSize(0.035);
      tex_yield->Draw("same");

      TLatex *tex_purity = new TLatex(labeloff, 1.35*maxCounts, Form("S/(S+B) = %6.2f", nSignal/nTotal));
      tex_purity->SetTextFont(42);
      tex_purity->SetTextSize(0.035);
      tex_purity->Draw("same");

      TLatex *tex_significance = new TLatex(labeloff, 1.2*maxCounts, Form("S/#sqrt{S+B} = %6.2f", nSignal/sqrt(nTotal)));
      tex_significance->SetTextFont(42);
      tex_significance->SetTextSize(0.035);
      tex_significance->Draw("same");
      
      TLatex *tex_signif = new TLatex(labeloff, 1.05*maxCounts, Form("#mu = %7.4f #pm %7.4f", IMcut_mean[imode][icent], IMcut_meanErr[imode][icent]));
      tex_signif->SetTextFont(42);
      tex_signif->SetTextSize(0.035);
      tex_signif->Draw("same");
      
      TLatex *tex_sigma = new TLatex(labeloff, 0.9*maxCounts, Form("#sigma_{mass} = %7.4f #pm %7.4f", (IMcut_3Sigma[imode][icent][1]-IMcut_3Sigma[imode][icent][0])/(2*nsigma), IMcut_SigmaErr[imode][icent]));
      tex_sigma->SetTextFont(42);
      tex_sigma->SetTextSize(0.035);
      tex_sigma->Draw("same");

      TLatex *tex_chi2ndf = new TLatex(labeloff, 0.75*maxCounts, Form("#frac{#chi2}{NDF} = #frac{%.2f}{%.2f}", chi2, ndf));
      tex_chi2ndf->SetTextFont(42);
      tex_chi2ndf->SetTextSize(0.035);
      tex_chi2ndf->Draw("same");

      TLatex *tex_yrange = new TLatex(labeloff+0.036, 0.5*maxCounts, yName[imode]);
      tex_yrange->SetTextFont(42);
      tex_yrange->SetTextSize(0.035);
      tex_yrange->Draw("same");

      TLatex *tex_ptrange = new TLatex(labeloff+0.036, 0.4*maxCounts, ptName[imode]);
      tex_ptrange->SetTextFont(42);
      tex_ptrange->SetTextSize(0.035);
      tex_ptrange->Draw("same");

      TLatex *tex_fitfunc = new TLatex(labeloff, -0.27*maxCounts, "Fitting function: Student_t + Pol1 function");
      tex_fitfunc->SetTextFont(42);
      tex_fitfunc->SetTextSize(0.035);
      tex_fitfunc->Draw("same");

      TLegend *legT1 = new TLegend(0.62, 0.65, 0.9, 0.92);
      legT1->SetFillColor(10);
      legT1->SetFillStyle(0);
      legT1->SetLineStyle(3004);
      legT1->SetLineColor(10);
      legT1->SetLineWidth(0.);
      legT1->SetTextFont(42);
      legT1->SetTextSize(0.04);
      legT1->SetHeader(Form("#Lambda %s Au+Au", centName[icent].Data()));
      legT1->AddEntry(invMassVsCent_sig[imode][icent], "total - comb. BG", "pe");
      legT1->AddEntry(fInvMassVsCent_tot[imode][icent], "total - comb. fit", "l");
      legT1->AddEntry(fInvMassVsCentReal_bg[imode][icent], "residual bkg", "l");
      legT1->AddEntry(fInvMassVsCentReal_sig[imode][icent], "signal", "l");
      legT1->Draw();

      if(imode == 2 && icent == nCent-1)
      {
        // invmass_ncent[imode]->Close();

        TFile *fhis_cent0540 = new TFile("~/workspace/3p0/kfp/Ld_his_cent0540_acc3GeV_P24ia.root", "Recreate");
        fhis_cent0540->cd();

        invMassVsCent_tot[imode][icent]->Write("hist_tot");
        invMassVsCent_bg[imode][icent]->Write("hist_bg");
        invMassVsCent_sig[imode][icent]->Write("hist_sig");
        invMassVsCentReal_sig[imode][icent]->Write("hist_realsig");
        fInvMassVsCentReal_bg[imode][icent]->Write("residual_bkg");
        fInvMassVsCent_tot[imode][icent]->Write("tot_minus_bk_fit");

        fhis_cent0540->Close();

        invmass_ncent[imode]->cd();
      }
//      if(imode == 1 && icent == nCent-2)
//      {
//        // invmass_ncent[imode]->Close();
//
//        TFile *fhis_cent080 = new TFile("~/workspace/3p0/kfp/Ld_his_cent080_accMax_P24ia.root", "Recreate");
//        fhis_cent080->cd();
//
//        invMassVsCent_tot[imode][icent]->Write("hist_tot");
//        invMassVsCent_bg[imode][icent]->Write("hist_bg");
//        invMassVsCent_sig[imode][icent]->Write("hist_sig");
//        invMassVsCentReal_sig[imode][icent]->Write("hist_realsig");
//        fInvMassVsCentReal_bg[imode][icent]->Write("residual_bkg");
//        fInvMassVsCent_tot[imode][icent]->Write("tot_minus_bk_fit");
//
//        fhis_cent080->Close();
//
//        invmass_ncent[imode]->cd();
//      }

      ca_invMass->Write(Form("ca_invmass_cent%d", icent));
      if(icent == nCent-1) ca_invMass->Print(Form("ca_invmass_mode%d_054.pdf", imode));
      p0_invMass->Update();
      p1_invMass->Update();
    }
    invmass_ncent[imode]->Close();
  }
  
  //=========================================================================================
  
  auto fs = new TFile("fcent.root", "recreate");
  fs->cd();
  TCanvas *cCent = new TCanvas("ca_cent", "cCent", 800, 700);
  for(int imode=0; imode<nMode; imode++){
    cCent->cd();
    gPad->SetTopMargin(0.04);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.12);
    gPad->SetGrid(1);
    
    hCountsCent[imode]->SetMarkerStyle(24);
    hCountsCent[imode]->GetXaxis()->SetTitle("Centrality");
    hCountsCent[imode]->GetYaxis()->SetTitle("Counts");
    hCountsCent[imode]->SetMarkerColor(kBlue);
    auto legT = new TLegend(0.5, 0.65, 0.9, 0.88);
    legT->SetFillStyle(0);
    legT->SetFillColor(0);
    legT->SetBorderSize(0);
    legT->SetHeader(Form("#Lambda Au+Au @ 3.0 GeV"));

    for(int i=0; i<nCent;i++)hCountsCent[imode]->GetXaxis()->SetBinLabel(i+1, centName[i]);
    hCountsCent[imode]->Draw("pe");
    legT->Draw("same");
    
    TLatex *tex_yptrange = new TLatex(.4, .5, "");
    tex_yptrange->SetTextSize(0.035);
    tex_yptrange->DrawLatexNDC(.5, .65, yName[imode]);
    tex_yptrange->DrawLatexNDC(.5, .6, ptName[imode]);
    //=========================================================================================
    hCountsCent[imode]->Write();
    cCent->Write(Form("counts_mode%d", imode));
    cCent->Print(Form("ca_counts_mode%d.pdf", imode));
    cCent->Update();
  }
  fs->Close();
  //============================================================================
  return; 
}

