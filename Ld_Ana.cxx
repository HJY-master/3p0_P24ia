// My header files
#include "Ld_Ana.h"
#include "EventClass.h"

// ROOT header files
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"

#include "badrun.h"

// C/C++ heder
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
float pi = TMath::Pi();
const float ycm = -1.045;

// EPD eABeDtB 1st res
Float_t ep_res[9] = {0.247145, 0.424986, 0.606138, 0.717986, 0.737337, 0.692473, 0.585713, 0.432261, 0.296892};
// 2nd res
Float_t ep_res21[9] = {0.0394116, 0.119806, 0.256308, 0.377406, 0.402134, 0.346645, 0.237622, 0.12414, 0.0572088};

float piMass = 0.13957;
float kMass = 0.493677;
float pMass = 0.938272;
float dMass = 1.875613;
float tMass = 2.80892;
float he3Mass = 2.80932;
float he4Mass = 3.727417;
float htriton_mass = 2.99131;

Ld_Ana::Ld_Ana(const Char_t *outfile, EventClass* in_p_eve) : mEventCounter(0), mFile(0)
{
    mFileName = outfile;
    gRandom->SetSeed(0);
    p_eve = in_p_eve;
}

Ld_Ana::~Ld_Ana() { /* noop */ }

Int_t Ld_Ana::Init()
{
    DefHistograms();

    // fqaCheck = new TNtuple("qaCheck", "qa for sys", "opt:centnumber:invM:p_ld:pt:rap:l:dl:ldl:chi2ndf:chi2topo:chi2primary_proton:chi2primary_pi:nhits_proton:nhits_pi:mPhiPsi:fResInv:gweight");

    return 0;
}

Int_t Ld_Ana::NFInit()
{
    // DefHistograms();

    // fqaCheck = new TNtuple("qaCheck", "qa for sys", "opt:centnumber:invM:p_ld:pt:rap:l:dl:ldl:chi2ndf:chi2topo:chi2primary_proton:chi2primary_pi:nhits_proton:nhits_pi:mPhiPsi:fResInv:gweight");

    return 0;
}

void Ld_Ana::Clear()
{
    p_eve->InitValue();
}

Int_t Ld_Ana::End()
{
    cout << "Ld_Ana::End()\n";
    cout << "\tProcessed " << mEventCounter << " events." << endl;

    // Output histograms
    mFile = new TFile(mFileName.c_str(), "RECREATE", "Output of Ld_Ana: 3-PC analysis", 9);
    cout << "\tHistograms will be stored in file '" << mFileName.c_str() << "'" << endl;

    // fqaCheck->Write();
    WriteHistograms();

    // Write histos to file and close it.
    if( mFile )
    {
        mFile->Write();
        mFile->Close();
    }

    // DeleteHistograms();
    
    return 0;
}

const int nrap = 18;
const float kRapMin[nrap] = {-1.8, -1.0, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
const float kRapMax[nrap] = {0.8, 0.6, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

const int npt = 4;
const float kPtMin[npt] = {0.4, 0.2, 0.2, 0.2};
const float kPtMax[npt] = {0.8, 3.0, 2.0, 1.0};


// This method is called every event.
Int_t Ld_Ana::ProcessEvent()
{
    mEventCounter++;
    if( (mEventCounter%100000) == 0 )
    {
        cout << "Processed events = " << mEventCounter << endl;
    }

    if( !p_eve )
    {
        cout << "Ld_Ana::ProcessEvent() : p_eve is NULL !!" << endl;
        return -1;
    }

    // ---------------------------------------------------

    // flow analysis
    int runid = p_eve->get_runId();
    for(int ii=0; ii<nbadrun; ii++)
    {
        if(runid == badrun[ii]) return 0;
    }
    float gweight = p_eve->get_gweight();
    gweight = 1;

    // define centrality bins
    int centnumber = p_eve->get_fCentrality();
    if(centnumber<0 || centnumber>8) return 0;
    hCentrality->Fill(centnumber);

    int centIndex = -1;
    if(centnumber==8 || centnumber==7)  centIndex = 0;  //0-10%
    if(centnumber==6)  centIndex = 1; //10-20%
    if(centnumber==5)  centIndex = 2; //20-30%
    if(centnumber==4)  centIndex = 3; //30-40%
    if(centnumber==3)  centIndex = 4; //40-50%
    if(centnumber==2)  centIndex = 5; //50-60%
    if(centnumber==1)  centIndex = 6; //60-70%
    if(centnumber==0)  centIndex = 7; //70-80%

    float psi1shift_EPD = p_eve->get_psi_1_EPD_4(); // EPD eABeCtB 1st
    h_psi1shift_EPD[centnumber]->Fill(psi1shift_EPD);

    float px = p_eve->get_px();
    float py = p_eve->get_py();
    float pz = p_eve->get_pz();
    float invM = p_eve->get_mass();
    float fResInv = 1.0 / ep_res[8-centnumber];
    float fRes21Inv = 1.0 / ep_res21[8-centnumber];
    float px_pi = p_eve->get_px_pi();
    float py_pi = p_eve->get_py_pi();
    float pz_pi = p_eve->get_pz_pi();
    float px_proton = p_eve->get_px_pro();
    float py_proton = p_eve->get_py_pro();
    float pz_proton = p_eve->get_pz_pro();
    float dedx_proton = p_eve->get_dedx_pro();
    float dedx_pi = p_eve->get_dedx_pi();

    float chi2primary_proton = p_eve->get_chi2primary_pro();
    float chi2primary_pi = p_eve->get_chi2primary_pi();

    //======================================================================
    float nhits_proton = p_eve->get_nhits_pro();
    float nhits_pi = p_eve->get_nhits_pi();
    if(nhits_proton<15 || nhits_pi<15) return -1;

    TLorentzVector ptv(px, py, pz, sqrt(px*px + py*py + pz*pz + invM*invM));
    float rap = ptv.Rapidity() - ycm;
    rap = -1.0*rap;

    float p_ld = sqrt(px*px + py*py + pz*pz);
    float p_pi = sqrt(px_pi*px_pi + py_pi*py_pi + pz_pi*pz_pi);
    float p_proton = sqrt(px_proton*px_proton + py_proton*py_proton + pz_proton*pz_proton);
    float pt = sqrt(px*px + py*py);
    float pt_pi = sqrt(px_pi*px_pi + py_pi*py_pi);
    float pt_proton = sqrt(px_proton*px_proton + py_proton*py_proton);
    float phi = atan2(py, px);
    float mPhiPsi = phi - psi1shift_EPD;
    
    TLorentzVector ptv_p(px_proton, py_proton, pz_proton, sqrt(px_proton*px_proton + py_proton*py_proton + pz_proton*pz_proton + pMass*pMass));
    float rap_p = ptv_p.Rapidity() - ycm;
    rap_p = -1.0*rap_p;
    
    TLorentzVector ptv_pi(px_pi, py_pi, pz_pi, sqrt(px_pi*px_pi + py_pi*py_pi + pz_pi*pz_pi + piMass*piMass));
    float rap_pi = ptv_pi.Rapidity() - ycm;
    rap_pi = -1.0*rap_pi;
    
    TVector3 mom_p(px_proton, py_proton, pz_proton);
    float eta_p = mom_p.PseudoRapidity();
    
    TVector3 mom_pi(px_pi, py_pi, pz_pi);
    float eta_pi = mom_pi.PseudoRapidity();

//    int eff_bin = hphasespace_mceff->FindBin(rap, pt);
//    double eff = hphasespace_mceff->GetBinContent(eff_bin);

    if(pt_proton < 0.1) return -1;
    if(pt_pi < 0.1) return -1;
    
    if(eta_p > 0 || eta_pi > 0) return -1;

    float l = p_eve->get_l();
    float ldl = p_eve->get_ldl();
    float chi2ndf = p_eve->get_chi2ndf();
    float chi2topo = p_eve->get_chi2topo();

    hdedx_p_p->Fill(p_proton, dedx_proton);
    hdedx_p_pi->Fill(p_pi, dedx_pi);
    
    heta_pt_p->Fill(eta_p, pt_proton);
    heta_pt_pi->Fill(eta_pi, pt_pi);
    
    hrap_pt_p->Fill(rap_p, pt_proton);
    hrap_pt_pi->Fill(rap_pi, pt_pi);

    if(l<=3 || ldl<=5 || chi2topo>=4 || chi2ndf>=4 || chi2primary_proton<=7 || chi2primary_pi<=10) return -1;
    
    hMassLdRapidityvsPt[8-(int)centnumber]->Fill(invM, rap, pt, gweight);
    hMassLd[0][8-(int)centnumber]->Fill(invM, gweight);
    if(isinpTBin(pt, 0) && isinRapBin(rap, 0)) { hMassLd[1][8-(int)centnumber]->Fill(invM, gweight); } // y:[-1.4, 0.2] pt:[0.4, 2.4]
    if(isinpTBin(pt, 1) && isinRapBin(rap, 1)) { hMassLd[2][8-(int)centnumber]->Fill(invM, gweight); } // y:[-1.0, 0.0] pt:[0.4, 0.8]

    hCentrality_ld->Fill(8-centnumber);
    hacc_ld->Fill(rap, pt);

    // to get Ld yield in small centraltiy bins, for later EP resolution in large cent.range
    if(centIndex < 0) return -1;

    mPhiPsi = ShiftAngle(mPhiPsi);
    
//    if(isinpTBin(pt, 2) && isinRapBin(rap, 0))
//    {
//        InvMassv2_y_cent[8-(int)centnumber]->Fill(rap, cos(2*mPhiPsi), invM, gweight);
//    }
//    
//    if(isinpTBin(pt, 1) && isinRapBin(rap, 1))
//    {
//        InvMassv2_pt_cent[8-(int)centnumber]->Fill(pt, cos(2*mPhiPsi), invM, gweight);
//    }
    
    for(int ii = 0; ii < 10; ii++)
    {
        if(isinpTBin(pt, 1) && isinRapBin(rap, 2+ii)) InvMassv2_pt_cent[8-(int)centnumber][ii]->Fill(pt, invM, cos(2*mPhiPsi)*fRes21Inv, gweight);
    }

    if(centnumber>=4 && centnumber<=7)//5-40%
    {
        // v2 res21
        
        if(isinpTBin(pt, 2) && isinRapBin(rap, 0)) FlowInvMEp[0]->Fill(rap, mPhiPsi, invM, fRes21Inv*gweight);
        if(isinpTBin(pt, 3) && isinRapBin(rap, 0)) FlowInvMEp[1]->Fill(rap, mPhiPsi, invM, fRes21Inv*gweight);
        
        if(isinpTBin(pt, 1) && isinRapBin(rap, 1)) FlowInvMEp_pt->Fill(pt, mPhiPsi, invM, fRes21Inv*gweight);
        
//        for(int ii = 0; ii < 16; ii++)
//        {
//            if(isinpTBin(pt, 1) && isinRapBin(rap, 2+ii)) FlowInvMEp_pt_diffy[ii]->Fill(pt, mPhiPsi, invM, fRes21Inv*gweight);
//        }
    }
    
    // centnumber = 8 (0-5%), 7 (5-10%), 6 (10-20%), 5 (20-30%), 4 (30-40%)
    
//    if(centnumber >= 4 && centnumber <= 8)
//    {
//        for(int ii = 0; ii < 10; ii++)
//        {
//            if(isinpTBin(pt, 1) && isinRapBin(rap, 2+ii)) FlowInvMEp_pt_diffy_cent[8-(int)centnumber][ii]->Fill(pt, mPhiPsi, invM, fRes21Inv*gweight);
//        }
//    }

    return 0;
}

void Ld_Ana::DefHistograms()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();

    hCentrality = new TH1D("hCentrality", "Centrality", 9, -0.5, 8.5);
    for(int i = 0; i < 10; i++)
    {
        h_psi1shift_EPD[i] = new TH1D(Form("EPDshift_cent%d", i), "", 314, -pi, pi);
    }
    hCentrality_ld = new TH1D("hCentrality_ld", "Centrality Lambda", 9, -0.5, 8.5);
    hacc_ld = new TH2D("hacc_ld", "htriton acceptance", 300, -1.5, 1.5, 500, 0, 10);
    hvz_ld = new TH1D("hvz_ld", "Vz;Vz(cm);Counts", 300, 0, 300);
    hcountrefmult_ld = new TH1D("hcountrefmult_ld", "hcountrefmult;countrefmult;N_{evt}", 600, 0, 600);
    htofmult_ld = new TH1D("htofmult_ld", "htofmult;tofmult;N_{evt}", 1000, 0, 1000);
    hTofMvsRM_ld = new TH2D("hTofMvsRM", "hTofMvsRM", 600, 0, 600, 1000, 0, 1000);

    hdedx_p_p = new TH2F("hdedx_p_p", "", 2000, -6., 14., 1500, 0., 150.);
    hdedx_p_pi = new TH2F("hdedx_p_pi", "", 1200, -6., 6., 1500, 0., 100.);
    hmass2_p_p = new TH2F("hmass2_p_p", "", 1200, -6., 6., 1600, -0.2, 15.8);
    hmass2_p_pi = new TH2F("hmass2_p_pi", "", 1200, -6., 6., 1600, -0.2, 15.8);
    
    heta_pt_p = new TH2F("heta_pt_p", "", 400, -3, 1, 500, 0, 5);
    hrap_pt_p = new TH2F("hrap_pt_p", "", 400, -3, 1, 500, 0, 5);
    heta_pt_pi = new TH2F("heta_pt_pi", "", 400, -3, 1, 500, 0, 5);
    hrap_pt_pi = new TH2F("hrap_pt_pi", "", 400, -3, 1, 500, 0, 5);

    for(int icent = 0; icent < 9; icent++){
        hMassLdRapidityvsPt[icent] = new TH3F(Form("hMassLdRapidityvsPt_cent%d", icent),Form("hMassLdRapidityvsPt_cent%d", icent),160,1.090,1.150,200,-2.,2.,200,0,10);
    }
    for(int imode = 0; imode < 3; imode++){
        for(int icent = 0; icent < 9; icent++){
            hMassLd[imode][icent] = new TH1F(Form("hMassLd_mode%d_cent%d", imode, icent),Form("hMassLd_mode%d_cent%d", imode, icent),150,1.050,1.200);
        }
    }
    
    for(int i = 0; i < 2; i++)
    {
        FlowInvMEp[i] = new TH3F(Form("FlowInvMEp_yptmode%d", i),"", 200,-2.0,2.0,6,0,TMath::Pi(), 100, 1.090, 1.150);
    }
    FlowInvMEp_pt = new TH3F("FlowInvMEp_pt","", 200,0.0,5.0,6,0,TMath::Pi(), 100, 1.090, 1.150);
    
    for(int icent = 0; icent < 9; icent++){
        for(int i = 0; i < 10; i++)
        {
            InvMassv2_pt_cent[icent][i] = new TProfile2D(Form("InvMassv2_yptmode%d_pt_diffy_cent%d", i, icent),"", 30, 0.0, 5.0, 120, 1.080, 1.150);
        }
    }
     
//    for(int i = 0; i < 16; i++)
//    {
//        FlowInvMEp_pt_diffy[i] = new TH3F(Form("FlowInvMEp_yptmode%d_pt_diffy", i),"", 200,0.0,5.0,6,0,TMath::Pi(), 100, 1.090, 1.150);
//    }
//    
//    for(int icent = 0; icent < 5; icent++){
//        for(int i = 0; i < 10; i++)
//        {
//            FlowInvMEp_pt_diffy_cent[icent][i] = new TH3F(Form("FlowInvMEp_yptmode%d_pt_diffy_cent%d", i, icent),"", 200,0.0,5.0,6,0,TMath::Pi(), 100, 1.090, 1.150);
//        }
//    }
    
}

void Ld_Ana::WriteHistograms()
{
    hCentrality->Write();
    for(int i = 0; i < 10; i++)
    {
        h_psi1shift_EPD[i]->Write();
    }
    hCentrality_ld->Write();
    hacc_ld->Write();
    hvz_ld->Write();
    hcountrefmult_ld->Write();
    htofmult_ld->Write();
    hTofMvsRM_ld->Write();

    hdedx_p_p->Write();
    hdedx_p_pi->Write();
    hmass2_p_p->Write();
    hmass2_p_pi->Write();
    
    heta_pt_p->Write();
    hrap_pt_p->Write();
    heta_pt_pi->Write();
    hrap_pt_pi->Write();

    for(int icent = 0; icent < 9; icent++){
        hMassLdRapidityvsPt[icent]->Write();
    }
    for(int imode = 0; imode < 3; imode++){
        for(int icent = 0; icent < 9; icent++){
            hMassLd[imode][icent]->Write();
        }
    }
    
    for(int i = 0; i < 2; i++)
    {
        FlowInvMEp[i]->Write();
    }
    FlowInvMEp_pt->Write();
    
    for(int icent = 0; icent < 9; icent++){
        for(int i = 0; i < 10; i++)
        {
            InvMassv2_pt_cent[icent][i]->Write();
        }
    }
    
//    for(int i = 0; i < 16; i++)
//    {
//        FlowInvMEp_pt_diffy[i]->Write();
//    }
//    
//    for(int icent = 0; icent < 5; icent++){
//        for(int i = 0; i < 10; i++)
//        {
//            FlowInvMEp_pt_diffy_cent[icent][i]->Write();
//        }
//    }
            
}

Double_t Ld_Ana::isinpTBin(float pt, int mode)
{
    return (pt <= kPtMax[mode] && pt >= kPtMin[mode]);
}

Double_t Ld_Ana::isinRapBin(float rap, int mode)
{
    return (rap <= kRapMax[mode] && rap >= kRapMin[mode]);
}

Float_t Ld_Ana::ShiftAngle(Float_t angle)
{
  double pi = TMath::Pi();
  // Float_t shift_angle = angle + pi/2.;
  Float_t shift_angle = angle;
  // if(shift_angle < 0.0) shift_angle += 2*pi;
  // if(shift_angle > pi) shift_angle -= pi;
  
  // TODO: !!!
  // if(shift_angle > pi.) shift_angle = 2*pi - shift_angle;

  for(int i = 0; i < 3; i++)
  {
    if(shift_angle < -1.0*pi) shift_angle += 2.0*pi;
    if(shift_angle > 1.0*pi) shift_angle -= 2.0*pi;
  }
  shift_angle = fabs(shift_angle);
  return shift_angle;
}
