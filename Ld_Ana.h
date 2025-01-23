/*!
 * \class  Ld_Ana
 * \brief
 *
 * -------------------------------------------------------------------------
 *
 * -------------------------------------------------------------------------
 */
#ifndef Ld_Ana_hh
#define Ld_Ana_hh

// Include files
//#include "RtypesCore.h"
#include "Rtypes.h"
//#include "StTrackDef.h"
#include <string>
#include <iostream>
using std::string;

class TFile;
class TH1D;
class TH1F;
class TH2D;
class TH2F;
class TH3F;
class TProfile;
class TProfile2D;
class TVector2;
class TF1;
class TNtuple;
class EventClass;

typedef unsigned int UINT;

class Ld_Ana{
public:
    Ld_Ana(); // constructor
    Ld_Ana(const Char_t *outfile, EventClass* in_p_eve);
    ~Ld_Ana();

    void Clear();
    Int_t Init();
    Int_t NFInit();
    Int_t ProcessEvent();
    Int_t End();

private:
    EventClass *p_eve;

    // data member
    int mEventCounter;
    string mFileName;
    TFile *mFile;

    void DefHistograms();
    void WriteHistograms();
    void DeleteHistograms();
    Int_t Centrality(int gRefMult);
    Float_t Calibration_p(double nsigma, double mP);
    Float_t Calibration_pi(double nsigma, double mP);
    Float_t ShiftAngle(Float_t angle);
    Double_t isinpTBin(float pt, int mode);
    Double_t isinRapBin(float rap, int mode);

    // def histograms
    TH1D *hCentrality;
    TH1D *h_psi1shift_EPD[10];
    TH1D *hCentrality_ld;
    TH2D *hacc_ld;
    TH1D *hvz_ld;
    TH1D *hcountrefmult_ld;
    TH1D *htofmult_ld;
    TH2D *hTofMvsRM_ld;

    TH2F *hdedx_p_p;
    TH2F *hdedx_p_pi;
    TH2F *hmass2_p_p;
    TH2F *hmass2_p_pi;
    
    TH2F *heta_pt_p;
    TH2F *heta_pt_pi;
    
    TH2F *hrap_pt_p;
    TH2F *hrap_pt_pi;

    TH3F *hMassLdRapidityvsPt[9];
    TH1F *hMassLd[3][9];
    
    TH3F *FlowInvMEp[2];
    TH3F *FlowInvMEp_pt;
    
//    TH3F *FlowInvMEp_pt_diffy[16];
//    
//    TH3F *FlowInvMEp_pt_diffy_cent[5][10];
    
//    TH3F *FlowInvMEp_new[1];
//    TH3F *FlowInvMEp_m1n2_pt_new;
    
//    TH3F *InvMassv2_y_cent[9];
    TH3F *InvMassv2_pt_cent[5][10];

    // TNtuple *fqaCheck;
};
#endif
