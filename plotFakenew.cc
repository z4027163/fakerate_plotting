
#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TH1.h"

#include "RooRealVar.h"
#include "TLatex.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgusBG.h"
#include "RooBernstein.h"
#include "RooProduct.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TChain.h"
#include "TFile.h"
#include "TArrow.h"
#include <map>
#include "TStyle.h"
#include "TLegend.h"

#include <memory>
#include <dirent.h>
#include <sys/types.h>

using namespace std;
using namespace RooFit;
class plotFakenew
{
public:
  void initPol1(double &, double &, TH1 *);
  TF1 *pol1gauss2c(TH1 *, double, double);
  double dRatio(double a, double ae, double b, double be)
  {
    return TMath::Sqrt(((ae * ae) / (b * b)) + ((a * a * be * be) / (b * b * b * b)));
  }
  void sbsDistributions(TChain *, string, string what = "nothing", string path = "ks/", TString binning = "pT");
  void list_dir(TChain *, const char *);
  void list_dir_file(TChain *, string ob = "bla");
  void loopOverChain(TChain *, string, string, string, TString);

  plotFakenew();
  ~plotFakenew();

private:
  bool fVerbose = true;
};

void plotFakenew::initPol1(double &p0, double &p1, TH1 *h)
{
  int EDG(4), NB(EDG + 1);
  int lbin(1), hbin(h->GetNbinsX() + 1);
  double fLo = 99.0, fHi = -99.0;
  if (fLo < fHi)
  {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }
  double xlo = h->GetBinLowEdge(lbin);

  double dx = h->GetBinLowEdge(hbin) - xlo;
  double ylo = h->Integral(lbin, lbin + EDG) / NB;
  double yhi = h->Integral(hbin - EDG - 1, hbin - 1) / NB;

  p1 = (yhi - ylo) / dx;
  p0 = ylo - p1 * xlo;
  if (1 || fVerbose)
  {
    cout << "lbin: " << lbin << " hbin: " << hbin
         << " ylo: " << ylo << " yhi: " << yhi << " dx: " << dx
         << " p0: " << p0 << " p1: " << p1 << endl;
  }
}
double iF_pol1(double *x, double *par)
{
  return par[0] + par[1] * x[0];
}
plotFakenew::~plotFakenew(void)
{
  // destructor
}
plotFakenew::plotFakenew(void)
{
  // constructor
}

void plotFakenew::sbsDistributions(TChain *tC, string sample = "bla", string what, string path, TString binning)
{
  cout << "plotFake::sbsDistributions(" << sample << ", " << what << ")" << endl;
  string sbsControlPlotsFileName = what.c_str();

  int NBINS = 50;
  double fMassLo(0.), fMassHi(0.);
  if (string::npos != sample.find("ks"))
  {
    fMassLo = 0.450;
    fMassHi = 0.550;
    NBINS = 50;
  }
  else if (string::npos != sample.find("d0"))
  {
    fMassLo = 1.8;
    fMassHi = 1.90;
  }
  else if (string::npos != sample.find("lambda"))
  {
    fMassLo = 1.095;
    NBINS = 45;
    fMassHi = 1.140;
  }
  else if (string::npos != sample.find("phi"))
  {
    fMassLo = 1.0;
    fMassHi = 1.04;
    NBINS = 40;
  }
  string name = "Ks_sbs_";

  int NREG_pt = 8;

  TString massbin_pt[NREG_pt];
  massbin_pt[0] = "bin0";
  massbin_pt[1] = "bin1";
  massbin_pt[2] = "bin2";
  massbin_pt[3] = "bin3";
  massbin_pt[4] = "bin4";
  massbin_pt[5] = "bin5";
  massbin_pt[6] = "bin6";
  massbin_pt[7] = "allbin";

  TH1D *hv0_kin_pt_bin[NREG_pt];
  TH1D *hv0_Mass_bin[NREG_pt];

  TH1D *hv0_kin_pt_bin_muid[NREG_pt];
  TH1D *hv0_Mass_bin_muid[NREG_pt];

  TH1D *hv0_kin_pt_bin_softid[NREG_pt];
  TH1D *hv0_Mass_bin_softid[NREG_pt];
  TH1D *hv0_kin_pt_bin_medid[NREG_pt];
  TH1D *hv0_Mass_bin_medid[NREG_pt];

  TH1D *hv0_kin_pt_bin_muid_30[NREG_pt];
  TH1D *hv0_Mass_bin_muid_30[NREG_pt];

  TH1D *hv0_kin_pt_bin_muid_40[NREG_pt];
  TH1D *hv0_Mass_bin_muid_40[NREG_pt];
  TH1D *hv0_kin_pt_bin_muid_45[NREG_pt];
  TH1D *hv0_Mass_bin_muid_45[NREG_pt];
  TH1D *hv0_kin_pt_bin_muid_20[NREG_pt];
  TH1D *hv0_Mass_bin_muid_20[NREG_pt];
  TH1D *hv0_kin_pt_bin_muid_50[NREG_pt];
  TH1D *hv0_Mass_bin_muid_50[NREG_pt];
  TH1D *hv0_kin_pt_bin_muid_55[NREG_pt];
  TH1D *hv0_Mass_bin_muid_55[NREG_pt];
  TH1D *hv0_kin_pt_bin_muid_60[NREG_pt];
  TH1D *hv0_Mass_bin_muid_60[NREG_pt];

  for (int k = 0; k < NREG_pt; k++)
  {
    // double var[]={0.,4.0,8.0,12.0,16.0,20.0, 30.0,50.};
    int ptbins = 150;
    hv0_kin_pt_bin[k] = new TH1D(Form("hv0_kin_pt_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);
    hv0_Mass_bin[k] = new TH1D(Form("hv0_Mass_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);

    hv0_Mass_bin_muid[k] = new TH1D(Form("hv0_Mass_muid_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid[k] = new TH1D(Form("hv0_kin_pt_muid_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muid;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);
    hv0_Mass_bin_softid[k] = new TH1D(Form("hv0_Mass_softid_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_softid", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_softid[k] = new TH1D(Form("hv0_kin_pt_softid_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_softid;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);
    hv0_Mass_bin_medid[k] = new TH1D(Form("hv0_Mass_medid_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_medid", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_medid[k] = new TH1D(Form("hv0_kin_pt_medid_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_medid;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);

    hv0_Mass_bin_muid_30[k] = new TH1D(Form("hv0_Mass_muid_30_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid_mumva30", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid_30[k] = new TH1D(Form("hv0_kin_pt_muid_30_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muidmumva30;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);

    hv0_Mass_bin_muid_40[k] = new TH1D(Form("hv0_Mass_muid_40_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid_mumva40", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid_40[k] = new TH1D(Form("hv0_kin_pt_muid_40_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muidmumva40;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);

    hv0_Mass_bin_muid_45[k] = new TH1D(Form("hv0_Mass_muid_45_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid_mumva45", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid_45[k] = new TH1D(Form("hv0_kin_pt_muid_45_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muidmumva45;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);

    hv0_Mass_bin_muid_20[k] = new TH1D(Form("hv0_Mass_muid_20_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid_mumva20", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid_20[k] = new TH1D(Form("hv0_kin_pt_muid_20_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muidmumva20;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);

    hv0_Mass_bin_muid_50[k] = new TH1D(Form("hv0_Mass_muid_50_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid_mumva50", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid_50[k] = new TH1D(Form("hv0_kin_pt_muid_50_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muidmumva50;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);
    hv0_Mass_bin_muid_55[k] = new TH1D(Form("hv0_Mass_muid_55_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid_mumva55", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid_55[k] = new TH1D(Form("hv0_kin_pt_muid_55_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muidmumva55;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);

    hv0_Mass_bin_muid_60[k] = new TH1D(Form("hv0_Mass_muid_60_%s", massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid_mumva60", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo, fMassHi);
    hv0_kin_pt_bin_muid_60[k] = new TH1D(Form("hv0_kin_pt_muid_60_%s", massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muidmumva60;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins, 0, 50.);
    hv0_Mass_bin[k]->SetMinimum(0.);
  }

  UInt_t run;
  UInt_t luminosityBlock;
  UInt_t nMuon;
  Bool_t HLT_Ele30 = true;
  Float_t Muon_pt[150];
  Bool_t Muon_softMvaId[150];
  Float_t Muon_softMva[150];
  Bool_t Muon_softId[150];
  Bool_t Muon_mediumId[150];
  Bool_t Muon_isGlobal[150];
  Bool_t Muon_isTracker[150];
  Bool_t Muon_looseId[150];
  Float_t MuonId_hlt_pt[150];

  UInt_t nV0;
  Float_t V0_doca[150];
  Float_t V0_kin_cosAlphaXY[150];
  Float_t V0_kin_eta[150];
  Float_t V0_kin_lxy[150];
  Float_t V0_kin_mass[150];
  Float_t V0_kin_massErr[150];
  Float_t V0_kin_phi[150];
  Float_t V0_kin_pt[150];
  Float_t V0_kin_slxy[150];
  Float_t V0_kin_vtx_chi2dof[150];
  Float_t V0_kin_vtx_prob[150];
  Float_t V0_mass[150];
  Float_t V0_trk1_eta[150];
  Float_t V0_trk1_phi[150];
  Float_t V0_trk1_pt[150];
  Float_t V0_trk2_eta[150];
  Float_t V0_trk2_phi[150];
  Float_t V0_trk2_pt[150];
  Int_t V0_kin_valid[150];
  Int_t V0_trk1_mu_index[150];
  Int_t V0_trk2_mu_index[150];

  Int_t V0_gen_pdgId[150];
  Int_t V0_gen_trk1_mpdgId[150];
  Int_t V0_gen_trk1_pdgId[150];
  Int_t V0_gen_trk2_mpdgId[150];
  Int_t V0_gen_trk2_pdgId[150];
  Float_t V0_kin_sipPV[150];
  Float_t V0_trk1_sip[150];
  Float_t V0_trk2_sip[150];

  TBranch *b_run;
  TBranch *b_luminosityBlock;
  TBranch *b_nMuon;
  TBranch *b_HLT_Ele30;
  TBranch *b_Muon_pt;
  TBranch *b_Muon_softMvaId;
  TBranch *b_Muon_softMva;
  TBranch *b_Muon_softId;
  TBranch *b_Muon_mediumId;
  TBranch *b_Muon_isGlobal;
  TBranch *b_Muon_isTracker;
  TBranch *b_Muon_looseId;
  TBranch *b_MuonId_hlt_pt;

  TBranch *b_nV0;
  TBranch *b_V0_doca;
  TBranch *b_V0_kin_cosAlphaXY;
  TBranch *b_V0_kin_eta;
  TBranch *b_V0_kin_lxy;
  TBranch *b_V0_kin_mass;
  TBranch *b_V0_kin_massErr;
  TBranch *b_V0_kin_phi;
  TBranch *b_V0_kin_pt;
  TBranch *b_V0_kin_slxy;
  TBranch *b_V0_kin_vtx_chi2dof;
  TBranch *b_V0_kin_vtx_prob;
  TBranch *b_V0_mass;
  TBranch *b_V0_trk1_eta;
  TBranch *b_V0_trk1_phi;
  TBranch *b_V0_trk1_pt;
  TBranch *b_V0_trk2_eta;
  TBranch *b_V0_trk2_phi;
  TBranch *b_V0_trk2_pt;
  TBranch *b_V0_kin_valid;
  TBranch *b_V0_trk1_mu_index;
  TBranch *b_V0_trk2_mu_index;
  TBranch *b_V0_kin_sipPV;
  TBranch *b_V0_trk1_sip;
  TBranch *b_V0_trk2_sip;

  TBranch *b_V0_gen_pdgId;
  TBranch *b_V0_gen_trk1_mpdgId;
  TBranch *b_V0_gen_trk1_pdgId;
  TBranch *b_V0_gen_trk2_mpdgId;
  TBranch *b_V0_gen_trk2_pdgId;

  tC->SetBranchAddress("run", &run, &b_run);
  tC->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
  tC->SetBranchAddress("nMuon", &nMuon, &b_nMuon);

  if (string::npos != sample.find("trigger"))
  {
    tC->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30, &b_HLT_Ele30);
  }

  tC->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
  tC->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
  tC->SetBranchAddress("Muon_softMva", Muon_softMva, &b_Muon_softMva);
  tC->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
  tC->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
  tC->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
  tC->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
  tC->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);

  if (string::npos != sample.find("DoubleMuMu"))
  {
    tC->SetBranchAddress("MuonId_hlt_pt", MuonId_hlt_pt, &b_MuonId_hlt_pt);
  }

  if (string::npos != sample.find("ks"))
  {
    tC->SetBranchAddress("nks", &nV0, &b_nV0);
    tC->SetBranchAddress("ks_doca", V0_doca, &b_V0_doca);
    tC->SetBranchAddress("ks_kin_cosAlphaXY", V0_kin_cosAlphaXY, &b_V0_kin_cosAlphaXY);
    tC->SetBranchAddress("ks_kin_eta", V0_kin_eta, &b_V0_kin_eta);
    tC->SetBranchAddress("ks_kin_lxy", V0_kin_lxy, &b_V0_kin_lxy);
    tC->SetBranchAddress("ks_kin_mass", V0_kin_mass, &b_V0_kin_mass);
    tC->SetBranchAddress("ks_kin_massErr", V0_kin_massErr, &b_V0_kin_massErr);
    tC->SetBranchAddress("ks_kin_phi", V0_kin_phi, &b_V0_kin_phi);
    tC->SetBranchAddress("ks_kin_pt", V0_kin_pt, &b_V0_kin_pt);
    tC->SetBranchAddress("ks_kin_slxy", V0_kin_slxy, &b_V0_kin_slxy);
    tC->SetBranchAddress("ks_kin_vtx_chi2dof", V0_kin_vtx_chi2dof, &b_V0_kin_vtx_chi2dof);
    tC->SetBranchAddress("ks_kin_vtx_prob", V0_kin_vtx_prob, &b_V0_kin_vtx_prob);
    tC->SetBranchAddress("ks_mass", V0_mass, &b_V0_mass);
    tC->SetBranchAddress("ks_trk1_eta", V0_trk1_eta, &b_V0_trk1_eta);
    tC->SetBranchAddress("ks_trk1_phi", V0_trk1_phi, &b_V0_trk1_phi);
    tC->SetBranchAddress("ks_trk1_pt", V0_trk1_pt, &b_V0_trk1_pt);
    tC->SetBranchAddress("ks_trk2_eta", V0_trk2_eta, &b_V0_trk2_eta);
    tC->SetBranchAddress("ks_trk2_phi", V0_trk2_phi, &b_V0_trk2_phi);
    tC->SetBranchAddress("ks_trk2_pt", V0_trk2_pt, &b_V0_trk2_pt);
    tC->SetBranchAddress("ks_kin_valid", V0_kin_valid, &b_V0_kin_valid);
    tC->SetBranchAddress("ks_trk1_mu_index", V0_trk1_mu_index, &b_V0_trk1_mu_index);
    tC->SetBranchAddress("ks_trk2_mu_index", V0_trk2_mu_index, &b_V0_trk2_mu_index);

    tC->SetBranchAddress("ks_kin_sipPV", V0_kin_sipPV, &b_V0_kin_sipPV);
    tC->SetBranchAddress("ks_trk1_sip", V0_trk1_sip, &b_V0_trk1_sip);
    tC->SetBranchAddress("ks_trk2_sip", V0_trk2_sip, &b_V0_trk2_sip);

    /*if(string::npos != what.find("MC")){
       tC->SetBranchAddress("ks_gen_pdgId", V0_gen_pdgId, &b_V0_gen_pdgId);
       tC->SetBranchAddress("ks_gen_trk1_mpdgId",V0_gen_trk1_mpdgId, &b_V0_gen_trk1_mpdgId);
       tC->SetBranchAddress("ks_gen_trk2_mpdgId",V0_gen_trk2_mpdgId, &b_V0_gen_trk2_mpdgId);
       tC->SetBranchAddress("ks_gen_trk1_pdgId",V0_gen_trk1_pdgId, &b_V0_gen_trk1_pdgId);
       tC->SetBranchAddress("ks_gen_trk2_pdgId",V0_gen_trk2_pdgId, &b_V0_gen_trk2_pdgId);
    }*/
  }
  else if (string::npos != sample.find("lambda"))
  {
    tC->SetBranchAddress("nlambda", &nV0, &b_nV0);
    tC->SetBranchAddress("lambda_doca", V0_doca, &b_V0_doca);
    tC->SetBranchAddress("lambda_kin_cosAlphaXY", V0_kin_cosAlphaXY, &b_V0_kin_cosAlphaXY);
    tC->SetBranchAddress("lambda_kin_eta", V0_kin_eta, &b_V0_kin_eta);
    tC->SetBranchAddress("lambda_kin_lxy", V0_kin_lxy, &b_V0_kin_lxy);
    tC->SetBranchAddress("lambda_kin_mass", V0_kin_mass, &b_V0_kin_mass);
    tC->SetBranchAddress("lambda_kin_massErr", V0_kin_massErr, &b_V0_kin_massErr);
    tC->SetBranchAddress("lambda_kin_phi", V0_kin_phi, &b_V0_kin_phi);
    tC->SetBranchAddress("lambda_kin_pt", V0_kin_pt, &b_V0_kin_pt);
    tC->SetBranchAddress("lambda_kin_slxy", V0_kin_slxy, &b_V0_kin_slxy);
    tC->SetBranchAddress("lambda_kin_vtx_chi2dof", V0_kin_vtx_chi2dof, &b_V0_kin_vtx_chi2dof);
    tC->SetBranchAddress("lambda_kin_vtx_prob", V0_kin_vtx_prob, &b_V0_kin_vtx_prob);
    tC->SetBranchAddress("lambda_mass", V0_mass, &b_V0_mass);
    tC->SetBranchAddress("lambda_pion_eta", V0_trk1_eta, &b_V0_trk1_eta);
    tC->SetBranchAddress("lambda_pion_phi", V0_trk1_phi, &b_V0_trk1_phi);
    tC->SetBranchAddress("lambda_pion_pt", V0_trk1_pt, &b_V0_trk1_pt);
    tC->SetBranchAddress("lambda_proton_eta", V0_trk2_eta, &b_V0_trk2_eta);
    tC->SetBranchAddress("lambda_proton_phi", V0_trk2_phi, &b_V0_trk2_phi);
    tC->SetBranchAddress("lambda_proton_pt", V0_trk2_pt, &b_V0_trk2_pt);
    tC->SetBranchAddress("lambda_kin_valid", V0_kin_valid, &b_V0_kin_valid);
    tC->SetBranchAddress("lambda_pion_mu_index", V0_trk1_mu_index, &b_V0_trk1_mu_index);
    tC->SetBranchAddress("lambda_proton_mu_index", V0_trk2_mu_index, &b_V0_trk2_mu_index);

    tC->SetBranchAddress("lambda_kin_sipPV", V0_kin_sipPV, &b_V0_kin_sipPV);
    tC->SetBranchAddress("lambda_pion_sip", V0_trk1_sip, &b_V0_trk1_sip);
    tC->SetBranchAddress("lambda_proton_sip", V0_trk2_sip, &b_V0_trk2_sip);

    // if(string::npos != what.find("MC")){
    //   tC->SetBranchAddress("lambda_gen_pdgId", V0_gen_pdgId, &b_V0_gen_pdgId);
    //   tC->SetBranchAddress("ks_gen_trk1_mpdgId",V0_gen_trk1_mpdgId, &b_V0_gen_trk1_mpdgId);
    //   tC->SetBranchAddress("ks_gen_trk2_mpdgId",V0_gen_trk2_mpdgId, &b_V0_gen_trk2_mpdgId);
    //   tC->SetBranchAddress("ks_gen_trk1_pdgId",V0_gen_trk1_pdgId, &b_V0_gen_trk1_pdgId);
    //   tC->SetBranchAddress("ks_gen_trk2_pdgId",V0_gen_trk2_pdgId, &b_V0_gen_trk2_pdgId);
    // }
  }
  else if (string::npos != sample.find("phi"))
  {
    tC->SetBranchAddress("nphi", &nV0, &b_nV0);
    tC->SetBranchAddress("phi_doca", V0_doca, &b_V0_doca);
    tC->SetBranchAddress("phi_kin_cosAlphaXY", V0_kin_cosAlphaXY, &b_V0_kin_cosAlphaXY);
    tC->SetBranchAddress("phi_kin_eta", V0_kin_eta, &b_V0_kin_eta);
    tC->SetBranchAddress("phi_kin_lxy", V0_kin_lxy, &b_V0_kin_lxy);
    tC->SetBranchAddress("phi_kin_mass", V0_kin_mass, &b_V0_kin_mass);
    tC->SetBranchAddress("phi_kin_massErr", V0_kin_massErr, &b_V0_kin_massErr);
    tC->SetBranchAddress("phi_kin_phi", V0_kin_phi, &b_V0_kin_phi);
    tC->SetBranchAddress("phi_kin_pt", V0_kin_pt, &b_V0_kin_pt);
    tC->SetBranchAddress("phi_kin_slxy", V0_kin_slxy, &b_V0_kin_slxy);
    tC->SetBranchAddress("phi_kin_vtx_chi2dof", V0_kin_vtx_chi2dof, &b_V0_kin_vtx_chi2dof);
    tC->SetBranchAddress("phi_kin_vtx_prob", V0_kin_vtx_prob, &b_V0_kin_vtx_prob);
    tC->SetBranchAddress("phi_mass", V0_mass, &b_V0_mass);
    tC->SetBranchAddress("phi_trk1_eta", V0_trk1_eta, &b_V0_trk1_eta);
    tC->SetBranchAddress("phi_trk1_phi", V0_trk1_phi, &b_V0_trk1_phi);
    tC->SetBranchAddress("phi_trk1_pt", V0_trk1_pt, &b_V0_trk1_pt);
    tC->SetBranchAddress("phi_trk2_eta", V0_trk2_eta, &b_V0_trk2_eta);
    tC->SetBranchAddress("phi_trk2_phi", V0_trk2_phi, &b_V0_trk2_phi);
    tC->SetBranchAddress("phi_trk2_pt", V0_trk2_pt, &b_V0_trk2_pt);
    tC->SetBranchAddress("phi_kin_valid", V0_kin_valid, &b_V0_kin_valid);
    tC->SetBranchAddress("phi_trk1_mu_index", V0_trk1_mu_index, &b_V0_trk1_mu_index);
    tC->SetBranchAddress("phi_trk2_mu_index", V0_trk2_mu_index, &b_V0_trk2_mu_index);

    tC->SetBranchAddress("phi_kin_sipPV", V0_kin_sipPV, &b_V0_kin_sipPV);
    tC->SetBranchAddress("phi_trk1_sip", V0_trk1_sip, &b_V0_trk1_sip);
    tC->SetBranchAddress("phi_trk2_sip", V0_trk2_sip, &b_V0_trk2_sip);

    // if(string::npos != what.find("MC")){
    //   tC->SetBranchAddress("phi_gen_pdgId", V0_gen_pdgId, &b_V0_gen_pdgId);
    //   tC->SetBranchAddress("ks_gen_trk1_mpdgId",V0_gen_trk1_mpdgId, &b_V0_gen_trk1_mpdgId);
    //   tC->SetBranchAddress("ks_gen_trk2_mpdgId",V0_gen_trk2_mpdgId, &b_V0_gen_trk2_mpdgId);
    //   tC->SetBranchAddress("ks_gen_trk1_pdgId",V0_gen_trk1_pdgId, &b_V0_gen_trk1_pdgId);
    //   tC->SetBranchAddress("ks_gen_trk2_pdgId",V0_gen_trk2_pdgId, &b_V0_gen_trk2_pdgId);
    // }
  }
  else if (string::npos != sample.find("d0"))
  {
    tC->SetBranchAddress("nd0", &nV0, &b_nV0);
    tC->SetBranchAddress("d0_doca", V0_doca, &b_V0_doca);
    tC->SetBranchAddress("d0_kin_cosAlphaXY", V0_kin_cosAlphaXY, &b_V0_kin_cosAlphaXY);
    tC->SetBranchAddress("d0_kin_eta", V0_kin_eta, &b_V0_kin_eta);
    tC->SetBranchAddress("d0_kin_lxy", V0_kin_lxy, &b_V0_kin_lxy);
    tC->SetBranchAddress("d0_kin_mass", V0_kin_mass, &b_V0_kin_mass);
    tC->SetBranchAddress("d0_kin_massErr", V0_kin_massErr, &b_V0_kin_massErr);
    tC->SetBranchAddress("d0_kin_phi", V0_kin_phi, &b_V0_kin_phi);
    tC->SetBranchAddress("d0_kin_pt", V0_kin_pt, &b_V0_kin_pt);
    tC->SetBranchAddress("d0_kin_slxy", V0_kin_slxy, &b_V0_kin_slxy);
    tC->SetBranchAddress("d0_kin_vtx_chi2dof", V0_kin_vtx_chi2dof, &b_V0_kin_vtx_chi2dof);
    tC->SetBranchAddress("d0_kin_vtx_prob", V0_kin_vtx_prob, &b_V0_kin_vtx_prob);
    tC->SetBranchAddress("d0_mass", V0_mass, &b_V0_mass);
    tC->SetBranchAddress("d0_pion_eta", V0_trk1_eta, &b_V0_trk1_eta);
    tC->SetBranchAddress("d0_pion_phi", V0_trk1_phi, &b_V0_trk1_phi);
    tC->SetBranchAddress("d0_pion_pt", V0_trk1_pt, &b_V0_trk1_pt);
    tC->SetBranchAddress("d0_kaon_eta", V0_trk2_eta, &b_V0_trk2_eta);
    tC->SetBranchAddress("d0_kaon_phi", V0_trk2_phi, &b_V0_trk2_phi);
    tC->SetBranchAddress("d0_kaon_pt", V0_trk2_pt, &b_V0_trk2_pt);
    tC->SetBranchAddress("d0_kin_valid", V0_kin_valid, &b_V0_kin_valid);
    tC->SetBranchAddress("d0_pion_mu_index", V0_trk1_mu_index, &b_V0_trk1_mu_index);
    tC->SetBranchAddress("d0_kaon_mu_index", V0_trk2_mu_index, &b_V0_trk2_mu_index);

    tC->SetBranchAddress("d0_kin_sipPV", V0_kin_sipPV, &b_V0_kin_sipPV);
    tC->SetBranchAddress("d0_pion_sip", V0_trk1_sip, &b_V0_trk1_sip);
    tC->SetBranchAddress("d0_kaon_sip", V0_trk2_sip, &b_V0_trk2_sip);

    // if(string::npos != what.find("MC")){
    //   tC->SetBranchAddress("d0_gen_pdgId", V0_gen_pdgId, &b_V0_gen_pdgId);
    //   tC->SetBranchAddress("ks_gen_trk1_mpdgId",V0_gen_trk1_mpdgId, &b_V0_gen_trk1_mpdgId);
    //   tC->SetBranchAddress("ks_gen_trk2_mpdgId",V0_gen_trk2_mpdgId, &b_V0_gen_trk2_mpdgId);
    //   tC->SetBranchAddress("ks_gen_trk1_pdgId",V0_gen_trk1_pdgId, &b_V0_gen_trk1_pdgId);
    //   tC->SetBranchAddress("ks_gen_trk2_pdgId",V0_gen_trk2_pdgId, &b_V0_gen_trk2_pdgId);
    // }
  }
  else if (string::npos != sample.find("ds"))
  {
    tC->SetBranchAddress("nphi", &nV0, &b_nV0);
    tC->SetBranchAddress("phi_doca", V0_doca, &b_V0_doca);
    tC->SetBranchAddress("phi_ds_cosAlphaXY", V0_kin_cosAlphaXY, &b_V0_kin_cosAlphaXY);
    tC->SetBranchAddress("phi_ds_eta", V0_kin_eta, &b_V0_kin_eta);
    tC->SetBranchAddress("phi_ds_lxy", V0_kin_lxy, &b_V0_kin_lxy);
    tC->SetBranchAddress("phi_ds_mass", V0_kin_mass, &b_V0_kin_mass);
    tC->SetBranchAddress("phi_ds_massErr", V0_kin_massErr, &b_V0_kin_massErr);
    tC->SetBranchAddress("phi_ds_phi", V0_kin_phi, &b_V0_kin_phi);
    tC->SetBranchAddress("phi_ds_pt", V0_kin_pt, &b_V0_kin_pt);
    tC->SetBranchAddress("phi_ds_slxy", V0_kin_slxy, &b_V0_kin_slxy);
    tC->SetBranchAddress("phi_ds_vtx_chi2dof", V0_kin_vtx_chi2dof, &b_V0_kin_vtx_chi2dof);
    tC->SetBranchAddress("phi_ds_vtx_prob", V0_kin_vtx_prob, &b_V0_kin_vtx_prob);
    tC->SetBranchAddress("phi_ds_pion_eta", V0_trk1_eta, &b_V0_trk1_eta);
    tC->SetBranchAddress("phi_ds_pion_phi", V0_trk1_phi, &b_V0_trk1_phi);
    tC->SetBranchAddress("phi_ds_pion_pt", V0_trk1_pt, &b_V0_trk1_pt);
    tC->SetBranchAddress("phi_trk2_eta", V0_trk2_eta, &b_V0_trk2_eta);
    tC->SetBranchAddress("phi_trk2_phi", V0_trk2_phi, &b_V0_trk2_phi);
    tC->SetBranchAddress("phi_trk2_pt", V0_trk2_pt, &b_V0_trk2_pt);
    tC->SetBranchAddress("ds_kin_valid", V0_kin_valid, &b_V0_kin_valid);
    tC->SetBranchAddress("phi_ds_pion_mu_index", V0_trk1_mu_index, &b_V0_trk1_mu_index);
    tC->SetBranchAddress("phi_trk2_mu_index", V0_trk2_mu_index, &b_V0_trk2_mu_index);

    tC->SetBranchAddress("phi_ds_sipPV", V0_kin_sipPV, &b_V0_kin_sipPV);
    tC->SetBranchAddress("phi_ds_sipBS", V0_trk1_sip, &b_V0_trk1_sip);
    tC->SetBranchAddress("phi_trk2_sip", V0_trk2_sip, &b_V0_trk2_sip);
  }
  Long64_t n = tC->GetEntries();
  int count0 = 0;
  int count1 = 0, count2 = 0;
  TString cut_print = "bla";
  TString cut_print_mid = "bla";

  int i_permille_old = 0;
  for (unsigned int i = 0; i < n; ++i)
  {
    int i_permille = (int)floor(100. * i / n);
    if (i_permille != i_permille_old)
    {
      printf("\015\033[32m ---> \033[1m\033[31m%d%%"
             "\033[0m\033[32m <---\033[0m\015",
             i_permille);
      fflush(stdout);
      i_permille_old = i_permille;
    }

    // Load a proper try and get relative even index
    Long64_t localEntry = tC->LoadTree(i);
    b_nV0->GetEntry(localEntry);
    if (nV0 == 0)
      continue;
    count0++;
    b_run->GetEntry(localEntry);
    b_luminosityBlock->GetEntry(localEntry);
    b_V0_kin_mass->GetEntry(localEntry);
    b_V0_doca->GetEntry(localEntry);
    b_V0_kin_cosAlphaXY->GetEntry(localEntry);
    b_V0_kin_eta->GetEntry(localEntry);
    b_V0_kin_lxy->GetEntry(localEntry);
    b_V0_kin_slxy->GetEntry(localEntry);
    b_V0_kin_phi->GetEntry(localEntry);
    b_V0_kin_pt->GetEntry(localEntry);
    b_V0_kin_vtx_chi2dof->GetEntry(localEntry);
    b_V0_kin_vtx_prob->GetEntry(localEntry);
    b_V0_trk1_eta->GetEntry(localEntry);
    b_V0_trk2_eta->GetEntry(localEntry);
    b_V0_trk2_pt->GetEntry(localEntry);
    b_V0_trk1_pt->GetEntry(localEntry);
    b_V0_trk1_phi->GetEntry(localEntry);
    b_V0_trk2_phi->GetEntry(localEntry);
    b_V0_trk1_mu_index->GetEntry(localEntry);
    b_V0_trk2_mu_index->GetEntry(localEntry);
    b_Muon_softMvaId->GetEntry(localEntry);
    b_Muon_softMva->GetEntry(localEntry);
    b_V0_kin_sipPV->GetEntry(localEntry);
    b_V0_trk1_sip->GetEntry(localEntry);
    b_V0_trk2_sip->GetEntry(localEntry);
    b_Muon_softId->GetEntry(localEntry);
    b_Muon_mediumId->GetEntry(localEntry);
    b_Muon_isGlobal->GetEntry(localEntry);
    b_Muon_isTracker->GetEntry(localEntry);
    b_Muon_looseId->GetEntry(localEntry);
    b_nMuon->GetEntry(localEntry);
    b_Muon_pt->GetEntry(localEntry);
    if (string::npos != sample.find("trigger"))
    {
      b_HLT_Ele30->GetEntry(localEntry);
    }
    if (string::npos != sample.find("DoubleMuMu"))
    {
      b_MuonId_hlt_pt->GetEntry(localEntry);
    }
    if (string::npos != what.find("MC"))
    {
      b_V0_gen_pdgId->GetEntry(localEntry);
      b_V0_gen_trk1_mpdgId->GetEntry(localEntry);
      b_V0_gen_trk2_mpdgId->GetEntry(localEntry);
      b_V0_gen_trk1_pdgId->GetEntry(localEntry);
      b_V0_gen_trk2_pdgId->GetEntry(localEntry);
    }
    for (unsigned int bhad = 0; bhad < nV0; ++bhad)
    {
      count1++;

      bool Finalcut_bin0, Finalcut_bin1, Finalcut_bin2, Finalcut_bin3;
      bool Finalcut_bin4, Finalcut_bin5, Finalcut_bin6, Finalcut_bin7;
      if (binning == "pT") {
         Finalcut_bin0 = V0_kin_pt[bhad] > 0.0 && V0_kin_pt[bhad] < 4.0;
         Finalcut_bin1 = V0_kin_pt[bhad] > 4.0 && V0_kin_pt[bhad] < 8.0;
         Finalcut_bin2 = V0_kin_pt[bhad] > 8.0 && V0_kin_pt[bhad] < 12.0;
         Finalcut_bin3 = V0_kin_pt[bhad] > 12.0 && V0_kin_pt[bhad] < 16.0;
         Finalcut_bin4 = V0_kin_pt[bhad] > 16.0 && V0_kin_pt[bhad] < 20.0;
         Finalcut_bin5 = V0_kin_pt[bhad] > 20.0 && V0_kin_pt[bhad] < 30.0;
         Finalcut_bin6 = V0_kin_pt[bhad] > 30.0 && V0_kin_pt[bhad] < 50.0;
         Finalcut_bin7 = V0_kin_pt[bhad] > 0.0 && V0_kin_pt[bhad] < 50.0;
      } else {
         Finalcut_bin0 = abs(V0_kin_lxy[bhad]) > 0.0 &&  abs(V0_kin_lxy[bhad]) < 1.;
         Finalcut_bin1 = abs(V0_kin_lxy[bhad]) > 1.0 && abs(V0_kin_lxy[bhad]) < 2.;
         Finalcut_bin2 = abs(V0_kin_lxy[bhad]) > 2.0 && abs(V0_kin_lxy[bhad]) < 4.0;
         Finalcut_bin3 = abs(V0_kin_lxy[bhad]) > 4.0 && abs(V0_kin_lxy[bhad]) < 8.0;
         Finalcut_bin4 = abs(V0_kin_lxy[bhad]) > 8.0 && abs(V0_kin_lxy[bhad]) < 12.;
         Finalcut_bin5 = abs(V0_kin_lxy[bhad]) > 12. && abs(V0_kin_lxy[bhad]) < 16;
         Finalcut_bin6 = abs(V0_kin_lxy[bhad]) > 16.0 && abs(V0_kin_lxy[bhad]) < 40.0;
         Finalcut_bin7 =abs(V0_kin_lxy[bhad]) > 0.0 && abs(V0_kin_lxy[bhad]) < 40.0;
      }

      bool cut_pre = false;
      bool cut_pre_softid = false;
      bool cut_pre_medid = false;
      bool cut_pre_mid = false;
      bool cut_pre_mid20 = false;
      bool cut_pre_mid30 = false;
      bool cut_pre_mid40 = false;
      bool cut_pre_mid45 = false;
      bool cut_pre_mid50 = false;
      bool cut_pre_mid55 = false;
      bool cut_pre_mid60 = false;

      int lxy_cut =  100;

      if (string::npos != sample.find("ks") && string::npos != sample.find("DoubleMuMu"))
      {

        cut_pre = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 2);

        cut_pre_softid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && ((Muon_isGlobal[V0_trk1_mu_index[bhad]] && Muon_isTracker[V0_trk1_mu_index[bhad]] && Muon_looseId[V0_trk1_mu_index[bhad]]) || (Muon_isGlobal[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]] && Muon_looseId[V0_trk2_mu_index[bhad]])));

        cut_pre_medid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_mediumId[V0_trk1_mu_index[bhad]] || Muon_mediumId[V0_trk2_mu_index[bhad]]));

        cut_pre_mid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMvaId[V0_trk1_mu_index[bhad]] || Muon_softMvaId[V0_trk2_mu_index[bhad]]));

        cut_pre_mid20 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMva[V0_trk1_mu_index[bhad]] > 0.20 || Muon_softMva[V0_trk2_mu_index[bhad]] > 0.20));

        cut_pre_mid30 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMva[V0_trk1_mu_index[bhad]] > 0.30 || Muon_softMva[V0_trk2_mu_index[bhad]] > 0.30));

        cut_pre_mid40 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMva[V0_trk1_mu_index[bhad]] > 0.40 || Muon_softMva[V0_trk2_mu_index[bhad]] > 0.40));

        cut_pre_mid45 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMva[V0_trk1_mu_index[bhad]] > 0.45 || Muon_softMva[V0_trk2_mu_index[bhad]] > 0.45));

        cut_pre_mid50 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMva[V0_trk1_mu_index[bhad]] > 0.50 || Muon_softMva[V0_trk2_mu_index[bhad]] > 0.50));

        cut_pre_mid55 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMva[V0_trk1_mu_index[bhad]] > 0.50 || Muon_softMva[V0_trk2_mu_index[bhad]] > 0.55));

        cut_pre_mid60 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && nMuon == 3 && (((V0_trk1_mu_index[bhad] == 0 ^ V0_trk2_mu_index[bhad] == 0) && MuonId_hlt_pt[1] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 1 ^ V0_trk2_mu_index[bhad] == 1) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[2] > 0) || ((V0_trk1_mu_index[bhad] == 2 ^ V0_trk2_mu_index[bhad] == 2) && MuonId_hlt_pt[0] > 0 && MuonId_hlt_pt[1] > 0)) && (Muon_softMva[V0_trk1_mu_index[bhad]] > 0.50 || Muon_softMva[V0_trk2_mu_index[bhad]] > 0.60));
      }
      /*      else if(string::npos != sample.find("ks") && string::npos != what.find("MC")){


        cut_pre =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && (((V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) && abs(V0_gen_trk1_mpdgId[bhad])==311) || ((V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) && abs(V0_gen_trk2_mpdgId[bhad])==311)));

        cut_pre_softid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_isGlobal[V0_trk1_mu_index[bhad]] && Muon_isTracker[V0_trk1_mu_index[bhad]] && Muon_looseId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_looseId[V0_trk2_mu_index[bhad]] && Muon_isGlobal[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

        cut_pre_medid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

        cut_pre_mid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

        cut_pre_mid20 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.20) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.20)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
        cut_pre_mid30 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.30) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.30)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
        cut_pre_mid40 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.40) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.40)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
        cut_pre_mid45 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.45) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.45)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
        cut_pre_mid50 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.50) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.50)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
        cut_pre_mid55 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.55) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.55)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
        cut_pre_mid60 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.60) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.60)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

            }
            else if(string::npos != sample.find("ks") && string::npos == what.find("MC")){
     */
      else if (string::npos != sample.find("ks"))
      {

        cut_pre = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);

        cut_pre_softid = V0_gen_trk1_pdgId[bhad] != 22 && (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_isGlobal[V0_trk1_mu_index[bhad]] && Muon_isTracker[V0_trk1_mu_index[bhad]] && Muon_looseId[V0_trk1_mu_index[bhad]]) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_looseId[V0_trk2_mu_index[bhad]] && Muon_isGlobal[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);

        cut_pre_medid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk1_mu_index[bhad]]) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);

        cut_pre_mid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk1_mu_index[bhad]]) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);

        cut_pre_mid20 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.20) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.20)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);
        cut_pre_mid30 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.30) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.30)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);
        cut_pre_mid40 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.40) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.40)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);
        cut_pre_mid45 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.45) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.45)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);
        cut_pre_mid50 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.50) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.50)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);
        cut_pre_mid55 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.55) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.55)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);
        cut_pre_mid60 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.60) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.60)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3 && HLT_Ele30);
      }
      else if (string::npos != sample.find("phi"))
      {
        cut_pre = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        // cut_pre_mid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]])) &&  V0_kin_vtx_prob[bhad]>0.3 && V0_kin_sipPV[bhad]<1&& V0_doca[bhad]<0.004 && V0_kin_lxy[bhad]<4);

        cut_pre_softid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_looseId[V0_trk1_mu_index[bhad]] && Muon_isTracker[V0_trk1_mu_index[bhad]] && Muon_isGlobal[V0_trk1_mu_index[bhad]]) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_looseId[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]] && Muon_isGlobal[V0_trk2_mu_index[bhad]])) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_medid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk1_mu_index[bhad]]) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]])) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk1_mu_index[bhad]]) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]])) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid20 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.20) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.20)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid30 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.30) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.30)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid40 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.40) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.40)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid45 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.45) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.45)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid50 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.50) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.50)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid55 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.55) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.55)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
        cut_pre_mid60 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]] > 0.60) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.60)) && V0_kin_vtx_prob[bhad] > 0.3 && V0_kin_sipPV[bhad] < 1 && V0_doca[bhad] < 0.004 && V0_kin_lxy[bhad] < 4);
      }
      else if (string::npos != sample.find("lambda"))
      {
        cut_pre = V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2;
        cut_pre_softid = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_looseId[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]] && Muon_isGlobal[V0_trk2_mu_index[bhad]] && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);

        cut_pre_medid = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]] && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        cut_pre_mid = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]] && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        // cut_pre_mid =( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]] && V0_kin_slxy[bhad]>3 &&  V0_kin_sipPV[bhad]<1 && V0_trk1_sip[bhad]>2 && V0_trk2_sip[bhad]>2);

        cut_pre_mid20 = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.20 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        cut_pre_mid30 = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.30 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        cut_pre_mid40 = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.40 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        cut_pre_mid45 = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.45 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        cut_pre_mid50 = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.50 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        cut_pre_mid55 = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.55 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
        cut_pre_mid60 = (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]] > 0.60 && V0_kin_slxy[bhad] > 3 && V0_kin_sipPV[bhad] < 1 && V0_trk1_sip[bhad] > 2 && V0_trk2_sip[bhad] > 2);
      }
      // else if(string::npos != sample.find("lambda")){
      // 	 cut_pre = V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 ;//&& V0_kin_slxy[bhad]>3 &&  V0_kin_sipPV[bhad]<1 && V0_trk1_sip[bhad]>2 && V0_trk2_sip[bhad]>2;
      // 	 cut_pre_mid =( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]]) ;//&& V0_kin_slxy[bhad]>3 &&  V0_kin_sipPV[bhad]<1 && V0_trk1_sip[bhad]>2 && V0_trk2_sip[bhad]>2);
      // 	 //cut_pre_mid =( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]] && V0_kin_slxy[bhad]>3 &&  V0_kin_sipPV[bhad]<1 && V0_trk1_sip[bhad]>2 && V0_trk2_sip[bhad]>2);

      // }

      if (cut_pre)
        count2++;
      for (int k = 0; k < NREG_pt; k++)
      {

        bool cut_tmp = false;
        if (k == 0)
          cut_tmp = Finalcut_bin0;
        if (k == 1)
          cut_tmp = Finalcut_bin1;
        if (k == 2)
          cut_tmp = Finalcut_bin2;
        if (k == 3)
          cut_tmp = Finalcut_bin3;
        if (k == 4)
          cut_tmp = Finalcut_bin4;
        if (k == 5)
          cut_tmp = Finalcut_bin5;
        if (k == 6)
          cut_tmp = Finalcut_bin6;
        if (k == 7)
          cut_tmp = Finalcut_bin7;
        if (cut_tmp)
        {
          if (binning=="pT") {
              if (cut_pre)
            {
              hv0_kin_pt_bin[k]->Fill(V0_kin_pt[bhad]);
              hv0_Mass_bin[k]->Fill(V0_kin_mass[bhad]);
            }
            if (cut_pre_softid)
            {

              hv0_Mass_bin_softid[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_softid[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_medid)
            {

              hv0_Mass_bin_medid[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_medid[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid)
            {

              hv0_Mass_bin_muid[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid20)
            {
              hv0_Mass_bin_muid_20[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_20[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid30)
            {
              hv0_Mass_bin_muid_30[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_30[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid40)
            {
              hv0_Mass_bin_muid_40[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_40[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid45)
            {
              hv0_Mass_bin_muid_45[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_45[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid50)
            {
              hv0_Mass_bin_muid_50[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_50[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid55)
            {
              hv0_Mass_bin_muid_55[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_55[k]->Fill(V0_kin_pt[bhad]);
            }
            if (cut_pre_mid60)
            {
              hv0_Mass_bin_muid_60[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_60[k]->Fill(V0_kin_pt[bhad]);
            }
          } else {
            if (cut_pre)
            {
              hv0_kin_pt_bin[k]->Fill(V0_kin_lxy[bhad]);
              hv0_Mass_bin[k]->Fill(V0_kin_mass[bhad]);
            }
            if (cut_pre_softid)
            {

              hv0_Mass_bin_softid[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_softid[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_medid)
            {

              hv0_Mass_bin_medid[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_medid[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid)
            {

              hv0_Mass_bin_muid[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid20)
            {
              hv0_Mass_bin_muid_20[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_20[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid30)
            {
              hv0_Mass_bin_muid_30[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_30[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid40)
            {
              hv0_Mass_bin_muid_40[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_40[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid45)
            {
              hv0_Mass_bin_muid_45[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_45[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid50)
            {
              hv0_Mass_bin_muid_50[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_50[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid55)
            {
              hv0_Mass_bin_muid_55[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_55[k]->Fill(V0_kin_lxy[bhad]);
            }
            if (cut_pre_mid60)
            {
              hv0_Mass_bin_muid_60[k]->Fill(V0_kin_mass[bhad]);
              hv0_kin_pt_bin_muid_60[k]->Fill(V0_kin_lxy[bhad]);
            }
          }
        }
      } // loop for pt bins
    }
  }

  string mode = "sideband";
  if ((string::npos != what.find("MC")))
    mode = "MC";

  cout << " count " << count0 << "\t" << count1 << "\t" << count2 << endl;

  cout << " pt " << endl;
  cout << hv0_kin_pt_bin[0]->Integral() << "\t" << hv0_kin_pt_bin[1]->Integral() << endl;
  cout << hv0_kin_pt_bin[2]->Integral() << "\t" << hv0_kin_pt_bin[3]->Integral() << endl;
  cout << hv0_kin_pt_bin[4]->Integral() << "\t" << hv0_kin_pt_bin[5]->Integral() << endl;

  cout << "muid  pt " << endl;
  cout << hv0_kin_pt_bin_muid[0]->Integral() << "\t" << hv0_kin_pt_bin_muid[1]->Integral() << endl;
  cout << hv0_kin_pt_bin_muid[2]->Integral() << "\t" << hv0_kin_pt_bin_muid[3]->Integral() << endl;
  cout << hv0_kin_pt_bin_muid[4]->Integral() << "\t" << hv0_kin_pt_bin_muid[5]->Integral() << endl;

  cout << " mass pt " << endl;
  cout << hv0_Mass_bin[0]->Integral() << "\t" << hv0_Mass_bin[1]->Integral() << endl;
  cout << hv0_Mass_bin[2]->Integral() << "\t" << hv0_Mass_bin[3]->Integral() << endl;
  cout << hv0_Mass_bin[4]->Integral() << "\t" << hv0_Mass_bin[5]->Integral() << endl;

  cout << "muid  pt " << endl;
  cout << hv0_Mass_bin_muid[0]->Integral() << "\t" << hv0_Mass_bin_muid[1]->Integral() << endl;
  cout << hv0_Mass_bin_muid[2]->Integral() << "\t" << hv0_Mass_bin_muid[3]->Integral() << endl;
  cout << hv0_Mass_bin_muid[4]->Integral() << "\t" << hv0_Mass_bin_muid[5]->Integral() << endl;

  cout << "muid  pt mva30 " << endl;
  cout << hv0_kin_pt_bin_muid_30[0]->Integral() << "\t" << hv0_kin_pt_bin_muid_30[1]->Integral() << endl;
  cout << hv0_kin_pt_bin_muid_30[2]->Integral() << "\t" << hv0_kin_pt_bin_muid_30[3]->Integral() << endl;
  cout << hv0_kin_pt_bin_muid_30[4]->Integral() << "\t" << hv0_kin_pt_bin_muid_30[5]->Integral() << endl;
  cout << "muid Mass mva30 " << endl;
  cout << hv0_Mass_bin_muid_30[0]->Integral() << "\t" << hv0_Mass_bin_muid_30[1]->Integral() << endl;
  cout << hv0_Mass_bin_muid_30[2]->Integral() << "\t" << hv0_Mass_bin_muid_30[3]->Integral() << endl;
  cout << hv0_Mass_bin_muid_30[4]->Integral() << "\t" << hv0_Mass_bin_muid_30[5]->Integral() << endl;

  cout << "muid  pt softid " << endl;
  cout << hv0_kin_pt_bin_softid[0]->Integral() << "\t" << hv0_kin_pt_bin_softid[1]->Integral() << endl;
  cout << hv0_kin_pt_bin_softid[2]->Integral() << "\t" << hv0_kin_pt_bin_softid[3]->Integral() << endl;
  cout << hv0_kin_pt_bin_softid[4]->Integral() << "\t" << hv0_kin_pt_bin_softid[5]->Integral() << endl;
  cout << "muid Mass softid " << endl;
  cout << hv0_Mass_bin_softid[0]->Integral() << "\t" << hv0_Mass_bin_softid[1]->Integral() << endl;
  cout << hv0_Mass_bin_softid[2]->Integral() << "\t" << hv0_Mass_bin_softid[3]->Integral() << endl;
  cout << hv0_Mass_bin_softid[4]->Integral() << "\t" << hv0_Mass_bin_softid[5]->Integral() << endl;

  TString new_filename = Form("%s/histo_%s.root", path.c_str(), sample.c_str());
  TFile *output = TFile::Open(new_filename.Data(), "RECREATE");
  for (int k(0); k < NREG_pt; k++)
  {
    hv0_kin_pt_bin[k]->Write();
    hv0_Mass_bin[k]->Write();

    hv0_kin_pt_bin_softid[k]->Write();
    hv0_Mass_bin_softid[k]->Write();
    hv0_kin_pt_bin_medid[k]->Write();
    hv0_Mass_bin_medid[k]->Write();
    hv0_kin_pt_bin_muid[k]->Write();
    hv0_Mass_bin_muid[k]->Write();

    hv0_kin_pt_bin_muid_30[k]->Write();
    hv0_kin_pt_bin_muid_40[k]->Write();
    hv0_kin_pt_bin_muid_45[k]->Write();
    hv0_kin_pt_bin_muid_20[k]->Write();
    hv0_kin_pt_bin_muid_50[k]->Write();
    hv0_kin_pt_bin_muid_55[k]->Write();
    hv0_kin_pt_bin_muid_60[k]->Write();

    hv0_Mass_bin_muid_30[k]->Write();
    hv0_Mass_bin_muid_40[k]->Write();
    hv0_Mass_bin_muid_45[k]->Write();
    hv0_Mass_bin_muid_20[k]->Write();
    hv0_Mass_bin_muid_50[k]->Write();
    hv0_Mass_bin_muid_55[k]->Write();
    hv0_Mass_bin_muid_60[k]->Write();
  }
  output->Write();
  output->Close();
  for (int k = 0; k < NREG_pt; k++)
  {

    hv0_kin_pt_bin[k]->Delete();
    hv0_Mass_bin[k]->Delete();
    hv0_kin_pt_bin_softid[k]->Delete();
    hv0_Mass_bin_softid[k]->Delete();
    hv0_kin_pt_bin_medid[k]->Delete();
    hv0_Mass_bin_medid[k]->Delete();
    hv0_kin_pt_bin_muid[k]->Delete();
    hv0_Mass_bin_muid[k]->Delete();
    hv0_kin_pt_bin_muid_30[k]->Delete();
    hv0_kin_pt_bin_muid_40[k]->Delete();
    hv0_kin_pt_bin_muid_45[k]->Delete();
    hv0_kin_pt_bin_muid_20[k]->Delete();
    hv0_kin_pt_bin_muid_50[k]->Delete();
    hv0_kin_pt_bin_muid_55[k]->Delete();
    hv0_kin_pt_bin_muid_60[k]->Delete();

    hv0_Mass_bin_muid_30[k]->Delete();
    hv0_Mass_bin_muid_40[k]->Delete();
    hv0_Mass_bin_muid_45[k]->Delete();
    hv0_Mass_bin_muid_20[k]->Delete();
    hv0_Mass_bin_muid_50[k]->Delete();
    hv0_Mass_bin_muid_55[k]->Delete();
    hv0_Mass_bin_muid_60[k]->Delete();
  }
  cout << "complete" << endl;

  if (string::npos != sample.find("ks"))
  {
    cut_print = " (((ks_trk1_pt[bhad] > 4 && ks_trk2_pt[bhad] > 1.0) || ( ks_trk2_pt[bhad] > 4 && ks_trk1_pt[bhad] > 1.0)) && ks_kin_slxy[bhad]>15 && ks_kin_vtx_chi2dof[bhad] < 3 && ks_kin_cosAlphaXY[bhad]>0.999 && ks_doca[bhad]<0.004 &&  ks_kin_sipPV[bhad]<3 && ks_trk1_sip[bhad]>5 && ks_trk2_sip[bhad]>5)";
    cut_print_mid = "(((ks_trk1_pt[bhad] > 4 && ks_trk2_pt[bhad] > 1.0 && ks_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[ks_trk1_mu_index[bhad]]) || ( ks_trk2_pt[bhad] > 4 && ks_trk1_pt[bhad] > 1.0 && ks_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[ks_trk2_mu_index[bhad]])) && ks_kin_slxy[bhad]>15 && ks_kin_vtx_chi2dof[bhad] < 3 && ks_kin_cosAlphaXY[bhad]>0.999 && ks_doca[bhad]<0.004 &&  ks_kin_sipPV[bhad]<3 && ks_trk1_sip[bhad]>5 && ks_trk2_sip[bhad]>5)";
  }
  if (string::npos != sample.find("phi"))
  {
    cut_print = "(((phi_trk1_pt[bhad] > 4 && phi_trk2_pt[bhad] > 3.0) || ( phi_trk2_pt[bhad] > 4 && phi_trk1_pt[bhad] > 3.0)) &&  phi_kin_vtx_prob[bhad]>0.3 && phi_kin_sipPV[bhad]<1  && phi_doca[bhad]<0.004 && phi_kin_lxy[bhad]<4)";
    cut_print_mid = "(((phi_trk1_pt[bhad] > 4 && phi_trk2_pt[bhad] > 1.0 && phi_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[phi_trk1_mu_index[bhad]]) || ( phi_trk2_pt[bhad] > 4 && phi_trk1_pt[bhad] > 1.0 && phi_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[phi_trk2_mu_index[bhad]])) &&  phi_kin_vtx_prob[bhad]>0.3 && phi_kin_sipPV[bhad]<1&& phi_doca[bhad]<0.004 && phi_kin_lxy[bhad]<4)";
  }
  if (string::npos != sample.find("lambda"))
  {

    // cut_print ="((lambda_pion_pt[bhad] > 4 && lambda_proton_pt[bhad] > 2.0) || (
    cut_print = "lambda_proton_pt[bhad] > 4 && lambda_pion_pt[bhad] > 1.0 && lambda_kin_slxy[bhad]>2 &&  lambda_kin_sipPV[bhad]<1 && lambda_pion_sip[bhad]>2 && lambda_proton_sip[bhad]>2";
    cut_print_mid = "( lambda_proton_pt[bhad] > 4 && lambda_pion_pt[bhad] > 1.0 && lambda_proton_mu_index[bhad] >= 0 && Muon_softMvaId[lambda_proton_mu_index[bhad]] && V0_kin_slxy[bhad]>2 &&  V0_kin_sipPV[bhad]<3 && lambda_pion_sip[bhad]>2 && lambda_proton_sip[bhad]>2)";
  }
  if (string::npos != sample.find("d0"))
  {
    cut_print = "(((d0_pion_pt[bhad] > 4 && d0_kaon_pt[bhad] > 3.0) || ( d0_kaon_pt[bhad] > 4 && d0_pion_pt[bhad] > 3.0)) && d0_kin_sipPV[bhad]<1 )";
    cut_print_mid = "(((d0_pion_pt[bhad] > 4 && d0_kaon_pt[bhad] > 1.0 && d0_pion_mu_index[bhad] >= 0 && Muon_softMvaId[d0_pion_mu_index[bhad]]) || ( d0_kaon_pt[bhad] > 4 && d0_pion_pt[bhad] > 1.0 && d0_kaon_mu_index[bhad] >= 0 && Muon_softMvaId[d0_kaon_mu_index[bhad]])) && V0_kin_sipPV[bhad]<1 )";
    // cut_print_mid = "(( d0_pion_mu_index[bhad] >= 0 && Muon_softMvaId[d0_pion_mu_index[bhad]])||(d0_kaon_mu_index[bhad] >= 0 && Muon_softMvaId[d0_kaon_mu_index[bhad]]))";
  }

  TCanvas *c = new TCanvas("c", "", 2000, 2000);
  TLatex latex;
  latex.SetTextSize(0.025);
  latex.SetTextAlign(13); // align at top
  latex.DrawLatex(.2, .9, "Preselection");
  int check_length = cut_print.Length();
  int number_ = ceil(check_length / 50);
  double pos_ = 0.83;
  for (int loop = 0; loop <= number_; loop++)
  {

    TString tmp_print = cut_print(loop * 50, loop + 50);
    latex.DrawLatex(.1, pos_, tmp_print);
    pos_ = pos_ - 0.03;
  }
  latex.SetTextAlign(12); // centered
  latex.DrawLatex(.2, .6, "Preselection + muid");

  int check_length_mid = cut_print_mid.Length();
  int number_mid = ceil(check_length_mid / 50);
  double pos_mid = 0.53;
  for (int loop = 0; loop <= number_mid; loop++)
  {

    TString tmp_print = cut_print_mid(loop * 50, loop + 50);
    latex.DrawLatex(.1, pos_mid, tmp_print);
    pos_mid = pos_mid - 0.03;
  }

  // latex.DrawLatex(.2,.5,cut_print_mid);
  c->SaveAs(Form("%s_cuts.png", path.c_str()));

  return;
}

// ----------------------------------------------------------------------
void plotFakenew::loopOverChain(TChain *tC, string basic = "Data", string modf = "noprod", string path = "ks/", TString binning = "pT")
{
  cout << "binning is " << binning << endl;
  // if(string::npos != modf.find("Data")) sbsDistributions_bin(tC, cutsks, "ks", modf);
  if (string::npos != modf.find("Histproduction"))
    sbsDistributions(tC, basic, modf, path, binning);
  // playfit(basic, path);
  //  playPhifit(basic, path);

  cout << "lopping overchain function closing" << endl;
}

void plotFakenew::list_dir(TChain *tC, const char *path)
{
  struct dirent *entry;
  DIR *dir = opendir(path);
  if (dir == NULL)
  {
    cout << "No file is in this drectory" << endl;
    return;
  }
  while ((entry = readdir(dir)) != NULL)
  {
    string filenamec = entry->d_name;
    string s1 = ".root";

    if (filenamec.find(s1) != string::npos)
    {
      TString filename_f = path + filenamec;
      // std::cout << "found! " << filename_f.c_str()<< '\n';
      TFile *file = TFile::Open(filename_f, "READ");
      if (!file)
        continue; // file does not exist
      if (file->IsZombie())
        continue; // file is unusable
      if (file->TestBit(TFile::kRecovered))
        continue; // file has been recovered
      delete file;
      string check_mu = "MuEnrichedPt";
      string check_HT = "ZH";
      string check_zz = "TTT";
      string data_ = "Run";
      string modename = Form("%s", filename_f.Data());
      // if((modename.find(check_mu) == string::npos) || (modename.find(check_HT) == string::npos) || (modename.find(check_zz) == string::npos)){
      tC->Add(filename_f);
      // cout<<modename<<endl;
      // }
      // if((modename.find(data_) != string::npos)){
      // 	cout<<modename<<endl;
      // 	tC->Add(filename_f);
      // }
      // tC->Add(filename_f);
      // cout<<"    checking entry "<<tC->GetEntries()<<endl;
    }
  }
  closedir(dir);
}

void plotFakenew::list_dir_file(TChain *tC, string filename)
{
  fstream newfile;
  newfile.open(filename, ios::in); // open a file to perform read operation using file object
  if (newfile.is_open())
  { // checking whether the file is open
    string tp;
    while (getline(newfile, tp))
    { // read data from file object and put it into string.
      // cout << tp << "\n"; //print the data of the string
      const char *c = tp.c_str();
      list_dir(tC, c);
    }
  }
}
int main()
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(1, false);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(Optimization);
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().Print();

  string path = "ks/";
  TString binning = "lxy";

  plotFakenew c1;
  TChain *tC1 = new TChain("Events");
  c1.list_dir_file(tC1, "2022_MC_DY.txt");
  cout << "Entry Data " << tC1->GetEntries() << endl;
  c1.loopOverChain(tC1, Form("MC_DY_522_2022_ks_%s", binning.Data()), "Histproduction", path, binning);
  delete tC1;

  plotFakenew c2;
  TChain* tC2 = new TChain("Events");
  c2.list_dir_file(tC2, "2022_data_ks_egamma.txt");
  cout<<"Entry Data"<<tC2->GetEntries()<<endl;
  c2.loopOverChain(tC2, Form("data_egamma_522_2022_ks_trigger_%s", binning.Data()), "Histproduction", path, binning);
  delete tC2;

  plotFakenew c3;
  TChain* tC3 = new TChain("Events");
  c3.list_dir_file(tC3, "2022_data_ks_parking.txt");
  cout<<"Entry Data"<<tC3->GetEntries()<<endl;
  c3.loopOverChain(tC3, Form("data_parking_522_2022_ks_%s", binning.Data()),"Histproduction", path, binning);
  delete tC3;

  plotFakenew c4;
  TChain* tC4 = new TChain("Events");
  c4.list_dir_file(tC4, "2022_MC_TT.txt");
  cout<<"Entry Data"<<tC4->GetEntries()<<endl;
  c4.loopOverChain(tC4, Form("MC_TT_522_2022_ks_%s", binning.Data()), "Histproduction", path, binning);
  delete tC4;

  plotFakenew c5;
  TChain* tC5 = new TChain("Events");
  c5.list_dir_file(tC5, "2022_MC_W.txt");
  cout<<"Entry Data"<<tC5->GetEntries()<<endl;
  c5.loopOverChain(tC5, Form("MC_W_522_2022_ks_%s", binning.Data()), "Histproduction", path, binning);
  delete tC5;

  plotFakenew c6;
  TChain* tC6 = new TChain("Events");
  c6.list_dir_file(tC6, "2022_MC_combined.txt");
  cout<<"Entry Data"<<tC6->GetEntries()<<endl;
  c6.loopOverChain(tC6, Form("MC_combined_522_2022_ks_%s", binning.Data()), "Histproduction", path, binning);
  delete tC6; 

  /*plotFakenew c7;
  TChain* tC7 = new TChain("Events");
  c7.list_dir_file(tC7, "2022_data_ks_DoubleMuMu.txt");
  cout<<"Entry Data"<<tC7->GetEntries()<<endl;
  c7.loopOverChain(tC7, Form("data_DoubleMuMu_521_2022_ks_%s", binning.Data()), "Histproduction", path, binning);
  */
  
  return 0;
}
