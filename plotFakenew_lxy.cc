
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
class plotFakenew_lxy {
public:
void initPol1(double &, double &, TH1 *);
  //void applyLimits(TF1 *, string );
  TF1* pol1gauss2c(TH1 *, double , double );
  TF1* phiKK(TH1 *);
  double dRatio(double a, double ae, double b, double be) {
    return TMath::Sqrt(((ae*ae)/(b*b)) + ((a*a*be*be)/(b*b*b*b)));
  }
  //  vector<double> fitV0Roofit(TH1D* , string what = "nothing", string path="/eos/home-c/ckar/");  
  
  //void fitV0(TH1D*, string what = "nothing", string path="/eos/home-c/ckar/");
  //vector<double> fitV0phikk(TH1D* , string what = "nothing", string path="/eos/home-c/ckar/");
  void sbsDistributions(TChain* , string , string what = "nothing", string path="ks/");      
  void playfit(string, string path="ks/"); 
  //vector<double> fitV0Roofit_altphi(TH1D* , TH1D*, string what = "nothing", string path="ks/");  
  void playPhifit(string, string path="ks/"); 
  void OverlayMCDATA(string ,string path="ks/");
  void FAKERATE(string , string path="ks/" );
  void list_dir(TChain * , const char *);
  void list_dir_file(TChain * , string ob="bla");
  void loopOverChain(TChain *, string , string, string);
  plotFakenew_lxy();
  ~plotFakenew_lxy();
  
private:

  bool fVerbose =true;
  int retry = 3 ;
  
  std::map<string, std::vector<TH1D*> >  hv0_kin_pt_Fake;
  std::map<string, std::vector<TH1D*> >  hv0_kin_eta_Fake;
  TH1D* hv0_kin_pt_Fake_D;
  TH1D* hv0_kin_pt_Fake_M;
};


 void plotFakenew_lxy::initPol1(double &p0, double &p1, TH1 *h) {
   int EDG(4), NB(EDG+1);
   int lbin(1), hbin(h->GetNbinsX()+1);
   double fLo =99.0, fHi =-99.0; 
   if (fLo < fHi) {
     lbin = h->FindBin(fLo);
     hbin = h->FindBin(fHi);
   }
   double xlo = h->GetBinLowEdge(lbin);

   double dx = h->GetBinLowEdge(hbin) - xlo;
   double ylo = h->Integral(lbin, lbin+EDG)/NB;
   double yhi = h->Integral(hbin-EDG-1, hbin-1)/NB;

   p1  = (yhi-ylo)/dx;
   p0  = ylo - p1*xlo;
   if (1 || fVerbose) {
     cout << "lbin: " << lbin << " hbin: " << hbin
	  << " ylo: " << ylo << " yhi: " << yhi << " dx: " << dx
	  << " p0: " << p0 << " p1: " << p1 << endl;
   }
 }
 // ----------------------------------------------------------------------                                                                                                    
 double iF_pol1(double *x, double *par) {
   return par[0] + par[1]*x[0];
 }
 plotFakenew_lxy::~plotFakenew_lxy(void){
   //destructor
 }
 plotFakenew_lxy::plotFakenew_lxy(void){
   //constructor
 }

 double iF_gauss2c(double *x, double *par) {
   // constrained to have the same mean in the second gaussian                             
   // par[0] -> const                                                                      
   // par[1] -> mean                                                                       
   // par[2] -> sigma                                                                      
   // par[3] -> fraction in second gaussian                                                
   // par[4] -> sigma of second gaussian                                                   
   Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
   if (par[2] > 0) {
     arg1 = (x[0] - par[1]) / par[2];
     fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
   }
   if (par[4] > 0.) {
     arg2 = (x[0] - par[1]) / par[4];
     fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
   }
   Double_t fitval = fitval1 + fitval2;
   return fitval;
 }

 double iF_pol1_gauss2c(double *x, double *par) {
   //   par[0] = norm of gaussian                                                          
   //   par[1] = mean of gaussian                                                          
   //   par[2] = sigma of gaussian                                                         
   //   par[3] = fraction in second gaussian                                               
   //   par[4] = sigma of gaussian                                                         
   //   par[5] = par 0 of pol1                                                             
   //   par[6] = par 1 of pol1                                                             
   return  (iF_pol1(x, &par[5]) + iF_gauss2c(x, &par[0]));
 }


 TF1* plotFakenew_lxy::pol1gauss2c(TH1 *h, double peak, double sigma) {
   string fName = "Ks";
   TF1 *f(0);
   while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1_gauss2c", fName.c_str())))) if (f) delete f;

   f = new TF1(Form("%s_pol1_gauss2c", fName.c_str()), iF_pol1_gauss2c, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 7);
   f->SetParNames("norm", "peak", "sigma", "fraction", "sigma2", "constant", "slope");
   //  f->SetLineColor(kBlue);                                                                                                                                                                     
   f->SetLineWidth(2);

   int lbin(1), hbin(h->GetNbinsX()+1);
   double fLo = 99., fHi= -99.;
   if (fLo < fHi) {
     lbin = h->FindBin(fLo);
     hbin = h->FindBin(fHi);
   }
   double xlo(h->GetBinLowEdge(lbin));
   double xhi(h->GetBinLowEdge(hbin));

   double p0, p1;
   initPol1(p0, p1, h);
   if (fVerbose) cout << "p0: " << p0 << " p1: " << p1 << endl;
   double A   = 0.5*p1*(xhi*xhi - xlo*xlo) + p0*(xhi - xlo);

   double sqrt2pi = 2.506628275;
   double gInt    = h->Integral(lbin, hbin) - A;

   double gaussN  = gInt/(2.*sqrt2pi*sigma)*h->GetBinWidth(1);

   if (fVerbose) cout << "initFunc> gaussN = " << gaussN << " peak = " << peak << " sigma = " << sigma << " p0 = " << p0 << " p1 = " << p1 << endl;
   f->SetParameters(gaussN, peak, sigma, 0.2, 1.8*sigma, p0, p1);
   f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
   f->ReleaseParameter(1);     f->SetParLimits(1, peak-0.1, peak+0.1);
   f->ReleaseParameter(2);     f->SetParLimits(2, sigma*0.4, sigma*1.3);
   f->ReleaseParameter(3);     f->SetParLimits(3, 0.01, 2.0);
   f->ReleaseParameter(4);     f->SetParLimits(4, sigma*1.2, sigma*10.0);
   //applyLimits(f, "pol1gauss2c");
   return f;


 }
 double iF_argus_gauss2(double *x, double *par) {
   // par[0] -> Gaussian const
   // par[1] -> Gaussian mean
   // par[2] -> Gaussian sigma
   // par[3] -> fraction in second gaussian
   // par[4] -> sigma of second gaussian
   // par[5] -> Argus normalization
   // par[6] -> Argus exponential factor
   // par[7] -> Argus endpoint (fixed!)
   //           > 0: normal situation (zero above endpoint)
   //           < 0: inverted situation (zero below endpoint)

   //  double ebeam = 10.58/2;
   double ebeam = par[7];
   double ebeam2 = ebeam*ebeam;
   double background = 0.;
   double x2 = x[0]*x[0];
   double ratio = x2/ebeam2;
   if (par[7] < 0) ratio = ebeam2/x2;
   if (ratio < 1.) {
     background = par[5]*x[0] * TMath::Sqrt(1. - ratio) * TMath::Exp(par[6] * (1. - ratio));
   } else {
     background = 0.;
   }


   Double_t arg1(0.), arg2(0.), sig1(0.), sig2(0.);
   if (par[2] > 0) {
     arg1 = (x[0] - par[1]) / par[2];
     sig1 = par[0]*TMath::Exp(-0.5*arg1*arg1);
   }
   if (par[5] > 0.) {
     arg2 = (x[0] - par[1]) / par[4];
     sig2 = par[0]*TMath::Exp(-0.5*arg2*arg2);
   }

   return ((1-par[3])*sig1 + par[3]*sig2 + background);
 }

 TF1* plotFakenew_lxy::phiKK(TH1 *h) {
   int npar(8);
   string fName= "test";
   TF1* f = new TF1(Form("%s_phiKK", fName.c_str()), iF_argus_gauss2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
   f->SetParNames("const", "peak", "sigma1", "fraction", "sigma2", "norm.", "expo.", "endpoint");
   f->SetLineWidth(2);

   f->SetParameter(0, 2010.23);
   f->SetParameter(1, 1.019);
   f->SetParameter(2, 0.003);
   f->SetParameter(3, 0.4);
   f->SetParameter(4, 0.010);
   f->SetParameter(5, 30.);
   f->SetParameter(6, 3.);
   
   f->SetParameter(7, -0.89);
   //f->FixParameter(7, -2.0*0.498);
   //applyLimits(f, "phiKK");
   return f;
 }

void plotFakenew_lxy::sbsDistributions(TChain *tC,  string sample="bla", string what, string path) {
   cout << "plotFake::sbsDistributions(" << sample << ", " << what << ")" << endl;
   string sbsControlPlotsFileName = what.c_str();

   int NBINS = 50;  
   double fMassLo(0.),fMassHi(0.);
   if (string::npos != sample.find("ks")) {
     fMassLo    = 0.450;
     fMassHi    = 0.550;
     NBINS = 100;
   } 
   else if (string::npos != sample.find("d0")) {
     fMassLo    = 1.8;
     fMassHi    = 1.90;
   } else if (string::npos != sample.find("lambda")) {
     fMassLo    = 1.095;
     NBINS = 45;
     fMassHi    = 1.140;
   } 
   else if (string::npos != sample.find("phi")) {    
     fMassLo    = 1.0;
     fMassHi    = 1.04;
     NBINS = 40;
   } 
   string name="Ks_sbs_";


   int NREG_pt=8;
   int NREG_eta=9;
   int NREG_lxy=9;
   TString massbin_pt[NREG_pt];
   massbin_pt[0] = "bin0";
   massbin_pt[1] = "bin1";
   massbin_pt[2] = "bin2";
   massbin_pt[3] = "bin3";
   massbin_pt[4] = "bin4";
   massbin_pt[5] = "bin5";
   massbin_pt[6] = "bin6";   
   massbin_pt[7] = "allbin";

   TString massbin_lxy[NREG_lxy];
   massbin_lxy[0] = "bin0";
   massbin_lxy[1] = "bin1";
   massbin_lxy[2] = "bin2";
   massbin_lxy[3] = "bin3";
   massbin_lxy[4] = "bin4";
   massbin_lxy[5] = "bin5";
   massbin_lxy[6] = "bin6";   
   massbin_lxy[7] = "bin7";
   massbin_lxy[8] = "allbin";

   TString massbineta[NREG_eta];
   massbineta[0] = "bin0";
   massbineta[1] = "bin1";
   massbineta[2] = "bin2";
   massbineta[3] = "bin3";
   massbineta[4] = "bin4";
   massbineta[5] = "bin5";
   massbineta[6] = "bin6";
   massbineta[7] = "bin7";
   massbineta[8] = "allbin";

   TH1D* hv0_kin_eta_bin[NREG_eta];
   TH1D* hv0_kin_pt_bin[NREG_pt];
   TH1D* hv0_Mass_bin[NREG_pt];
   TH1D* hv0_Mass_bin_eta[NREG_eta];

   TH1D* hv0_kin_pt_bin_muid[NREG_pt];
   TH1D* hv0_Mass_bin_muid[NREG_pt];
   TH1D* hv0_Mass_bin_eta_muid[NREG_eta];
   TH1D* hv0_kin_eta_bin_muid[NREG_eta];

   TH1D* hv0_kin_bin_lxy_muid[NREG_lxy];
   TH1D* hv0_Mass_bin_lxy_muid[NREG_lxy];

   TH1D* hv0_kin_bin_lxy_softid[NREG_lxy];
   TH1D* hv0_Mass_bin_lxy_softid[NREG_lxy];

   TH1D* hv0_kin_bin_lxy_medid[NREG_lxy];
   TH1D* hv0_Mass_bin_lxy_medid[NREG_lxy];


   TH1D* hv0_kin_bin_lxy_mva20[NREG_lxy];
   TH1D* hv0_Mass_bin_lxy_mva20[NREG_lxy];
   TH1D* hv0_kin_bin_lxy_mva30[NREG_lxy];
   TH1D* hv0_Mass_bin_lxy_mva30[NREG_lxy];
   TH1D* hv0_kin_bin_lxy_mva40[NREG_lxy];
   TH1D* hv0_Mass_bin_lxy_mva40[NREG_lxy];
   TH1D* hv0_kin_bin_lxy_mva50[NREG_lxy];
   TH1D* hv0_Mass_bin_lxy_mva50[NREG_lxy];

   TH1D* hv0_Mass_bin_lxy[NREG_lxy];
   TH1D* hv0_kin_bin_lxy[NREG_lxy];

   for(int k=0;k<NREG_pt;k++){
     double var[]={0.,4.0,8.0,12.0,16.0,20.0, 30.0,50.};
     int ptbins=150;    
     hv0_kin_pt_bin[k] = new TH1D(Form("hv0_kin_pt_%s",massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins,0,50. );
     hv0_Mass_bin[k]   = new TH1D(Form("hv0_Mass_%s",massbin_pt[k].Data()), Form("%s_%s_v0_Mass", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo,fMassHi);

     hv0_Mass_bin_muid[k]   = new TH1D(Form("hv0_Mass_muid_%s",massbin_pt[k].Data()), Form("%s_%s_v0_Mass_muid", name.c_str(), massbin_pt[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_pt_bin_muid[k] = new TH1D(Form("hv0_kin_pt_muid_%s",massbin_pt[k].Data()), Form("%s_%s_v0_kin_pt_muid;p_{T}(GeV)", name.c_str(), massbin_pt[k].Data()), ptbins,0,50.);

     hv0_Mass_bin[k]->SetMinimum(0.);
   }

   for(int k=0;k<NREG_lxy;k++){
     int lxybins=150;    
     hv0_kin_bin_lxy[k] = new TH1D(Form("hv0_kin_lxy_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50. );
     hv0_Mass_bin_lxy[k]   = new TH1D(Form("hv0_Mass_lxy_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);

     hv0_Mass_bin_lxy_muid[k]   = new TH1D(Form("hv0_Mass_lxy_muid_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass_muid", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_bin_lxy_muid[k] = new TH1D(Form("hv0_kin_lxy_muid_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy_muid;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50.);
     
     hv0_Mass_bin_lxy_softid[k]   = new TH1D(Form("hv0_Mass_lxy_softid_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass_softid", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_bin_lxy_softid[k] = new TH1D(Form("hv0_kin_lxy_softid_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy_softid;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50.);
     hv0_Mass_bin_lxy_medid[k]   = new TH1D(Form("hv0_Mass_lxy_medid_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass_medid", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_bin_lxy_medid[k] = new TH1D(Form("hv0_kin_lxy_medid_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy_medid;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50.);
     
     
     hv0_Mass_bin_lxy_mva20[k]   = new TH1D(Form("hv0_Mass_lxy_mva20_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass_mva20", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_bin_lxy_mva20[k] = new TH1D(Form("hv0_kin_lxy_mva20_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy_mva20;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50.);

     hv0_Mass_bin_lxy_mva30[k]   = new TH1D(Form("hv0_Mass_lxy_mva30_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass_mva30", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_bin_lxy_mva30[k] = new TH1D(Form("hv0_kin_lxy_mva30_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy_mva30;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50.);

     hv0_Mass_bin_lxy_mva40[k]   = new TH1D(Form("hv0_Mass_lxy_mva40_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass_mva40", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_bin_lxy_mva40[k] = new TH1D(Form("hv0_kin_lxy_mva40_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy_mva40;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50.);

     hv0_Mass_bin_lxy_mva50[k]   = new TH1D(Form("hv0_Mass_lxy_mva50_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_Mass_mva50", name.c_str(), massbin_lxy[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_bin_lxy_mva50[k] = new TH1D(Form("hv0_kin_lxy_mva50_%s",massbin_lxy[k].Data()), Form("%s_%s_v0_kin_lxy_mva50;p_{T}(GeV)", name.c_str(), massbin_lxy[k].Data()), lxybins,0,50.);
     
     
     //hv0_Mass_bin[k]->SetMinimum(0.);
   }

   for(int k=0;k<NREG_eta;k++){

     hv0_kin_eta_bin[k] = new TH1D(Form("hv0_kin_eta_%s",massbineta[k].Data()), Form("%s_%s_v0_kin_eta", name.c_str(), massbineta[k].Data()), 100, -2.,2. );
     hv0_Mass_bin_eta[k]   = new TH1D(Form("hv0_Mass_eta_%s",massbineta[k].Data()), Form("%s_%s_v0_Mass", name.c_str(), massbineta[k].Data()), NBINS, fMassLo,fMassHi);

     hv0_Mass_bin_eta_muid[k]   = new TH1D(Form("hv0_Mass_eta_muid_%s",massbineta[k].Data()), Form("%s_%s_v0_Mass_muid", name.c_str(), massbineta[k].Data()), NBINS, fMassLo,fMassHi);
     hv0_kin_eta_bin_muid[k] = new TH1D(Form("hv0_kin_eta_muid_%s",massbineta[k].Data()), Form("%s_%s_v0_kin_eta_muid", name.c_str(), massbineta[k].Data()), 100, -2.,2. );

   }


   UInt_t          run;
   UInt_t          luminosityBlock;
   UInt_t          nMuon;
   Bool_t          HLT_Ele30;
   Float_t         Muon_pt[150];
   Bool_t          Muon_softMvaId[100];
   Float_t         Muon_softMva[100];
   Bool_t          Muon_softId[100];
   Bool_t          Muon_mediumId[150];
   Bool_t          Muon_isGlobal[150];
   Bool_t          Muon_isTracker[150];
   Bool_t          Muon_looseId[150];
   Float_t         MuonId_hlt_pt[150];

   UInt_t          nV0;
   Float_t         V0_doca[100];                                                                                                                                  
   Float_t         V0_kin_cosAlphaXY[100];                                                                                                                        
   Float_t         V0_kin_eta[100];                                                                                                                               
   Float_t         V0_kin_lxy[100];                                                                                                                               
   Float_t         V0_kin_mass[100];                                                                                                                              
   Float_t         V0_kin_massErr[100];                                                                                                                           
   Float_t         V0_kin_phi[100];                                                                                                                               
   Float_t         V0_kin_pt[100];                                                                                                                                
   Float_t         V0_kin_slxy[100];                                                                                                                              
   Float_t         V0_kin_vtx_chi2dof[100];                                                                                                                       
   Float_t         V0_kin_vtx_prob[100];                                                                                                                          
   Float_t         V0_mass[100];                                                                                                                                  
   Float_t         V0_trk1_eta[100];                                                                                                                              
   Float_t         V0_trk1_phi[100];                                                                                                                              
   Float_t         V0_trk1_pt[100];                                                                                                                               
   Float_t         V0_trk2_eta[100];                                                                                                                              
   Float_t         V0_trk2_phi[100];                                                                                                                              
   Float_t         V0_trk2_pt[100];                                                                                                                               
   Int_t           V0_kin_valid[100];                                                                                                                             
   Int_t           V0_trk1_mu_index[100];                                                                                                                         
   Int_t           V0_trk2_mu_index[100];   

   Int_t           V0_gen_pdgId[100];                         
   Int_t           V0_gen_trk1_mpdgId[100];                     
   Int_t           V0_gen_trk1_pdgId[100];                      
   Int_t           V0_gen_trk2_mpdgId[100];                     
   Int_t           V0_gen_trk2_pdgId[100];  
   Float_t         V0_kin_sipPV[100];
   Float_t         V0_trk1_sip[100];
   Float_t         V0_trk2_sip[100];

   TBranch        *b_run;
   TBranch        *b_luminosityBlock;
   TBranch        *b_nMuon;
   TBranch        *b_HLT_Ele30;
   TBranch        *b_Muon_pt;
   TBranch        *b_Muon_softMvaId;
   TBranch        *b_Muon_softMva;
   TBranch        *b_Muon_softId;
   TBranch        *b_Muon_mediumId;
   TBranch        *b_Muon_isGlobal;
   TBranch        *b_Muon_isTracker;
   TBranch        *b_Muon_looseId;
   TBranch        *b_MuonId_hlt_pt;

   TBranch        *b_nV0;                                               
   TBranch        *b_V0_doca;                                           
   TBranch        *b_V0_kin_cosAlphaXY;                                 
   TBranch        *b_V0_kin_eta;                                        
   TBranch        *b_V0_kin_lxy;                                        
   TBranch        *b_V0_kin_mass;                                       
   TBranch        *b_V0_kin_massErr;                                    
   TBranch        *b_V0_kin_phi;                                        
   TBranch        *b_V0_kin_pt;                                         
   TBranch        *b_V0_kin_slxy;                                       
   TBranch        *b_V0_kin_vtx_chi2dof;                                
   TBranch        *b_V0_kin_vtx_prob;                                   
   TBranch        *b_V0_mass;                                           
   TBranch        *b_V0_trk1_eta;                                       
   TBranch        *b_V0_trk1_phi;                                       
   TBranch        *b_V0_trk1_pt;                                        
   TBranch        *b_V0_trk2_eta;                                       
   TBranch        *b_V0_trk2_phi;                                       
   TBranch        *b_V0_trk2_pt;                                        
   TBranch        *b_V0_kin_valid;                                      
   TBranch        *b_V0_trk1_mu_index;                                  
   TBranch        *b_V0_trk2_mu_index;
   TBranch        *b_V0_kin_sipPV;
   TBranch        *b_V0_trk1_sip;
   TBranch        *b_V0_trk2_sip;


   TBranch        *b_V0_gen_pdgId;
   TBranch        *b_V0_gen_trk1_mpdgId;
   TBranch        *b_V0_gen_trk1_pdgId;
   TBranch        *b_V0_gen_trk2_mpdgId;
   TBranch        *b_V0_gen_trk2_pdgId;

   tC->SetBranchAddress("run", &run, &b_run);
   tC->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   tC->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   if(string::npos != sample.find("trigger")) {
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

   if(string::npos != sample.find("DoubleMuMu")) {
      tC->SetBranchAddress("MuonId_hlt_pt", MuonId_hlt_pt, &b_MuonId_hlt_pt);
   }

   if(string::npos != sample.find("ks")){
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

     if(string::npos != what.find("MC")){
       tC->SetBranchAddress("ks_gen_pdgId", V0_gen_pdgId, &b_V0_gen_pdgId);
       tC->SetBranchAddress("ks_gen_trk1_mpdgId",V0_gen_trk1_mpdgId, &b_V0_gen_trk1_mpdgId);
       tC->SetBranchAddress("ks_gen_trk2_mpdgId",V0_gen_trk2_mpdgId, &b_V0_gen_trk2_mpdgId);
       tC->SetBranchAddress("ks_gen_trk1_pdgId",V0_gen_trk1_pdgId, &b_V0_gen_trk1_pdgId);
       tC->SetBranchAddress("ks_gen_trk2_pdgId",V0_gen_trk2_pdgId, &b_V0_gen_trk2_pdgId);
     }
   }
   else if(string::npos != sample.find("lambda")){
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
   else if(string::npos != sample.find("phi")){
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
   else if(string::npos != sample.find("d0")){
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
   Long64_t n = tC->GetEntries();
   int count0=0; int count1=0, count2=0, count3(0);
   TString cut_print="bla";
   TString cut_print_mid="bla";

   int i_permille_old = 0;
   for (unsigned int i=0; i < n; ++i){
     int i_permille = (int)floor(100. * i / n);
     if (i_permille != i_permille_old) {
       printf("\015\033[32m ---> \033[1m\033[31m%d%%"
	      "\033[0m\033[32m <---\033[0m\015", i_permille);
       fflush(stdout);
       i_permille_old = i_permille;
     }

     // Load a proper try and get relative even index
     Long64_t localEntry = tC->LoadTree(i);
     b_nV0->GetEntry(localEntry);
     if (nV0 == 0) continue;
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
     if(string::npos != sample.find("trigger")) {
         b_HLT_Ele30->GetEntry(localEntry);
     }

     if(string::npos != sample.find("DoubleMuMu")) {
         b_MuonId_hlt_pt->GetEntry(localEntry);
     }
     if(string::npos != what.find("MC")){
       b_V0_gen_pdgId->GetEntry(localEntry);    
       b_V0_gen_trk1_mpdgId->GetEntry(localEntry);
       b_V0_gen_trk2_mpdgId->GetEntry(localEntry);
       b_V0_gen_trk1_pdgId->GetEntry(localEntry);
       b_V0_gen_trk2_pdgId->GetEntry(localEntry);
     }
     for (unsigned int bhad=0; bhad < nV0; ++bhad){
       count1++;

       // number of pt bins
       bool Finalcut_bin0 = V0_kin_pt[bhad] > 0.0 &&  V0_kin_pt[bhad] < 4.0;
       bool Finalcut_bin1 = V0_kin_pt[bhad] > 4.0 && V0_kin_pt[bhad] < 8.0;
       bool Finalcut_bin2 = V0_kin_pt[bhad] > 8.0 && V0_kin_pt[bhad] < 12.0;
       bool Finalcut_bin3 = V0_kin_pt[bhad] > 12.0 && V0_kin_pt[bhad] < 16.0;
       bool Finalcut_bin4 = V0_kin_pt[bhad] > 16.0 && V0_kin_pt[bhad] < 20.0;
       bool Finalcut_bin5 = V0_kin_pt[bhad] > 20.0 && V0_kin_pt[bhad] < 30.0;
       bool Finalcut_bin6 = V0_kin_pt[bhad] > 30.0 && V0_kin_pt[bhad] < 50.0;
       bool Finalcut_bin7 = V0_kin_pt[bhad] > 0.0 && V0_kin_pt[bhad] < 50.0;
       bool cut_pre = false;
       bool cut_pre_mid = false;
       bool cut_pre_medid = false;
       bool cut_pre_softid = false;
       bool cut_pre_mid20 = false;
       bool cut_pre_mid30 = false;
       bool cut_pre_mid40 = false;
       bool cut_pre_mid50 = false;
       

     int lxy_cut = 40;
     if(string::npos != sample.find("ks") && string::npos != sample.find("DoubleMuMu")){
           cut_pre = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==2);

           cut_pre_softid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==3 && (((V0_trk1_mu_index[bhad]==0 ^ V0_trk2_mu_index[bhad]==0) && MuonId_hlt_pt[1]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==1 ^ V0_trk2_mu_index[bhad]==1) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==2 ^ V0_trk2_mu_index[bhad]==2) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[1]>0)) && ((Muon_isGlobal[V0_trk1_mu_index[bhad]] && Muon_isTracker[V0_trk1_mu_index[bhad]] && Muon_looseId[V0_trk1_mu_index[bhad]]) || (Muon_isGlobal[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]] && Muon_looseId[V0_trk2_mu_index[bhad]])));

            cut_pre_medid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==3 && (((V0_trk1_mu_index[bhad]==0 ^ V0_trk2_mu_index[bhad]==0) && MuonId_hlt_pt[1]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==1 ^ V0_trk2_mu_index[bhad]==1) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==2 ^ V0_trk2_mu_index[bhad]==2) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[1]>0)) && (Muon_mediumId[V0_trk1_mu_index[bhad]] || Muon_mediumId[V0_trk2_mu_index[bhad]]));

            cut_pre_mid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==3 && (((V0_trk1_mu_index[bhad]==0 ^ V0_trk2_mu_index[bhad]==0) && MuonId_hlt_pt[1]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==1 ^ V0_trk2_mu_index[bhad]==1) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==2 ^ V0_trk2_mu_index[bhad]==2) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[1]>0)) && (Muon_softMvaId[V0_trk1_mu_index[bhad]] || Muon_softMvaId[V0_trk2_mu_index[bhad]]));


            cut_pre_mid20 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==3 && (((V0_trk1_mu_index[bhad]==0 ^ V0_trk2_mu_index[bhad]==0) && MuonId_hlt_pt[1]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==1 ^ V0_trk2_mu_index[bhad]==1) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==2 ^ V0_trk2_mu_index[bhad]==2) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[1]>0)) && (Muon_softMva[V0_trk1_mu_index[bhad]]>0.20 || Muon_softMva[V0_trk2_mu_index[bhad]]>0.20));

            cut_pre_mid30 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==3 && (((V0_trk1_mu_index[bhad]==0 ^ V0_trk2_mu_index[bhad]==0) && MuonId_hlt_pt[1]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==1 ^ V0_trk2_mu_index[bhad]==1) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==2 ^ V0_trk2_mu_index[bhad]==2) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[1]>0)) && (Muon_softMva[V0_trk1_mu_index[bhad]]>0.30 || Muon_softMva[V0_trk2_mu_index[bhad]]>0.30));

            cut_pre_mid40 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==3 && (((V0_trk1_mu_index[bhad]==0 ^ V0_trk2_mu_index[bhad]==0) && MuonId_hlt_pt[1]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==1 ^ V0_trk2_mu_index[bhad]==1) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==2 ^ V0_trk2_mu_index[bhad]==2) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[1]>0)) && (Muon_softMva[V0_trk1_mu_index[bhad]]>0.40 || Muon_softMva[V0_trk2_mu_index[bhad]]>0.40));


            cut_pre_mid50 = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && nMuon==3 && (((V0_trk1_mu_index[bhad]==0 ^ V0_trk2_mu_index[bhad]==0) && MuonId_hlt_pt[1]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==1 ^ V0_trk2_mu_index[bhad]==1) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[2]>0) || ((V0_trk1_mu_index[bhad]==2 ^ V0_trk2_mu_index[bhad]==2) && MuonId_hlt_pt[0]>0 && MuonId_hlt_pt[1]>0)) && (Muon_softMva[V0_trk1_mu_index[bhad]]>0.50 || Muon_softMva[V0_trk2_mu_index[bhad]]>0.50));


       }
/*
     if(string::npos != sample.find("ks") && string::npos != what.find("MC")){

	 //cout << "trk1: " << V0_gen_trk1_pdgId[bhad] <<", mother: " << V0_gen_trk1_mpdgId[bhad] << endl;
         //cout << "trk2: " << V0_gen_trk2_pdgId[bhad] <<", mother: " << V0_gen_trk2_mpdgId[bhad] << endl;

         cut_pre =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && (((V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) && abs(V0_gen_trk1_mpdgId[bhad])==311) || ((V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) && abs(V0_gen_trk2_mpdgId[bhad])==311)));

         cut_pre_softid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_isGlobal[V0_trk1_mu_index[bhad]] && Muon_isTracker[V0_trk1_mu_index[bhad]] && Muon_looseId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_looseId[V0_trk2_mu_index[bhad]] && Muon_isGlobal[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

         cut_pre_medid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

         cut_pre_mid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

         cut_pre_mid20 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.20) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.20)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
         cut_pre_mid30 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.30) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.30)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
         cut_pre_mid40 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.40) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.40)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));
         cut_pre_mid50 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.50) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.50)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3 && ((( V0_gen_trk1_pdgId[bhad]==310 || V0_gen_trk1_pdgId[bhad]==-310) &&  abs(V0_gen_trk1_mpdgId[bhad])==311) || (( V0_gen_trk2_pdgId[bhad]==310 || V0_gen_trk2_pdgId[bhad]==-310) &&  abs(V0_gen_trk2_mpdgId[bhad])==311)));

       }
       else if(string::npos != sample.find("ks") && string::npos == what.find("MC")){
*/
       else if(string::npos != sample.find("ks")) {

         cut_pre =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);

         cut_pre_softid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_isGlobal[V0_trk1_mu_index[bhad]] && Muon_isTracker[V0_trk1_mu_index[bhad]] && Muon_looseId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_looseId[V0_trk2_mu_index[bhad]] && Muon_isGlobal[V0_trk2_mu_index[bhad]] && Muon_isTracker[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);

         cut_pre_medid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_mediumId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);

         cut_pre_mid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]])) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);

         cut_pre_mid20 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.20) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.20)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);
         cut_pre_mid30 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.30) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.30)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);
         cut_pre_mid40 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.40) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.40)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);
         cut_pre_mid50 =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMva[V0_trk1_mu_index[bhad]]>0.50) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMva[V0_trk2_mu_index[bhad]]>0.50)) && abs(V0_kin_lxy[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad]>0.999 && V0_doca[bhad]<0.004 &&  V0_kin_sipPV[bhad]<3 && V0_trk1_sip[bhad]>5 && V0_trk2_sip[bhad]>5 && abs(V0_kin_slxy[bhad])>3);

       }
       else if(string::npos != sample.find("phi")){
	 cut_pre =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 3.0) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 3.0)) &&  V0_kin_vtx_prob[bhad]>0.3 && V0_kin_sipPV[bhad]<1  && V0_doca[bhad]<0.004 && abs(V0_kin_lxy[bhad])>0);
	 cut_pre_mid =(((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk1_mu_index[bhad]]) || ( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]])) &&  V0_kin_vtx_prob[bhad]>0.3 && V0_kin_sipPV[bhad]<1&& V0_doca[bhad]<0.004 && abs(V0_kin_lxy[bhad])>0);
		       
   
       }
       else if(string::npos != sample.find("lambda")){
	 cut_pre = V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_kin_slxy[bhad]>3 &&  V0_kin_sipPV[bhad]<1 && V0_trk1_sip[bhad]>2 && V0_trk2_sip[bhad]>2;
	 cut_pre_mid =( V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]] && V0_kin_slxy[bhad]>3 &&  V0_kin_sipPV[bhad]<1 && V0_trk1_sip[bhad]>2 && V0_trk2_sip[bhad]>2);
	 
       }

       if(cut_pre)count2++;
       for(int k=0; k<NREG_pt; k++){
	 
	 bool cut_tmp = false;
	 if(k == 0) cut_tmp = Finalcut_bin0 ;
	 if(k == 1) cut_tmp = Finalcut_bin1 ;
	 if(k == 2) cut_tmp = Finalcut_bin2 ;
	 if(k == 3) cut_tmp = Finalcut_bin3 ;
	 if(k == 4) cut_tmp = Finalcut_bin4 ;
	 if(k == 5) cut_tmp = Finalcut_bin5 ;
	 if(k == 6) cut_tmp = Finalcut_bin6 ;
	 if(k == 7) cut_tmp = Finalcut_bin7 ;
	 if(cut_tmp){	                       
	   if(cut_pre){	     	   
	     hv0_kin_pt_bin[k]->Fill(V0_kin_pt[bhad]);
	     hv0_Mass_bin[k]->Fill(V0_kin_mass[bhad]);	    	   	  	       
	   }
	   if( cut_pre_mid){
	       
	       hv0_Mass_bin_muid[k]->Fill(V0_kin_mass[bhad]);
	       hv0_kin_pt_bin_muid[k]->Fill(V0_kin_pt[bhad]);	     
	   }	  
	   
	 }
       }// loop for pt bins
       // number of Lxy bins
       bool Finalcut_lxy_bin0 = abs(V0_kin_lxy[bhad]) > 0.0 &&  abs(V0_kin_lxy[bhad]) < 1.;
       bool Finalcut_lxy_bin1 = abs(V0_kin_lxy[bhad]) > 1.0 && abs(V0_kin_lxy[bhad]) < 2.;
       bool Finalcut_lxy_bin2 = abs(V0_kin_lxy[bhad]) > 2.0 && abs(V0_kin_lxy[bhad]) < 4.0;
       bool Finalcut_lxy_bin3 = abs(V0_kin_lxy[bhad]) > 4.0 && abs(V0_kin_lxy[bhad]) < 8.0;
       
       bool Finalcut_lxy_bin4 = abs(V0_kin_lxy[bhad]) > 8.0 && abs(V0_kin_lxy[bhad]) < 12.;
       bool Finalcut_lxy_bin5 = abs(V0_kin_lxy[bhad]) > 12. && abs(V0_kin_lxy[bhad]) < 16;
       bool Finalcut_lxy_bin6 = abs(V0_kin_lxy[bhad]) > 16.0 && abs(V0_kin_lxy[bhad]) < 30.0;
       bool Finalcut_lxy_bin7 = abs(V0_kin_lxy[bhad]) > 30.0 && abs(V0_kin_lxy[bhad]) < 50;
       bool Finalcut_lxy_bin8 =abs(V0_kin_lxy[bhad]) > 0.0 && abs(V0_kin_lxy[bhad]) < 50.0;
       
       double w=1.0;
       //if(V0_trk1_pt[bhad]>4.0 && V0_trk2_pt[bhad]>4.0 && abs(V0_trk1_eta[bhad])<1.4 && abs(V0_trk2_eta[bhad])<1.4) w=2.0;
      
       for(int k=0; k<NREG_lxy; k++){
	 bool cut_tmp = false;
	 if(k == 0) cut_tmp = Finalcut_lxy_bin0 ;
	 if(k == 1) cut_tmp = Finalcut_lxy_bin1 ;
	 if(k == 2) cut_tmp = Finalcut_lxy_bin2 ;
	 if(k == 3) cut_tmp = Finalcut_lxy_bin3 ;
	 if(k == 4) cut_tmp = Finalcut_lxy_bin4 ;
	 if(k == 5) cut_tmp = Finalcut_lxy_bin5 ;
	 if(k == 6) cut_tmp = Finalcut_lxy_bin6 ;
	 if(k == 7) cut_tmp = Finalcut_lxy_bin7 ;
	 if(k == 8) cut_tmp = Finalcut_lxy_bin8 ;
	 
	   if(cut_tmp){	    
	     if(cut_pre){
	       hv0_kin_bin_lxy[k]->Fill(V0_kin_lxy[bhad],w);
	       hv0_Mass_bin_lxy[k]->Fill(V0_kin_mass[bhad],w);
	     }
	     if(cut_pre_medid){
               hv0_Mass_bin_lxy_medid[k]->Fill(V0_kin_mass[bhad]);
               hv0_kin_bin_lxy_medid[k]->Fill(V0_kin_lxy[bhad]);
             }
	     if(cut_pre_softid){
               hv0_Mass_bin_lxy_softid[k]->Fill(V0_kin_mass[bhad]);
               hv0_kin_bin_lxy_softid[k]->Fill(V0_kin_lxy[bhad]);
             }
	     if(cut_pre_mid){		
	       hv0_Mass_bin_lxy_muid[k]->Fill(V0_kin_mass[bhad]);
	       hv0_kin_bin_lxy_muid[k]->Fill(V0_kin_lxy[bhad]);	      
	     }
	     if(cut_pre_mid20){		
	       hv0_Mass_bin_lxy_mva20[k]->Fill(V0_kin_mass[bhad]);
	       hv0_kin_bin_lxy_mva20[k]->Fill(V0_kin_lxy[bhad]);	      
	     }
	     if(cut_pre_mid30){		
	       hv0_Mass_bin_lxy_mva30[k]->Fill(V0_kin_mass[bhad]);
	       hv0_kin_bin_lxy_mva30[k]->Fill(V0_kin_lxy[bhad]);	      
	     }
	     if(cut_pre_mid40){		
	       hv0_Mass_bin_lxy_mva40[k]->Fill(V0_kin_mass[bhad]);
	       hv0_kin_bin_lxy_mva40[k]->Fill(V0_kin_lxy[bhad]);	      
	     }
	     if(cut_pre_mid50){		
	       hv0_Mass_bin_lxy_mva50[k]->Fill(V0_kin_mass[bhad]);
	       hv0_kin_bin_lxy_mva50[k]->Fill(V0_kin_lxy[bhad]);	      
	     }
	   }
       }// loop for lxy bins	        
     }
   }


   string mode = "sideband";
   if((string::npos != what.find("MC"))) mode = "MC";

   cout<<" count "<<count0<<"\t"<<count1<<"\t"<<count2<<endl;
   
   cout<< " lxy "<<endl;
   cout<< hv0_kin_bin_lxy[0]->Integral()<<"\t"<< hv0_kin_bin_lxy[1]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy[2]->Integral()<<"\t"<< hv0_kin_bin_lxy[3]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy[4]->Integral()<<"\t"<< hv0_kin_bin_lxy[5]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy[6]->Integral()<<"\t"<< hv0_kin_bin_lxy[7]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy[0]->Integral()<<"\t"<< hv0_kin_bin_lxy[8]->Integral()<<endl;
   cout<< "muid  lxy "<<endl;
   cout<< hv0_kin_bin_lxy_muid[0]->Integral()<<"\t"<< hv0_kin_bin_lxy_muid[1]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy_muid[2]->Integral()<<"\t"<< hv0_kin_bin_lxy_muid[3]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy_muid[4]->Integral()<<"\t"<< hv0_kin_bin_lxy_muid[5]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy_muid[6]->Integral()<<"\t"<< hv0_kin_bin_lxy_muid[7]->Integral()<<endl;
   cout<< hv0_kin_bin_lxy_muid[0]->Integral()<<"\t"<< hv0_kin_bin_lxy_muid[8]->Integral()<<endl;
   cout<< " Mass"<<endl;
   cout<< hv0_Mass_bin_lxy[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy[1]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy[2]->Integral()<<"\t"<< hv0_Mass_bin_lxy[3]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy[4]->Integral()<<"\t"<< hv0_Mass_bin_lxy[5]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy[6]->Integral()<<"\t"<< hv0_Mass_bin_lxy[7]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy[8]->Integral()<<endl;
   cout<< "muid Mass "<<endl;
   cout<< hv0_Mass_bin_lxy_muid[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy_muid[1]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_muid[2]->Integral()<<"\t"<< hv0_Mass_bin_lxy_muid[3]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_muid[4]->Integral()<<"\t"<< hv0_Mass_bin_lxy_muid[5]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_muid[6]->Integral()<<"\t"<< hv0_Mass_bin_lxy_muid[7]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_muid[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy_muid[8]->Integral()<<endl;
   cout<< "mva 30 Mass "<<endl;
   cout<< hv0_Mass_bin_lxy_mva30[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva30[1]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva30[2]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva30[3]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva30[4]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva30[5]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva30[6]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva30[7]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva30[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva30[8]->Integral()<<endl;
   cout<< "mva 50 Mass "<<endl;
   cout<< hv0_Mass_bin_lxy_mva50[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva50[1]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva50[2]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva50[3]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva50[4]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva50[5]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva50[6]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva50[7]->Integral()<<endl;
   cout<< hv0_Mass_bin_lxy_mva50[0]->Integral()<<"\t"<< hv0_Mass_bin_lxy_mva50[8]->Integral()<<endl;
 
 
   cout<< " pt "<<endl;
   cout<< hv0_kin_pt_bin[0]->Integral()<<"\t"<< hv0_kin_pt_bin[1]->Integral()<<endl;
   cout<< hv0_kin_pt_bin[2]->Integral()<<"\t"<< hv0_kin_pt_bin[3]->Integral()<<endl;
   cout<< hv0_kin_pt_bin[4]->Integral()<<"\t"<< hv0_kin_pt_bin[5]->Integral()<<endl;

   cout<< "muid  pt "<<endl;
   cout<< hv0_kin_pt_bin_muid[0]->Integral()<<"\t"<< hv0_kin_pt_bin_muid[1]->Integral()<<endl;
   cout<< hv0_kin_pt_bin_muid[2]->Integral()<<"\t"<< hv0_kin_pt_bin_muid[3]->Integral()<<endl;
   cout<< hv0_kin_pt_bin_muid[4]->Integral()<<"\t"<< hv0_kin_pt_bin_muid[5]->Integral()<<endl;

   cout<< " mass pt "<<endl;
   cout<< hv0_Mass_bin[0]->Integral()<<"\t"<< hv0_Mass_bin[1]->Integral()<<endl;
   cout<< hv0_Mass_bin[2]->Integral()<<"\t"<< hv0_Mass_bin[3]->Integral()<<endl;
   cout<< hv0_Mass_bin[4]->Integral()<<"\t"<< hv0_Mass_bin[5]->Integral()<<endl;

   cout<< "muid  pt "<<endl;
   cout<< hv0_Mass_bin_muid[0]->Integral()<<"\t"<< hv0_Mass_bin_muid[1]->Integral()<<endl;
   cout<< hv0_Mass_bin_muid[2]->Integral()<<"\t"<< hv0_Mass_bin_muid[3]->Integral()<<endl;
   cout<< hv0_Mass_bin_muid[4]->Integral()<<"\t"<< hv0_Mass_bin_muid[5]->Integral()<<endl;

   
   TString new_filename= Form("%shisto_Weight_%s.root",path.c_str(), sample.c_str());
   TFile *output = TFile::Open(new_filename.Data(),"RECREATE");
   for(int k(0);k<NREG_pt;k++ ){  
     hv0_kin_pt_bin[k]->Write();
     hv0_Mass_bin[k]->Write();   
     hv0_kin_pt_bin_muid[k]->Write();
     hv0_Mass_bin_muid[k]->Write();
   }
   for(int ik(0.);ik<NREG_eta;ik++){
      hv0_kin_eta_bin[ik]->Write();
      hv0_Mass_bin_eta[ik]->Write();
      hv0_Mass_bin_eta_muid[ik]->Write();
      hv0_kin_eta_bin_muid[ik]->Write();
   }
   for(int ik(0.);ik<NREG_lxy;ik++){
     hv0_kin_bin_lxy[ik]->Write();
     hv0_Mass_bin_lxy[ik]->Write();
     hv0_Mass_bin_lxy_muid[ik]->Write();
     hv0_kin_bin_lxy_muid[ik]->Write();
     hv0_Mass_bin_lxy_softid[ik]->Write();
     hv0_kin_bin_lxy_softid[ik]->Write();
     hv0_Mass_bin_lxy_medid[ik]->Write();
     hv0_kin_bin_lxy_medid[ik]->Write();
     hv0_Mass_bin_lxy_mva20[ik]->Write();
     hv0_kin_bin_lxy_mva20[ik]->Write();
     hv0_Mass_bin_lxy_mva30[ik]->Write();
     hv0_kin_bin_lxy_mva30[ik]->Write();
     hv0_Mass_bin_lxy_mva40[ik]->Write();
     hv0_kin_bin_lxy_mva40[ik]->Write();
     hv0_Mass_bin_lxy_mva50[ik]->Write();
     hv0_kin_bin_lxy_mva50[ik]->Write();

   }

   output->Write();
   output->Close();
   for(int ik(0.);ik<NREG_lxy;ik++){
     hv0_kin_bin_lxy[ik]->Delete();
     hv0_Mass_bin_lxy[ik]->Delete();
     hv0_Mass_bin_lxy_muid[ik]->Delete();
     hv0_kin_bin_lxy_muid[ik]->Delete();
     hv0_Mass_bin_lxy_softid[ik]->Delete();
     hv0_kin_bin_lxy_softid[ik]->Delete();
     hv0_Mass_bin_lxy_medid[ik]->Delete();
     hv0_kin_bin_lxy_medid[ik]->Delete();
  
     hv0_Mass_bin_lxy_mva20[ik]->Delete();
     hv0_kin_bin_lxy_mva20[ik]->Delete();
     hv0_Mass_bin_lxy_mva30[ik]->Delete();
     hv0_kin_bin_lxy_mva30[ik]->Delete();
     hv0_Mass_bin_lxy_mva40[ik]->Delete();
     hv0_kin_bin_lxy_mva40[ik]->Delete();
     hv0_Mass_bin_lxy_mva50[ik]->Delete();
     hv0_kin_bin_lxy_mva50[ik]->Delete();
   }

   for(int k=0;k<NREG_pt;k++){
     hv0_kin_pt_bin[k]->Delete(); 
     hv0_Mass_bin[k]->Delete();

     hv0_Mass_bin_muid[k]->Delete();
     hv0_kin_pt_bin_muid[k]->Delete(); 

   }
   cout<<" partially "<<endl;
   for(int k=0;k<NREG_eta;k++){
     hv0_kin_eta_bin[k]->Delete(); 
     hv0_Mass_bin_eta[k]->Delete();

     hv0_Mass_bin_eta_muid[k]->Delete();
     hv0_kin_eta_bin_muid[k]->Delete(); 

   }
   cout<<"complete"<<endl;

   return;
 }
 // //.........

 // ----------------------------------------------------------------------
void plotFakenew_lxy::loopOverChain(TChain *tC, string basic="Data", string modf="noprod", string path="ks/") {
      
  //if(string::npos != modf.find("Data")) sbsDistributions_bin(tC, cutsks, "ks", modf);
  if(string::npos != modf.find("Histproduction")) sbsDistributions(tC, basic, modf, path);
  //playfit(basic, path);
  // playPhifit(basic, path);
  
  cout<<"lopping overchain function closing"<<endl;
}

void plotFakenew_lxy::list_dir(TChain *tC, const char *path) {
  struct dirent *entry;
  DIR *dir = opendir(path);
  if (dir == NULL) {
    cout<<"No file is in this drectory"<<endl;
    return;
  }
  while ((entry = readdir(dir)) != NULL) {    
    string filenamec = entry->d_name;
    string s1 = ".root";
    
    if (filenamec.find(s1) != string::npos) {
      TString filename_f = path+filenamec;
      //std::cout << "found! " << filename_f.c_str()<< '\n';
      TFile *file = TFile::Open(filename_f,"READ"); 
      if (!file) continue ; //file does not exist                                                                                            
      if (file->IsZombie()) continue; //file is unusable                                                                                                               
      if (file->TestBit(TFile::kRecovered)) continue; //file has been recovered                                                                                        
      delete file;                                                                                                                                       
      tC->Add(filename_f); 
    }
    
    
   }
  closedir(dir);
}
void plotFakenew_lxy::list_dir_file(TChain *tC, string filename) {
  fstream newfile;
  newfile.open(filename,ios::in); //open a file to perform read operation using file object                                                                             
  if (newfile.is_open()){   //checking whether the file is open                                                                                                         
    string tp;
    while(getline(newfile, tp)){ //read data from file object and put it into string.                                                                                    
      //cout << tp << "\n"; //print the data of the string                                                                                                              
      const char * c = tp.c_str();
      list_dir(tC, c);
    }
    newfile.close(); //close the file object.                                                                                                                            
  }
}



 int main()
 {
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   RooMsgService::instance().setSilentMode(kTRUE);
   RooMsgService::instance().setStreamStatus(1,false);
   RooMsgService::instance().getStream(1).removeTopic(Integration) ;
   RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
   RooMsgService::instance().getStream(1).removeTopic(Fitting) ;
   RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
   RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
   RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
   RooMsgService::instance().getStream(1).removeTopic(Eval) ;
   RooMsgService::instance().Print() ;

   
   plotFakenew_lxy c1;  
   TChain* tC_22= new TChain("Events");                                                                                                                            
   c1.list_dir_file(tC_22, "2022_MC.txt");                                                                                                                    
   cout<<"Entry Data "<<tC_22->GetEntries()<<endl;                                                                                                                   
   c1.loopOverChain(tC_22, "Data_521_2022_ks_MCall_lxy_220821","Histproduction", "ks_MC_new/");        
   delete tC_22;                                           

                                                                                                       
   /*plotFakenew_lxy c2;
   TChain* tC_22_MC= new TChain("Events");                                                                                                                    
   c2.list_dir_file(tC_22_MC, "2022_data_ks_parking_small.txt");                                                                                                            
   cout<<"Entry MC "<<tC_22_MC->GetEntries()<<endl;                                                                                                           
   c2.loopOverChain(tC_22_MC, "Data_521_2022_ks_parking_notrig_lxy_220821","Histproduction", "ks_triggers/");
   delete tC_22_MC; */
   /*
   plotFakenew_lxy c3;
   TChain* tC_22_DoubleMuMu1 = new TChain("Events");
   c3.list_dir_file(tC_22_DoubleMuMu1, "2022_data_ks_DoubleMuMu1.txt");
   cout << "Entry Data DoubleMuMu" << tC_22_DoubleMuMu1->GetEntries() << endl;
   c3.loopOverChain(tC_22_DoubleMuMu1, "Data_521_2022_ks_DoubleMuMu_slxy3_lxy8_v1", "Histproduction", "ks_DoubleMuMu/");
   delete tC_22_DoubleMuMu1;

   plotFakenew_lxy c4;
   TChain* tC_22_DoubleMuMu2 = new TChain("Events");
   c4.list_dir_file(tC_22_DoubleMuMu2, "2022_data_ks_DoubleMuMu2.txt");
   cout << "Entry Data DoubleMuMu" << tC_22_DoubleMuMu2->GetEntries() << endl;
   c4.loopOverChain(tC_22_DoubleMuMu2, "Data_521_2022_ks_DoubleMuMu_slxy3_lxy8_v2", "Histproduction", "ks_DoubleMuMu/");
   delete tC_22_DoubleMuMu2;

   plotFakenew_lxy c5;
   TChain* tC_22_DoubleMuMu3 = new TChain("Events");
   c5.list_dir_file(tC_22_DoubleMuMu3, "2022_data_ks_DoubleMuMu3.txt");
   cout << "Entry Data DoubleMuMu" << tC_22_DoubleMuMu3->GetEntries() << endl;
   c5.loopOverChain(tC_22_DoubleMuMu3, "Data_521_2022_ks_DoubleMuMu_slxy3_lxy8_v3", "Histproduction", "ks_DoubleMuMu/");
   delete tC_22_DoubleMuMu3;

   plotFakenew_lxy c6;
   TChain* tC_22_DoubleMuMu4 = new TChain("Events");
   c6.list_dir_file(tC_22_DoubleMuMu4, "2022_data_ks_DoubleMuMu4.txt");
   cout << "Entry Data DoubleMuMu" << tC_22_DoubleMuMu4->GetEntries() << endl;
   c6.loopOverChain(tC_22_DoubleMuMu4, "Data_521_2022_ks_DoubleMuMu_slxy3_lxy8_v4", "Histproduction", "ks_DoubleMuMu/");
   delete tC_22_DoubleMuMu4;  
 */
   return 0;
 }





