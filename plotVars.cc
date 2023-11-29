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

#include "TLatex.h"
#include "TAxis.h"
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

#include "RooRealVar.h"
#include "RooProduct.h"
#include "RooDataHist.h"
#include "RooBernstein.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"

using namespace RooFit;
using namespace std;

class plotVars {
    public:
        TH1D* make_histogram(TChain *tC, TString var_name);
        void list_dir(TChain *, const char *);
        void list_dir_file(TChain *, string ob = "bla");
        void overlay(TH1D* h1,TH1D* h2, TH1D* h3, TString label1, TString label2, TString label3, 
                                      TString path, TString opt_point, TString which_file, TString binning);

        plotVars();
        ~plotVars();
};

plotVars::~plotVars(void)
{
  // destructor
}
plotVars::plotVars(void)
{
  // constructor
}

TH1D* plotVars::make_histogram(TChain *tC, TString var_name) {
    int NBINS = 50;
    int bin_min = 5;
    int bin_max = 30;

    TH1D *hist = new TH1D(Form("hist_%s", var_name.Data()), Form("hist_%s", var_name.Data()), NBINS, bin_min, bin_max);

    Float_t var[150];
    UInt_t nV0;
    Float_t V0_trk1_pt[150], V0_trk2_pt[150];
    Int_t V0_trk1_mu_index[150], V0_trk2_mu_index[150];
    Float_t V0_kin_slxy[150];
    //Float_t V0_kin_lxy[150];
    Float_t V0_trk1_sip[150], V0_trk2_sip[150];
    Float_t V0_kin_sipPV[150];
    Float_t V0_kin_vtx_chi2dof[150];
    Float_t V0_kin_cosAlphaXY[150];
    Float_t V0_doca[150];
    Bool_t Muon_softMvaId[150];
    Float_t V0_kin_pt[150];

    TBranch *b_var;
    TBranch *b_nV0;
    TBranch *b_V0_trk1_pt;
    TBranch *b_V0_trk2_pt;
    TBranch *b_V0_trk1_mu_index;
    TBranch *b_V0_trk2_mu_index;
    //TBranch *b_V0_kin_lxy;
    TBranch *b_V0_kin_slxy;
    TBranch *b_V0_trk1_sip;
    TBranch *b_V0_trk2_sip;
    TBranch *b_V0_kin_sipPV;
    TBranch *b_V0_kin_vtx_chi2dof;
    TBranch *b_V0_kin_cosAlphaXY;
    TBranch *b_V0_doca;
    TBranch *b_Muon_softMvaId;
    TBranch *b_V0_kin_pt;


    tC->SetBranchAddress("nks", &nV0, &b_nV0);
    tC->SetBranchAddress("ks_trk1_pt", V0_trk1_pt, &b_V0_trk1_pt);
    tC->SetBranchAddress("ks_trk2_pt", V0_trk2_pt, &b_V0_trk2_pt);
    tC->SetBranchAddress("ks_trk1_mu_index", V0_trk1_mu_index, &b_V0_trk1_mu_index);
    tC->SetBranchAddress("ks_trk2_mu_index", V0_trk2_mu_index, &b_V0_trk2_mu_index);
    //tC->SetBranchAddress("ks_kin_lxy", V0_kin_lxy, &b_V0_kin_lxy);
    tC->SetBranchAddress("ks_kin_slxy", V0_kin_slxy, &b_V0_kin_slxy);
    tC->SetBranchAddress("ks_trk1_sip", V0_trk1_sip, &b_V0_trk1_sip);
    tC->SetBranchAddress("ks_trk2_sip", V0_trk2_sip, &b_V0_trk2_sip);
    tC->SetBranchAddress("ks_kin_sipPV", V0_kin_sipPV, &b_V0_kin_sipPV);\
    tC->SetBranchAddress("ks_kin_vtx_chi2dof", V0_kin_vtx_chi2dof, &b_V0_kin_vtx_chi2dof);
    tC->SetBranchAddress("ks_kin_cosAlphaXY", V0_kin_cosAlphaXY, &b_V0_kin_cosAlphaXY);
    tC->SetBranchAddress("ks_doca", V0_doca, &b_V0_doca);
    tC->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
    tC->SetBranchAddress("ks_kin_pt", V0_kin_pt, &b_V0_kin_pt);
    tC->SetBranchAddress(var_name.Data(), var, &b_var);

    Long64_t n = tC->GetEntries();
    for (unsigned int i = 0; i < n; ++i) {
        
        Long64_t localEntry = tC->LoadTree(i);

        b_nV0->GetEntry(localEntry);

        float lxy_cut = 100;

        for (unsigned int bhad = 0; bhad < nV0; ++bhad) {

            b_var->GetEntry(localEntry);
            b_V0_trk1_pt->GetEntry(localEntry);
            b_V0_trk2_pt->GetEntry(localEntry);
            b_V0_trk1_mu_index->GetEntry(localEntry);
            b_V0_trk2_mu_index->GetEntry(localEntry);
            //cout << "here" << endl;
            //b_V0_kin_lxy->GetEntry(localEntry);
            //cout << "here3" << endl;
            b_V0_kin_slxy->GetEntry(localEntry);
            b_V0_trk1_sip->GetEntry(localEntry);
            b_V0_trk2_sip->GetEntry(localEntry);
            b_V0_kin_sipPV->GetEntry(localEntry);
            b_V0_kin_vtx_chi2dof->GetEntry(localEntry);
            b_V0_kin_cosAlphaXY->GetEntry(localEntry);
            b_V0_doca->GetEntry(localEntry);
            b_Muon_softMvaId->GetEntry(localEntry);
            b_V0_kin_pt->GetEntry(localEntry);
            
            bool cut_pre = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0)) && abs(var[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3);
            bool cut_pre_mid = (((V0_trk1_pt[bhad] > 4 && V0_trk2_pt[bhad] > 1.0 && V0_trk1_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk1_mu_index[bhad]]) || (V0_trk2_pt[bhad] > 4 && V0_trk1_pt[bhad] > 1.0 && V0_trk2_mu_index[bhad] >= 0 && Muon_softMvaId[V0_trk2_mu_index[bhad]])) && abs(var[bhad]) < lxy_cut && V0_kin_vtx_chi2dof[bhad] < 3 && V0_kin_cosAlphaXY[bhad] > 0.999 && V0_doca[bhad] < 0.004 && V0_kin_sipPV[bhad] < 3 && V0_trk1_sip[bhad] > 5 && V0_trk2_sip[bhad] > 5 && abs(V0_kin_slxy[bhad]) > 3);
            

            //if (cut_pre_mid) hist->Fill(var[bhad]);
            if (cut_pre) hist->Fill(V0_kin_pt[bhad]);
        }
    }

    TCanvas *canvas = new TCanvas("canvas","", 600, 600);
    canvas->SetMargin(0.14,0.06,0.13,0.07);

    RooRealVar roo_var("var","var",bin_min,bin_max);
    RooDataHist *rdh = new RooDataHist("rdh","",roo_var,hist);
    RooPlot *frame1 = roo_var.frame(Title(" "));

    rdh->plotOn(frame1, Name("data"));

    frame1->GetXaxis()->SetTitle(var_name);
    frame1->GetXaxis()->CenterTitle(); 
    frame1->Draw();

    TString title = "#font[61]{CMS}";
    TString title2 = "#font[52]{Preliminary}";
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.SetTextAlign(11);
    tex.SetNDC();
    tex.SetTextSize(0.035);
    tex.DrawLatex(0.21,0.94,title);
    tex.SetTextSize(0.028);
    tex.DrawLatex(0.29,0.94,title2);

    canvas->SaveAs(Form("plot_%s_test.png", var_name.Data()));

    return hist;

}

void plotVars::list_dir(TChain *tC, const char *path)
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


void plotVars::list_dir_file(TChain *tC, string filename)
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

void plotVars::overlay(TH1D* h1,TH1D* h2, TH1D* h3, TString label1, TString label2, TString label3, 
                                      TString path, TString opt_point, TString which_file, TString binning){
  
  float n_h1 = h1->Integral();
  float n_h2 = h2->Integral();
  float n_h3 = h3->Integral();
  h1->Scale(1/n_h1);
  h2->Scale(1/n_h2);
  h3->Scale(1/n_h3);
  
  TLatex tex;
  tex.SetTextFont(42);
  tex.SetTextSize(0.035);
  tex.SetTextAlign(11);
  tex.SetNDC();
  TString title = "#font[61]{CMS}";
  TString title2= "#font[52]{Preliminary}";
  TString title3 = "2022";
  if(which_file.Contains("MC")) title2 = "#font[52]{Simulation}";
  if(which_file.Contains("2017"))title3 = "2017";
  if(which_file.Contains("2016"))title3 = "2016";

  TCanvas* c1=new TCanvas("c1","", 600,600);
  c1->SetMargin(0.17,0.06,0.13,0.07);
  c1->Clear();
  gStyle->SetOptStat(0);
  if (binning=="lxy") {
    h1->SetTitle(Form(";%s[cm];Events (normalized)", binning.Data()));
  } else {
    h1->SetTitle(Form(";%s[GeV];Events (normalized)", binning.Data()));
  }
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h3->SetLineColor(kGreen);
  h1->SetMarkerColor(kRed);
  h2->SetMarkerColor(kBlue);
  h3->SetMarkerColor(kGreen);
  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(22);
  h3->SetMarkerStyle(22);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");

  TLegend* leg =new TLegend(0.77,0.78,0.9,0.9);
  //gStyle->SetLegendTextSize(0.043);
  leg->AddEntry(h1, label1, "ep");
  leg->AddEntry(h2, label2, "ep");
  leg->AddEntry(h3, label3, "ep");
  leg->Draw();
  tex.DrawLatex(0.18,0.94,title);
  tex.SetTextSize(0.028);
  tex.DrawLatex(0.26,0.94,title2);
  tex.DrawLatex(0.80,0.94,title3);

  TPaveText* paveText = new TPaveText(0.17+.25,0.77,0.53+.25,0.90,"NDC");
  //TPaveText* paveText = new TPaveText(0.17,0.77,0.53,0.90,"NDC");
  paveText->SetBorderSize(0.0);
  paveText->SetFillColor(kWhite);
  paveText->SetFillStyle(0);
  paveText->SetTextSize(0.02);
  paveText->AddText(Form("%s: n_events = %.0f", label1.Data(), n_h1));
  paveText->AddText(Form("%s: n_events = %.0f", label2.Data(), n_h2));
  paveText->AddText(Form("%s: n_events = %.0f", label3.Data(), n_h3));
  paveText->Draw();


  c1->SaveAs("overlay.png");
}

int main() {
    plotVars c1;
    TChain *tC1 = new TChain("Events");
    c1.list_dir_file(tC1, "2022_data_ks_parking.txt");
    TH1D* h1 = c1.make_histogram(tC1, "ks_kin_lxy");
    TString label1 = "EGamma Parking";

    plotVars c2;
    TChain *tC2 = new TChain("Events");
    c2.list_dir_file(tC2, "2022_data_ks_egamma.txt");
    TH1D* h2 = c2.make_histogram(tC2, "ks_kin_lxy");
    TString label2 = "Egamma";

    plotVars c3;
    TChain *tC3 = new TChain("Events");
    c3.list_dir_file(tC3, "2022_MC_DY.txt");
    TH1D* h3 = c3.make_histogram(tC3, "ks_kin_lxy");
    TString label3 = "MC";

    c3.overlay(h1, h2, h3, label1, label2, label3, "blah", "blah", "2022", "pT");

    return 0;
}