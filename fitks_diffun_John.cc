
#include "RooRealVar.h"
#include "RooProduct.h"
#include "RooDataHist.h"
#include "RooBernstein.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
//#include <nlohmann/json.hpp>

using namespace RooFit;
using namespace std;
//using json = nlohmann::json; 

vector<double> do_fit(vector<double>&param,TString option, TString name, TString ref, TString nameof, TString outname, TString pathto, TString opt_point, TString binning_type) {
  //TFile *f1 = new TFile(Form("/eos/home-c/ckar/www/plot/Fakerate/Files_UL/histo_WDY_dup_%s.root",nameof.Data()));
  TFile *f1 = new TFile(Form("%shisto_%s.root", pathto.Data(), nameof.Data()));
  TH1D *h1 = (TH1D*)f1->Get(name);
  TH1D *h2 = (TH1D*)f1->Get(ref);
  cout<<ref<<" binentry "<<h2->GetSumOfWeights()<<endl;
  cout<<name<<" binentry "<<h1->GetSumOfWeights()<<endl;
  double bin_min=0.45;
  double bin_max=0.55;
  RooRealVar m("m","",bin_min,bin_max);
  //RooRealVar m("m","",0.47,0.53);
  RooDataHist *rdh = new RooDataHist("rdh","",m,h1);
  RooDataHist *rdh_ref = new RooDataHist("rdh_ref","",m,h2);

  RooRealVar G1_mean("G1_mean","",0.4977,0.485,0.505);
  RooRealVar G1_sigma("G1_sigma","",0.003,0.001,0.01);
  RooRealVar G2_scale("G2_scale","",1.5,1.,4.5);
  RooRealVar G3_scale("G3_scale","",2.0,0.5,6);
  RooProduct G2_sigma("G2_sigma","",RooArgList(G1_sigma,G2_scale));
  RooProduct G3_sigma("G3_sigma","",RooArgList(G1_sigma,G3_scale));

  RooGaussian G1("G1","",m,G1_mean,G1_sigma);
  RooGaussian G2("G2","",m,G1_mean,G2_sigma);
  RooGaussian G3("G3","",m,G1_mean,G3_sigma);
  RooRealVar G2_fract("G2_fract","",0.2,0.0,1.0);
  RooRealVar G3_fract("G3_fract","",0.2,0.0,1.0);
  //RooAddPdf pdf_sig1 ("pdf_sig1","",RooArgSet(G2,G1),RooArgSet(G2_fract1));

  RooRealVar sig_mu("sig_mu", "sig_mu", 0.4977,0.485,0.505);
  RooRealVar sig_lambda("sig_lambda", "sig_lambda", 0.5, 0, 10);
  RooRealVar sig_gamma("sig_gamma", "sig_gamma", 0., -10, 10);
  RooRealVar sig_delta("sig_delta", "sig_delta", 1., 0, 200);
  RooJohnson Pdf_John("Pdf_John","",m,sig_mu,sig_lambda, sig_gamma, sig_delta);
  RooRealVar sigmaG_John("sigmaG_John","",0.01, 0.001, 1.);
  RooGaussian gaus_John("gaus_John","",m,sig_mu, sigmaG_John);
  RooRealVar JohnG_frac("JohnG_frac","",0.3,0.,1.0);    
  RooAddPdf massJG("massJG","mass2", RooArgList(Pdf_John,gaus_John), JohnG_frac);

  RooAddPdf pdf_sig ("pdf_sig","",RooArgSet(G3,G2,G1),RooArgSet(G2_fract,G3_fract));
 
  RooRealVar C1("C1","",0.5,0.,10.);
  RooRealVar C2("C2","",1.0,0.,10.);
  RooRealVar C3("C3","",1.0,0.,10.);
  RooRealVar C4("C4","",1.0,0.,10.);
  RooBernstein pdf_cmb("pdf_cmb","",m,RooArgList(RooConst(1.0),C1,C2));
 
  RooRealVar mean0("mean0","meam",0.495,0.485,0.505);
  RooRealVar sigma0("sigma0","m1",0.003,0.001,0.01);
  RooRealVar tail0("tail0","",5,2,20);
  RooRealVar pow0("pow0","",2,1,10);
  RooCBShape cbs0("cbs0","Signal Lineshape",m,G1_mean,sigma0,tail0,pow0);
  RooRealVar fg("fg","",.5,.1,.9);
  RooAddPdf mass_comb("mass_comb","mass1", RooArgList(cbs0,G1), fg);
  //RooAddPdf mass_comb("mass_comb","mass1", RooArgList(cbs0), fg);

  string sbsName = Form("%s",outname.Data());
  string sbsName_f = Form("%s",ref.Data());  


  cout<<"name "<<sbsName<<endl;
  RooRealVar nsig("nsig","",h2->GetSumOfWeights()*0.9,0.,h2->GetSumOfWeights());
  RooRealVar ncmb("ncmb","",h2->GetSumOfWeights()*0.2,0.,h2->GetSumOfWeights());
  RooAddPdf* model = NULL;
  //model= new RooAddPdf("model","",RooArgList(mass_comb, pdf_cmb),RooArgList(nsig, ncmb));
  // if(nameof.Contains("2016"))model= new RooAddPdf("model","",RooArgList(pdf_sig, pdf_cmb),RooArgList(nsig, ncmb));
  // else model= new RooAddPdf("model","",RooArgList(Pdf_John, pdf_cmb),RooArgList(nsig, ncmb)); 
  //  model= new RooAddPdf("model","",RooArgList(Pdf_John, pdf_cmb),RooArgList(nsig, ncmb));
  //model= new RooAddPdf("model","",RooArgList(massJG, pdf_cmb),RooArgList(nsig, ncmb));
  //model= new RooAddPdf("model","",RooArgList(pdf_sig, pdf_cmb),RooArgList(nsig, ncmb));
  model = new RooAddPdf("model", "", RooArgList(mass_comb, pdf_cmb), RooArgList(nsig, ncmb));
  //if(string::npos == sbsName.find("muid"))model= new RooAddPdf("model","",RooArgList(pdf_sig, pdf_cmb),RooArgList(nsig, ncmb));
  //else model= new RooAddPdf("model","",RooArgList(G1, pdf_cmb),RooArgList(nsig, ncmb));
  RooFitResult* fitout=model->fitTo(*rdh_ref, Save());
  fitout->Print("v");
  int retry =3;
  /*if(nameof.Contains("2016")){
    sig_gamma.setVal(0.3);
    sig_delta.setVal(1.2);
    fitout=model->fitTo(*rdh_ref, Save());
    fitout->Print("v");
  }

  while(retry>0){    
    if(fitout->status() !=0){
      cout<<"retrying"<<endl;    
      C3.setVal(1.5);
      fitout=model->fitTo(*rdh_ref, Save());
      fitout->Print("v");   
    }
    retry--;
  }*/

  cout << "Fitting" << endl;
  RooFitResult* best_fitout = NULL;
  RooFitResult* current_fitout = NULL;
  RooAddPdf* best_model = NULL;
  double best_chi2 = 1e10; //large number
  double current_chi2;
  while(retry>0){
    current_fitout=model->fitTo(*rdh_ref, Save());
    current_fitout->Print("v"); 
    RooPlot *test_frame1 = m.frame(Title(" "));
    rdh_ref->plotOn(test_frame1, Name("data"));
    model->plotOn(test_frame1, Name("fit"));

    int n_param1 = fitout->floatParsFinal().getSize();
    current_chi2 = test_frame1->chiSquare("fit", "data",n_param1);

    cout << current_chi2 << endl;
    
    if(current_chi2 < best_chi2) {
      current_chi2 = best_chi2;
      best_fitout = current_fitout;
      best_model = model;
    }
    retry --;
  }
  model = best_model;



  TCanvas *canvas = new TCanvas("canvas","", 600, 600);
  canvas->SetMargin(0.14,0.06,0.13,0.07);
  RooPlot *frame1 = m.frame(Title(" "));//outname));
  rdh_ref->plotOn(frame1, Name("data"));
  model->plotOn(frame1, Name("fit"), LineColor(6));
  model->plotOn(frame1, Components("pdf_cmb"),LineStyle(7), LineColor(kRed));  
  //model->plotOn(frame1, Components("Pdf_John"),LineStyle(7),LineColor(kGreen));
  // model->plotOn(frame1, Components("G1"),LineStyle(7),LineColor(kGreen));
  // model->plotOn(frame1, Components("G2"),LineStyle(7),LineColor(kViolet));  
  cout<<"name of "<<nameof<<endl;
  // if(nameof.Contains("2016"))model->plotOn(frame1, Components("pdf_sig"),LineStyle(7),LineColor(kGreen));  
  // else model->plotOn(frame1, Components("Pdf_John"),LineStyle(7),LineColor(kGreen));
  model->plotOn(frame1, Components("mass_comb"),LineStyle(7),LineColor(kGreen));
  model->plotOn(frame1, Components("cbs0"),LineStyle(7),LineColor(kPink));
  model->plotOn(frame1, Components("G1"),LineStyle(7),LineColor(kViolet));
  frame1->GetXaxis()->SetTitle("m(K_{s} #rightarrow #mu#pi) GeV");
  string opt_str = Form("%s",opt_point.Data());
  if (string::npos != opt_str.find("noid"))   frame1->GetXaxis()->SetTitle("m(K_{s} #rightarrow #pi^{+}#pi^{-}) GeV");
  frame1->GetXaxis()->CenterTitle(); 
  frame1->Draw();
  //canvas->SetLogy();
  TPaveText* paveText = new TPaveText(0.62,0.67,0.88,0.90,"NDC");
  paveText->SetBorderSize(0.0);
  paveText->SetFillColor(kWhite);
  paveText->SetFillStyle(0);
  paveText->SetTextSize(0.02);
  paveText->AddText(Form("Signal = %.1f #pm %.1f ", nsig.getVal(),nsig.getError()));
  paveText->AddText(Form("Background = %.1f #pm %.1f ", ncmb.getVal(),ncmb.getError()));
  paveText->AddText(Form("Mean = %.4f #pm %.6f ", G1_mean.getVal(),G1_mean.getError()));
  //paveText->AddText(Form("Sigma1 = %.4f #pm %.6f ", G1_sigma.getVal(),G1_sigma.getError()));
  paveText->Draw();
  cout<<"nsig "<<nsig.getVal()<<endl;
  cout<<"ncomb "<<ncmb.getVal()<<endl;
  string changename=Form("%s",nameof.Data());
  TString title = "#font[61]{CMS}";
  TString title2 = "#font[52]{Preliminary}";
  if(string::npos != changename.find("MC")) title2 = "#font[52]{Simulation}";
  TLatex tex;
  tex.SetTextFont(42);
  tex.SetTextSize(0.035);
  tex.SetTextAlign(11);
  tex.SetNDC();
  TString title3= Form("MC     [Muon ID = %s]",opt_point.Data());
  if(string::npos != changename.find("data") && string::npos != changename.find("2018") )title3= "#font[42]{60 fb^{-1}(13 TeV)}";
  if(string::npos != changename.find("data") && string::npos != changename.find("2017") )title3= "#font[42]{47 fb^{-1}(13 TeV)}";
  if(string::npos != changename.find("data") && string::npos != changename.find("2022") )title3= Form("#font[42]{(13.6 TeV)}     [Muon ID = %s]",opt_point.Data());

  TString title4= "bla";
  TString title_f = "bla";
  if (binning_type == "pT") {
    if(string::npos != sbsName.find("bin0")) title4=Form("p_{T} = (%.2f, %.2f)",0.0,4.0) ;
    if(string::npos != sbsName.find("bin1")) title4=Form("p_{T} = (%.2f, %.2f)",4.0,8.0) ;
    if(string::npos != sbsName.find("bin2")) title4=Form("p_{T} = (%.2f, %.2f)",8.0,12.0) ;
    if(string::npos != sbsName.find("bin3")) title4=Form("p_{T} = (%.2f, %.2f)",12.0,16.0) ;
    if(string::npos != sbsName.find("bin4")) title4=Form("p_{T} = (%.2f, %.2f)",16.0,20.0) ;
    if(string::npos != sbsName.find("bin5")) title4=Form("p_{T} = (%.2f, %.2f)",20.0,30.0) ;
    if(string::npos != sbsName.find("bin6")) title4=Form("p_{T} = (%.2f, %.2f)",30.0,50.0) ;
    if(string::npos != sbsName_f.find("allbin")) title_f=Form("p_{T} = (%.2f, %.2f)",0.0,50.0) ;
  } else {
    if(string::npos != sbsName.find("bin0")) title4=Form("L_{xy} = (%.2f, %.2f)",0.0,1.0) ;
    if(string::npos != sbsName.find("bin1")) title4=Form("L_{xy} = (%.2f, %.2f)",1.0,2.0) ;
    if(string::npos != sbsName.find("bin2")) title4=Form("L_{xy} = (%.2f, %.2f)",2.0,4.0) ;
    if(string::npos != sbsName.find("bin3")) title4=Form("L_{xy} = (%.2f, %.2f)",4.0,8.0) ;
    if(string::npos != sbsName.find("bin4")) title4=Form("L_{xy} = (%.2f, %.2f)",8.0,12.0) ;
    if(string::npos != sbsName.find("bin5")) title4=Form("L_{xy} = (%.2f, %.2f)",12.0,16.0) ;
    if(string::npos != sbsName.find("bin6")) title4=Form("L_{xy} = (%.2f, %.2f)",16.0,40.0) ;
    if(string::npos != sbsName_f.find("allbin")) title_f=Form("L_{xy} = (%.2f, %.2f)",0.0,40.0) ;
  }
  tex.SetTextSize(0.035);
  tex.DrawLatex(0.21,0.94,title);
  tex.SetTextSize(0.028);
  tex.DrawLatex(0.29,0.94,title2);
  tex.DrawLatex(0.54,0.94,title3);
  tex.DrawLatex(0.16,0.84,title_f);
  RooArgSet norm0(m);
  double chisq0 = 0.;
  int binning0= h2->GetSize()-2;
  TH1D* hist1 = (TH1D*)rdh_ref->createHistogram("hist1",m,Binning(binning0,bin_min,bin_max));//0.47,0.53));
  for(int bin=1;bin<=binning0;++bin) {
    double x = hist1->GetBinCenter(bin);
    m.setVal(x);
    double f = model->getVal(&norm0)*(nsig.getVal()+ncmb.getVal())*((bin_max-bin_min)/binning0);
    double v = hist1->GetBinContent(bin);
    double err = sqrt(v);
    if (v>0.) chisq0 += pow((v-f)/err,2);
  }
  int n_param1 = fitout->floatParsFinal().getSize();
  cout << "chi^2(by RooPlot) = " << frame1->chiSquare(n_param1) << endl;
  cout << "chi^2(by hand)    = " << chisq0/((h2->GetSize()-2)-n_param1) << endl;
  cout << "chi^2(by RooPlot) = " << frame1->chiSquare("fit", "data",n_param1)<<endl;
  string outname_g=Form("%s",nameof.Data());
  string modf= "bla";
  if(string::npos != outname_g.find("Data")) modf="Data";
  if(string::npos != outname_g.find("MC")) modf="MC";
  //int n_param1 = fitout->floatParsFinal().getSize();
  TString title5= Form("#chi^{2}/Ndf = %0.2f ",chisq0/((h2->GetSize()-2)-n_param1));
  tex.DrawLatex(0.16,0.80,title5);
  canvas->SaveAs(Form("%sCombined_%s_%s_%s.png",pathto.Data(),nameof.Data(),ref.Data(),modf.c_str()));
  canvas->SaveAs(Form("%sCombined_%s_%s_%s.pdf",pathto.Data(),nameof.Data(),ref.Data(),modf.c_str()));

  cout<<"now  bin fit"<<endl;
  nsig.setMax(h1->GetSumOfWeights());
  ncmb.setMax(h1->GetSumOfWeights());
  nsig.setVal(h1->GetSumOfWeights()*0.1);
  ncmb.setVal(h1->GetSumOfWeights()*0.9);

  C3.setVal(-3.);
  G2_scale.setConstant(true);
  G3_scale.setConstant(true);
  // G1_sigma.setConstant(true);
  //G2_fract.setConstant(true);
  //G3_scale.setConstant(true);
  //G3_fract.setConstant(true);
  string tmp_name=Form("%s",name.Data());
  string tmp_name2=Form("%s",option.Data());
  if(string::npos != tmp_name2.find("fix")){
    if(name.Contains("muid") || name.Contains("medid") || name.Contains("softid")){
      cout<<param.at(0)<<"\t"<<param.at(1)<<endl;
      
      sig_mu.setVal(param.at(0));
      sig_lambda.setVal(param.at(1));
      sig_gamma.setVal(param.at(2));
      sig_delta.setVal(param.at(3));
      sigmaG_John.setVal(param.at(4));
      JohnG_frac.setVal(param.at(5));
      G1_sigma.setVal(param.at(6));
      G2_fract.setVal(param.at(7));
      G3_fract.setVal(param.at(8));
      G2_scale.setVal(param.at(9));
      G3_scale.setVal(param.at(10));
      G1_mean.setVal(param.at(11));

      G1_sigma.setConstant(true);
      G2_fract.setConstant(true);
      G3_fract.setConstant(true);
      G2_scale.setConstant(true);
      G3_scale.setConstant(true);
      G1_mean.setConstant(true);
      sig_mu.setConstant(true);
      sig_lambda.setConstant(true);    
      sig_gamma.setConstant(true);
      sig_delta.setConstant(true);   
      sigmaG_John.setConstant(true);
      JohnG_frac.setConstant(true);
    }
  }


  fitout=model->fitTo(*rdh, Save());
  if(string::npos != tmp_name.find("Mass_bin")){
    //g1sigma
    RooRealVar g1sigma("g1sigma","",0.03);
    //g1sigma.setVal(g1sigmaparam.at(0));  
    RooRealVar g1mean("g1mean","",0.49);
    //g1mean.setVal(param.at(5));
    
    RooGaussian cons_sig("cons_sig","", G1_sigma, RooConst(g1sigma.getVal()),RooConst(0.005));
    RooGaussian cons_mean("cons_mean","",G1_mean, RooConst(g1mean.getVal()),RooConst(0.005));
    //fitout=model->fitTo(*rdh, Save(),ExternalConstraints(RooArgSet(cons_sig,cons_mean)));//,Minos(true));
  }


  fitout->Print("v");
  if(fitout->status() !=0){
    C3.setVal(-12.);
    fitout=model->fitTo(*rdh, Save());
    fitout->Print("v");
  }

  //canvas->Clear();//SetMargin(0.14,0.06,0.13,0.07);
  // gPad->SetLogy();
  TCanvas *cmas2ds=new TCanvas("cmas2ds","2d gausian1",600,600);
  TPad *pad5 = new TPad("pad5","The pad with the histograms",0,0.35,1.0,1.0);
  TPad *pad6 = new TPad("pad6","The pad with the ratio",0,0,1.0,0.35);
  pad5->Draw();
  pad6->Draw();
  pad5->cd();
  pad5->SetBottomMargin(0);
  RooPlot *frame = m.frame(Title(" "));//outname));
  rdh->plotOn(frame, Name("data"));
  model->plotOn(frame, Name("fit"));
  model->plotOn(frame, Components("pdf_cmb"),LineStyle(7), LineColor(kRed));
  model->plotOn(frame, Components("mass_comb"),LineStyle(7),LineColor(kGreen));
  model->plotOn(frame, Components("G1"),LineStyle(7),LineColor(kViolet));
  model->plotOn(frame, Components("cbs0"),LineStyle(7),LineColor(kPink));
  frame->GetXaxis()->SetTitle("m(K_{s} #rightarrow #mu#pi) GeV");
  if (string::npos != opt_str.find("noid"))   frame->GetXaxis()->SetTitle("m(K_{s} #rightarrow #pi^{+}#pi^{-}) GeV");
  frame->GetXaxis()->CenterTitle();
 
  frame->Draw();
  //gPad->SetLogy();  
  TPaveText* paveText2 = new TPaveText(0.62,0.67,0.88,0.90,"NDC");
  paveText2->SetBorderSize(0.0);
  paveText2->SetFillColor(kWhite);
  paveText2->SetFillStyle(0);
  paveText2->SetTextSize(0.02);
  paveText2->AddText(Form("Signal = %.1f #pm %.1f ", nsig.getVal(),nsig.getError()));
  paveText2->AddText(Form("Background = %.1f #pm %.1f ", ncmb.getVal(),ncmb.getError()));
  paveText2->AddText(Form("Mean = %.4f #pm %.6f ", G1_mean.getVal(),G1_mean.getError()));
  //paveText2->AddText(Form("Tail = %.4f #pm %.6f ", tail0.getVal(),tail0.getError()));
  //paveText2->AddText(Form("Power = %.4f #pm %.6f ", pow0.getVal(),pow0.getError()));
  //paveText2->AddText(Form("Sigma1 = %.4f #pm %.6f ", G1_sigma.getVal(),G1_sigma.getError()));

  paveText2->Draw();

  cout<<"nsig "<<nsig.getVal()<<endl;
  cout<<"ncomb "<<ncmb.getVal()<<endl;
  //TString title3= "MC";
  tex.SetTextSize(0.035);
  tex.DrawLatex(0.18,0.91,title);
  tex.SetTextSize(0.028);
  tex.DrawLatex(0.24,0.91,title2);
  tex.DrawLatex(0.54,0.91,title3);
  tex.DrawLatex(0.16,0.84,title4);
  RooArgSet norm(m);
  double chisq = 0.;
  int binning= h1->GetSize()-2;
  TH1D* hist = (TH1D*)rdh->createHistogram("hist",m,Binning(binning,bin_min,bin_max));//0.47,0.53));
  for(int bin=1;bin<=binning;++bin) {
    double x = hist->GetBinCenter(bin);
    m.setVal(x);
    double f = model->getVal(&norm)*(nsig.getVal()+ncmb.getVal())*((bin_max-bin_min)/binning);
    double v = hist->GetBinContent(bin);
    double err = sqrt(v);
    if (v>0.) chisq += pow((v-f)/err,2);
  }
  int n_param = fitout->floatParsFinal().getSize();
  cout << "chi^2(by RooPlot) = " << frame->chiSquare(n_param) << endl;
  cout << "chi^2(by hand)    = " << chisq/((h1->GetSize()-2)-n_param) << endl;
  cout << "chi^2(by RooPlot) = " << frame->chiSquare("fit", "data",n_param)<<endl;
 
  //int n_param = fitout->floatParsFinal().getSize();
  TString title6= Form("#chi^{2}/Ndf = %0.2f ",chisq/((h2->GetSize()-2)-n_param));//frame->chiSquare(n_param));
  tex.DrawLatex(0.16,0.80,title6);
  pad6->cd();
  pad6->SetBottomMargin(0.25);
  pad6->SetTopMargin(0);
  pad6->SetGridx(0);
  pad6->SetGridy(0);
  //cmas2ds->Clear();
  //cmas2ds->SetLogy(kFALSE);
  RooPlot* pull = m.frame(Title(" "));
  frame->Print("V");
  auto dataHist  = (RooHist*) frame->getHist("data");
  auto curve1 = (RooCurve*) frame->getObject(1);
  auto curve2 = (RooCurve*) frame->getObject(2);
  auto hresid1 =  dataHist->makePullHist(*curve1,true);

  RooExponential exp("exp", "exp", m, RooConst(1.0));
  exp.plotOn(pull,LineColor(kBlack));

  pull->GetXaxis()->SetTitle("m(K_{s} #rightarrow #mu#pi) GeV");
  if (string::npos != opt_str.find("noid"))   pull->GetXaxis()->SetTitle("m(K_{s} #rightarrow #pi^{+}#pi^{-}) GeV");
  pull->addPlotable(hresid1,"P") ;
  pull->GetXaxis()->SetTitleOffset(1.35);
  pull->GetXaxis()->SetLabelOffset(0.02);
  pull->GetXaxis()->SetTitleSize(0.08);
  pull->GetXaxis()->SetLabelSize(0.06);
  pull->GetYaxis()->SetLabelSize(0.06);
  pull->GetYaxis()->SetTitle("Pull");
  //pull->GetYaxis()->SetCenter
  pull->SetTitleSize(0.05,"Y");
  pull->GetYaxis()->CenterTitle(true);
  //ht1->GetXaxis()->SetLabelSize(0.06);
  pull->Draw();
  //canvas->SaveAs(Form("%s/pull_%s_%s.png",pathto.Data(),nameof.Data(),outname.Data()));
  //canvas->SaveAs(Form("%s/pull_%s_%s.pdf",pathto.Data(),nameof.Data(),outname.Data()));
  cmas2ds->SaveAs(Form("%s%s_%s_%s.png",pathto.Data(),nameof.Data(),outname.Data(),opt_point.Data()));
  cmas2ds->SaveAs(Form("%s%s_%s_%s.pdf",pathto.Data(),nameof.Data(),outname.Data(),opt_point.Data()));
  vector< double > fitvalue;
  fitvalue.push_back(nsig.getVal());
  fitvalue.push_back(nsig.getError());
  if(string::npos != sbsName.find("pt_bin")){
    param.push_back(sig_mu.getVal());
    param.push_back(sig_lambda.getVal());
    param.push_back(sig_gamma.getVal());
    param.push_back(sig_delta.getVal());
    param.push_back(sigmaG_John.getVal());
    param.push_back(JohnG_frac.getVal());
    param.push_back(G1_sigma.getVal());
    param.push_back(G2_fract.getVal());
    param.push_back(G3_fract.getVal());
    param.push_back(G2_scale.getVal());
    param.push_back(G3_scale.getVal());
    param.push_back(G1_mean.getVal());

  }
  canvas->Clear();
  cmas2ds->Clear();
  cout << endl;
  return fitvalue;
}
double dRatio(double a, double ae, double b, double be) {
    return TMath::Sqrt(((ae*ae)/(b*b)) + ((a*a*be*be)/(b*b*b*b)));
}

TH1D* playfit(TString which_file, TString path, TString option, TString opt_point, TString binning) {
  
  cout << "Playing fit with option: " << option << " and opt point " << opt_point << endl;

  // double var[]={0.,4.0,8.0,12.0,16.0,20.0, 30., 50.};
  // int n_bin=7;
  double var[8];
  if (binning == "pT") {
    var[0] = 0.0;
    var[1] = 4.0;
    var[2] = 8.0;
    var[3] = 12.0;
    var[4] = 16.0;
    var[5] = 20.0;
    var[6] = 30.0;
    var[7] = 50.0;
  } else {
    var[0] = 0.0;
    var[1] = 1.0;
    var[2] = 2.0;
    var[3] = 4.0;
    var[4] = 8.0;
    var[5] = 12.0;
    var[6] = 16.0;
    var[7] = 40.0;
  }
  int n_bin=6;
  TH1D* hv0_fakeid_pt = new TH1D("hv0_fakeid_pt_", ";p_{T}[GeV];Fakerate", 6,var );
  vector<double> meanval;
  vector<double> errval;
  for(int bin=0;bin<n_bin;++bin){
    vector<double> paramval;
    TString name = Form("pt_bin%d", bin);
    vector<double> pt_yield =do_fit(paramval,option, Form("hv0_Mass_bin%d",bin),"hv0_Mass_allbin", which_file, name, path,"noid", binning);
    cout<<" muon id fitting"<<endl;
    cout<<paramval.at(0) <<"\t"<<paramval.at(1)<<endl;
    vector<double> pt_yield_muid;
    string opt_point_=Form("%s",opt_point.Data());

    if(opt_point.Contains("bdt20")){
      name = Form("pt_muid20_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_20_bin%d",bin),"hv0_Mass_muid_20_allbin", which_file, name, path, opt_point, binning);
    }              
    else if(opt_point.Contains("bdt30")){
      name = Form("pt_muid30_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_30_bin%d",bin),"hv0_Mass_muid_30_allbin", which_file, name, path, opt_point, binning);
    }else if(opt_point.Contains("bdt40")){
      name = Form("pt_muid40_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_40_bin%d",bin),"hv0_Mass_muid_40_allbin", which_file, name, path, opt_point, binning);
    }else if(opt_point.Contains("bdt45")){
      name = Form("pt_muid45_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_45_bin%d",bin),"hv0_Mass_muid_45_allbin", which_file, name, path, opt_point, binning);
    }else if(opt_point.Contains("bdt50")){
      name = Form("pt_muid50_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_50_bin%d",bin),"hv0_Mass_muid_50_allbin", which_file, name, path, opt_point, binning);
    }
    else if(opt_point.Contains("bdt55")){
      name = Form("pt_muid55_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_55_bin%d",bin),"hv0_Mass_muid_55_allbin", which_file, name, path, opt_point, binning);
    }    
    else if(opt_point.Contains("bdt60")){
      name = Form("pt_muid60_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_60_bin%d",bin),"hv0_Mass_muid_60_allbin", which_file, name, path, opt_point, binning);
    }
    else if(opt_point.Contains("softid")){
      name = Form("pt_softid_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_softid_bin%d",bin),"hv0_Mass_softid_allbin", which_file, name, path, opt_point, binning);
    }
    else if(opt_point.Contains("mediumid")){
      name = Form("pt_medid_bin%d", bin);
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_medid_bin%d",bin),"hv0_Mass_medid_allbin", which_file, name, path, opt_point, binning);
    }
    else{
      name = Form("pt_muid_bin%d", bin);    
      pt_yield_muid = do_fit(paramval,option,Form("hv0_Mass_muid_bin%d",bin),"hv0_Mass_muid_allbin", which_file, name, path, opt_point, binning);
    }
        

    cout<<" fake rate in ptbin "<<bin<<" => "<<pt_yield_muid.at(0)<<" /t "<<pt_yield.at(0)<<"\t "<<pt_yield_muid.at(0)/pt_yield.at(0)<<endl;
    double error = 0;
    if(pt_yield_muid.at(0)>0 || pt_yield_muid.at(1)>0)error = dRatio(pt_yield_muid.at(0), pt_yield_muid.at(1), pt_yield.at(0), pt_yield.at(1));
    cout<<" Error "<<error<<endl;
    if(error>100)error = 1;
    hv0_fakeid_pt->SetBinContent(bin+1, pt_yield_muid.at(0)/pt_yield.at(0));
    hv0_fakeid_pt->SetBinError(bin+1, error);
    meanval.push_back(pt_yield_muid.at(0)/pt_yield.at(0));
    errval.push_back(error);
    if(pt_yield_muid.at(0)<2.0){
      hv0_fakeid_pt->SetBinContent(bin+1,0.0);
      hv0_fakeid_pt->SetBinError(bin+1,0.0);  
    }
    pt_yield.clear();
    pt_yield_muid.clear();
    paramval.clear();
  }
  cout<<"bin1 "<<hv0_fakeid_pt->GetBinContent(1)<<"\t err "<<hv0_fakeid_pt->GetBinError(1)<<endl;;
  cout<<"bin2 "<<hv0_fakeid_pt->GetBinContent(2)<<"\t err "<<hv0_fakeid_pt->GetBinError(2)<<endl;;
  cout<<"bin3 "<<hv0_fakeid_pt->GetBinContent(3)<<"\t err "<<hv0_fakeid_pt->GetBinError(3)<<endl;;
  cout<<"bin4 "<<hv0_fakeid_pt->GetBinContent(4)<<"\t err "<<hv0_fakeid_pt->GetBinError(4)<<endl;;
  cout<<"bin5 "<<hv0_fakeid_pt->GetBinContent(5)<<"\t err "<<hv0_fakeid_pt->GetBinError(5)<<endl;;
  cout<<"bin6 "<<hv0_fakeid_pt->GetBinContent(6)<<"\t err "<<hv0_fakeid_pt->GetBinError(6)<<endl;;
  // cout<<"bin7 "<<hv0_fakeid_pt->GetBinContent(7)<<"\t err "<<hv0_fakeid_pt->GetBinError(7)<<endl;;

   TCanvas* c1=new TCanvas("c1","", 600,600);
   c1->SetMargin(0.17,0.06,0.13,0.07);
   c1->Clear();
   gStyle->SetOptStat(0);

   TPaveText* paveText2 = new TPaveText(0.17,0.77,0.53,0.90,"NDC");
   paveText2->SetBorderSize(0.0);
   paveText2->SetFillColor(kWhite);
   paveText2->SetFillStyle(0);
   paveText2->SetTextSize(0.02);
   if (binning == "pT") {
     paveText2->AddText(Form("p_{T} [4.0, 8.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(2), hv0_fakeid_pt->GetBinError(2)));
     paveText2->AddText(Form("p_{T} [8.0, 12.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(3), hv0_fakeid_pt->GetBinError(3)));
     paveText2->AddText(Form("p_{T} [12.0, 16.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(4), hv0_fakeid_pt->GetBinError(4)));
   } else {
     paveText2->AddText(Form("L_{xy} [0.0, 1.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(1), hv0_fakeid_pt->GetBinError(1)));
     paveText2->AddText(Form("L_{xy} [1.0, 2.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(2), hv0_fakeid_pt->GetBinError(2)));
     paveText2->AddText(Form("L_{xy} [2.0, 4.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(3), hv0_fakeid_pt->GetBinError(3)));
   }

   TPaveText* pave2 = new TPaveText(0.56,0.77,0.98,0.90,"NDC");
   pave2->SetBorderSize(0.0);
   pave2->SetFillColor(kWhite);
   pave2->SetFillStyle(0);
   pave2->SetTextSize(0.02);
   if (binning == "pT") {
     pave2->AddText(Form("p_{T} [16.0, 20.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(5), hv0_fakeid_pt->GetBinError(5)));
     pave2->AddText(Form("p_{T} [20.0, 30.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(6), hv0_fakeid_pt->GetBinError(6)));
     //pave2->AddText(Form("p_{T} [30.0, 50.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(7), hv0_fakeid_pt->GetBinError(7)));
   }
   else {
    pave2->AddText(Form("L_{xy} [4.0, 8.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(4), hv0_fakeid_pt->GetBinError(4)));
    pave2->AddText(Form("L_{xy} [8.0, 12.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(5), hv0_fakeid_pt->GetBinError(5)));
    pave2->AddText(Form("L_{xy} [12.0, 16.0] = %.5f #pm %.5f ", hv0_fakeid_pt->GetBinContent(6), hv0_fakeid_pt->GetBinError(6)));
   }

   TString filenamechange="bla";
   string tmpstring=Form("%s",which_file.Data());   
   hv0_fakeid_pt->SetMinimum(0.0);
   // if(opt_point.Contains("mediumid"))hv0_fakeid_pt->SetMaximum(0.01);
   // else if (opt_point.Contains("softid"))hv0_fakeid_pt->SetMaximum(0.01);
   // else
   hv0_fakeid_pt->SetMaximum(0.005);
   // string checkname=Form("%s",which_file.Data());
   // if(string::npos != checkname.find("Data")){
   //   hv0_fakeid_pt->SetMinimum(0.0005);
   //   hv0_fakeid_pt->SetMaximum(0.007);
   // }
   TLatex tex;
   tex.SetTextFont(42);
   tex.SetTextSize(0.035);
   tex.SetTextAlign(11);
   tex.SetNDC();
   TString title = "#font[61]{CMS}";
   TString title2= "#font[52]{Preliminary}";
   TString title3 = "2022";
   if(which_file.Contains("MC")) title2 = "#font[52]{Simulation}";
   if(which_file.Contains("2018"))title3 = "2018";
   if(which_file.Contains("2017"))title3 = "2017";
   TString title4 = "muon mva = bmm4";
   if(opt_point.Contains("bdt20"))title4="muon mva>0.20";
   if(opt_point.Contains("bdt30"))title4="muon mva>0.30";
   if(opt_point.Contains("bdt40"))title4="muon mva>0.40";
   if(opt_point.Contains("bdt45"))title4="muon mva>0.45";
   if(opt_point.Contains("bdt50"))title4="muon mva>0.50";
   if(opt_point.Contains("bdt55"))title4="muon mva>0.55";
   if(opt_point.Contains("bdt60"))title4="muon mva>0.60";
   if(opt_point.Contains("softid"))title4="muon mva = looseid";
   if(opt_point.Contains("mediumid"))title4="muon mva = mediumid";
   
   if(string::npos != tmpstring.find("Data") && (string::npos != tmpstring.find("2018")))filenamechange="Data_516_2018_ks";
   if(string::npos != tmpstring.find("MC") && (string::npos != tmpstring.find("2018")))filenamechange="MC_516_2018_ks";

   if(string::npos != tmpstring.find("Data") && (string::npos != tmpstring.find("2017")))filenamechange="Data_516_2017_ks";
   if(string::npos != tmpstring.find("MC") && (string::npos != tmpstring.find("2017")))filenamechange="MC_516_2017_ks";

   if(string::npos != tmpstring.find("Data") && (string::npos != tmpstring.find("2016")))filenamechange="Data_516_2016_ks";
   if(string::npos != tmpstring.find("MC") && (string::npos != tmpstring.find("2016")))filenamechange="MC_516_2016_ks";

   //if(string::npos != tmpstring.find("momentum"))hv0_fakeid_pt->SetTitle("Fake Rate from P bin;P; Fake rate");
   //if(string::npos == tmpstring.find("momentum"))hv0_fakeid_pt->SetTitle("Fake Rate from p_{T} bin;p_{T}; Fake rate");
   hv0_fakeid_pt->GetYaxis()->SetTitleOffset(2.3);
   hv0_fakeid_pt->GetXaxis()->CenterTitle();
   hv0_fakeid_pt->GetYaxis()->CenterTitle();
   hv0_fakeid_pt->Draw();
   paveText2->Draw();
   pave2->Draw();
   tex.DrawLatex(0.18,0.94,title);
   tex.SetTextSize(0.028);
   tex.DrawLatex(0.26,0.94,title2);
   tex.DrawLatex(0.80,0.94,title3);
   tex.DrawLatex(0.45,0.94,title4);

   // c1->SaveAs(Form("%splayV0-ks_kin_pt_muid_byfitting_%s.pdf",path.Data(),which_file.Data()));
   // c1->SaveAs(Form("%splayV0-ks_kin_pt_muid_byfitting_%s.png",path.Data(),which_file.Data()));

   c1->SaveAs(Form("%splayV0-ks_kin_%s_muid%s_byfitting_%s.pdf",path.Data(),binning.Data(),opt_point.Data(),filenamechange.Data()));
   c1->SaveAs(Form("%splayV0-ks_kin_%s_muid%s_byfitting_%s.png",path.Data(),binning.Data(),opt_point.Data(),filenamechange.Data()));

   return hv0_fakeid_pt;
}
std::map<std::string, double> overlay(TH1D* h1,TH1D* h2, TH1D* h3, TString label1, TString label2, TString label3, 
                                      TString path, TString opt_point, TString which_file, TString binning){
  std::map<std::string, double>   nEff;
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
  TString title4 = "muon mva = bmm4";
  if(opt_point.Contains("bdt20"))title4="muon mva>0.20";
  if(opt_point.Contains("bdt30"))title4="muon mva>0.30";
  if(opt_point.Contains("bdt40"))title4="muon mva>0.40";
  if(opt_point.Contains("bdt45"))title4="muon mva>0.45";
  if(opt_point.Contains("bdt50"))title4="muon mva>0.50";
  if(opt_point.Contains("bdt55"))title4="muon mva>0.55";
  if(opt_point.Contains("bdt60"))title4="muon mva>0.60";
  if(opt_point.Contains("softid"))title4="muon mva = looseid";
  if(opt_point.Contains("mediumid"))title4="muon mva = mediumid";

  TCanvas* c1=new TCanvas("c1","", 600,600);
  c1->SetMargin(0.17,0.06,0.13,0.07);
  c1->Clear();
  gStyle->SetOptStat(0);
  if (binning == "pT") {
    h1->SetTitle(Form(";%s[GeV];Fakerate", binning.Data()));
  } else {
    h1->SetTitle(Form(";%s[cm];Fakerate", binning.Data()));
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
  tex.DrawLatex(0.45,0.94,title4);
  c1->SaveAs(Form("%splayV0-ks_kin_%s_muid%s_%s_overlay.pdf",path.Data(),binning.Data(),opt_point.Data(),which_file.Data()));
  c1->SaveAs(Form("%splayV0-ks_kin_%s_muid%s_%s_overlay.png",path.Data(),binning.Data(),opt_point.Data(),which_file.Data()));

  
  return nEff;

}

void fitks_loop(TString syear, TString binning){
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

   ////Final set of files
   
   TString path = "ks/";
   TString whichfile1 = Form("data_parking_522_2022_ks_%s", binning.Data());
   TString label1 = "Parking Data";
   TString whichfile2 = Form("data_egamma_522_2022_ks_trigger_%s", binning.Data());
   TString label2 = "Egamma Data";
   TString whichfile3 = Form("MC_DY_522_2022_ks_%s", binning.Data());
   TString label3 = "DY, W, TT";


   std::map<std::string, std::map<std::string, double>> flatDMap;
   std::string fJsonFileName;
   std::ofstream fJson;
   fJsonFileName = Form("pion_fakerate_%s.json",syear.Data());
   fJson.open(fJsonFileName.c_str(), ios::app);

   TString option = "fix";
   TString opt_point = "default";
  
   int num_opt_points = 8;
   TString opt_points[num_opt_points];
   TH1D* h1_arr[num_opt_points], h2_arr[num_opt_points], h3_arr[num_opt_points];
   opt_points[0] = "default";
   opt_points[1] = "bdt20";
   opt_points[2] = "bdt30";
   opt_points[3] = "bdt40";
   opt_points[4] = "bdt45";
   opt_points[5] = "bdt50";
   opt_points[6] = "mediumid";
   opt_points[7] = "softid";

   //TH1D h1, h2, h3;
   std::map<std::string, double> bmm4;
   std::string tag;
   for (int i=0; i<8; i++) {
      cout << "BDT is " << opt_points[i] << endl;
      //if (opt_points[i]=="bdt20" || opt_points[i]=="bdt30" || opt_points[i]=="bdt40" || opt_points[i]=="bdt45" || opt_points[i]=="bdt50") continue; //TODO: remove
      TH1D* h1_arr = playfit(whichfile1, path, option, opt_points[i], binning);
      TH1D* h2_arr = playfit(whichfile2, path, option, opt_points[i], binning);
      TH1D* h3_arr = playfit(whichfile3, path, option, opt_points[i], binning);
      bmm4 = overlay( h1_arr, h2_arr, h3_arr, label1, label2, label3, path, opt_points[i],syear, binning); 
      tag = Form("%s_%s",opt_points[i].Data(),syear.Data());
      flatDMap[tag]=bmm4;
   }

   
   
   //json effJson(flatDMap);
   //fJson << effJson.dump(4) << endl;
   //fJson.close();
   
}
void fitks_diffun_John(){
  
  // fitks_loop("2016");
  // fitks_loop("2017");
  // fitks_loop("2018");
  //fitks_loop("2022", "lxy");
  fitks_loop("2022", "pT");
  fitks_loop("2022", "lxy");
}
