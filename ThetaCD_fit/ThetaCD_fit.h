#ifndef _THETACD_FIT_H_
#define _THETACD_FIT_H_

#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TPaveLabel.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TPaveText.h>
#include <TInterpreter.h>
#include <TStyle.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLegend.h>

TFile* f[3];        //ROOT files with data and MC simulation results
TFile* th_q[50];    //ROOT files with ThetaCD distributions
TFile* MyFile;      //create new file.root with ThetaCD distributions

Double_t Integral_DATA;
Double_t Integral_signal;
Double_t Integral_bkgnd;

Int_t sign_param,bkgnd_param,total_param;

Double_t begin_fit_sign,end_fit_sign;
Double_t begin_draw_sign,end_draw_sign;

Double_t begin_fit_bkgnd,end_fit_bkgnd;
Double_t begin_draw_bkgnd,end_draw_bkgnd;

Double_t begin_fit_total,end_fit_total;
Double_t begin_draw_total,end_draw_total;

TF1 *signalFunct[20],*signalFunct_end[20];
TF1 *backgroundFunct[20],*backgroundFunct_end[20];
TF1 *totalFunct[20];

TH1D *hThetaCD_DATA_Q,*hThetaCD_DATA_exp;  //experimental data from file.root

TH1D* hThetaCD_DATA = new TH1D("hThetaCD_DATA","",180,0.,180.);             //experimental signal+bkgnd
TH1D* hThetaCD_DATA_sign = new TH1D("hThetaCD_DATA_sign","",180,0.,180.);   //experimental signal
TH1D* hThetaCD_DATA_bkgnd = new TH1D("hThetaCD_DATA_bkgnd","",180,0.,180.); //experimental background

TH1D* hThetaCD_DATA_sign_cut = new TH1D("hThetaCD_DATA_sign_cut","",180,0.,180.);

TH1D* hThetaCD_DATA_i[50];          //experimental signal+bkgnd for i'th bin
TH1D* hThetaCD_DATA_sign_i[50];     //experimental signal for i'th bin
TH1D* hThetaCD_DATA_bkgnd_i[50];    //background for i'th bin

TH1D* hThetaCD_MC  = new TH1D("hThetaCD_MC","",180,0.,180.);    //simulated events
TH1D* hThetaCD_MC_i[50];                                        //simulated events for i'th bin

TH1D* hThetaCD_MC_qf;       //pd->ppn_spec simulated events from file.root
TH1D* hThetaCD_MC_pd;       //pd->pd simulated events from file.root

TH1D* hExcessEnergy_DATA = new TH1D("hExcessEnergy_DATA","",40,-70.,30.);   //experimental excess energy distribution
TH1D* hExcessEnergy_sign = new TH1D("hExcessEnergy_sign","",40,-70.,30.);   //signal excess energy distribution
TH1D* hExcessEnergy_bkgnd = new TH1D("hExcessEnergy_bkgnd","",40,-70.,30.); //background excess energy distribution
TH1D* hExcessEnergy_diff = new TH1D("hExcessEnergy_diff","",40,-70.,30.);   //hExcessEnergy_DATA - hExcessEnergy_sign
TH1D* hExcessEnergy_exp;    //experimental events from file

#endif
