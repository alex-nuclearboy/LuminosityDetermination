/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2017-07
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-05
***********************************************/

//Macro used to merge files with fitted ThetaCD distributions

#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TVector3.h>
#include <TLorentzVector.h>
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
#include "ThetaCD_fit.h"

void ThetaCD_merge() {

    cout<<fixed;
    cout.precision(2);

    //open the ROOT file with data
    f[0] = new TFile("../input/DATA_ppn_qf_offset.root","READ");
    f[1] = new TFile("../input/MC-ppn_qf-PARIS-x6.root","READ");
    f[2] = new TFile("../input/MC-pd-momcut.root","READ");

    f[0]->cd("Histograms/DATA_lev4");
    hThetaCD_DATA_exp = (TH1D*)gDirectory->Get("hThetaCD_lev4");
    //hThetaCD_DATA_exp->Rebin(2);
    f[0]->cd("Histograms/DATA_lev5");
    hExcessEnergy_exp = (TH1D*)gDirectory->Get("hQ_lev5");

    f[1]->cd("Histograms/DATA_lev4");
    hThetaCD_MC_qf = (TH1D*)gDirectory->Get("hThetaCD_lev4");
    //hThetaCD_MC_qf->Rebin(2);

    f[2]->cd("Histograms/DATA_lev4");
    hThetaCD_MC_pd = (TH1D*)gDirectory->Get("hThetaCD_lev4");
    //hThetaCD_MC_pd->Rebin(2);

    for (Int_t i = 1; i < 41; i++) {

        th_q[i] = new TFile(Form("files/ThetaCD_Q%d.root",i),"READ");
        th_q[i]->cd();

        hThetaCD_DATA_i[i] = (TH1D*)gDirectory->Get(Form("hThetaCD_DATA_Q%d",i));
        hThetaCD_DATA_sign_i[i] = (TH1D*)gDirectory->Get(Form("hThetaCD_DATA_sign_Q%d",i));
        hThetaCD_DATA_bkgnd_i[i] = (TH1D*)gDirectory->Get(Form("hThetaCD_DATA_bckr_Q%d",i));

        hThetaCD_MC_i[i] = (TH1D*)gDirectory->Get(Form("hThetaCD_MC_Q%d",i));

        hThetaCD_DATA->Add(hThetaCD_DATA_i[i]);
        hThetaCD_DATA_sign->Add(hThetaCD_DATA_sign_i[i]);
        hThetaCD_DATA_bkgnd->Add(hThetaCD_DATA_bkgnd_i[i]);

        hThetaCD_MC->Add(hThetaCD_MC_i[i]);

        Integral_DATA = hThetaCD_DATA_i[i]->Integral();
        hExcessEnergy_DATA->SetBinContent(i,Integral_DATA);

        Integral_signal = hThetaCD_DATA_sign_i[i]->Integral();
        hExcessEnergy_sign->SetBinContent(i,Integral_signal);

        Integral_bkgnd = hThetaCD_DATA_bkgnd_i[i]->Integral();
        hExcessEnergy_bkgnd->SetBinContent(i,Integral_bkgnd);

    }

    cout<<"background = "<<(Integral_bkgnd/Integral_DATA)*100<<"%"<<endl;
    cout<<"signal = "<<(Integral_signal/Integral_DATA)*100<<"%"<<endl;
    cout<<"bacground/signal = "<<(Integral_bkgnd/Integral_signal)*100<<"%"<<endl;

    Double_t bin_content_cut;

    for(Int_t l = 1; l < 181; l++) {
        bin_content_cut = hThetaCD_DATA_sign->GetBinContent(l);

        if((l > 40) && (l <= 100)) {
            hThetaCD_DATA_sign_cut->SetBinContent(l,bin_content_cut);
        }
        else {
            hThetaCD_DATA_sign_cut->SetBinContent(l,0.);
        }
    }

    hExcessEnergy_diff = (TH1D*)hExcessEnergy_DATA->Clone("hExcessEnergy_diff");
    hExcessEnergy_diff->Add(hExcessEnergy_sign,-1);

    //create new file.root
    MyFile = new TFile("../output/files/ThetaCD_merged.root","RECREATE");

    MyFile->cd();
    hThetaCD_DATA->Write("");
    hThetaCD_DATA_sign->Write("");
    hThetaCD_DATA_bkgnd->Write("");
    hExcessEnergy_DATA->Write("");
    hExcessEnergy_sign->Write("");
    hExcessEnergy_bkgnd->Write("");
    MyFile->Close();

//////////////////////////////////////////////////////////////////////////////////////////////////////////

    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");

    //
    TCanvas* MyCanvas01 = new TCanvas;

    Double_t es = (hThetaCD_DATA_bkgnd->GetMaximum())/(hThetaCD_DATA_sign->GetMaximum());
    Double_t scaleMC01 = es*(hThetaCD_MC_qf->GetMaximum())/(hThetaCD_MC_pd->GetMaximum());
    //Double_t scaleMC021 = (hThetaCD_DATA->GetMaximum())/(hThetaCD_MC_pd->GetMaximum());

    double Ymax01 = 1.1*hThetaCD_MC_qf->GetMaximum();

    hThetaCD_MC_qf->SetTitle("");
    hThetaCD_MC_qf->GetXaxis()->SetTitle("#theta_{CD} [deg]");
    hThetaCD_MC_qf->GetXaxis()->SetTitleOffset(1.);
    hThetaCD_MC_qf->GetXaxis()->SetTitleSize(0.06);
    hThetaCD_MC_qf->GetXaxis()->SetLabelSize(0.05);
    hThetaCD_MC_qf->GetXaxis()->SetRangeUser(20.,170.);
    hThetaCD_MC_qf->GetYaxis()->SetTitle("counts");
    hThetaCD_MC_qf->GetYaxis()->SetTitleOffset(1.3);
    hThetaCD_MC_qf->GetYaxis()->SetTitleSize(0.06);
    hThetaCD_MC_qf->GetYaxis()->SetLabelSize(0.05);
    hThetaCD_MC_qf->GetYaxis()->SetRangeUser(0.,Ymax01);

    hThetaCD_MC_qf->SetLineColor(kAzure-2);
    hThetaCD_MC_qf->SetLineStyle(1);
    hThetaCD_MC_qf->SetLineWidth(2);
    //hThetaCD_MC_qf->SetMarkerColor(kAzure-2);
    //hThetaCD_MC_qf->SetMarkerStyle(2);
    //hThetaCD_MC_qf->SetMarkerSize(0.7);
    //hThetaCD_MC_qf->SetFillColor(kAzure-2);
    //hThetaCD_MC_qf->SetFillStyle(3002);
    hThetaCD_MC_qf->Draw("same C");

    hThetaCD_MC_pd->SetLineColor(kSpring-6);
    hThetaCD_MC_pd->SetLineStyle(1);
    hThetaCD_MC_pd->SetLineWidth(2);
    //hThetaCD_MC_pd->SetMarkerColor(kSpring-6);
    //hThetaCD_MC_pd->SetMarkerStyle(3);
    //hThetaCD_MC_pd->SetMarkerSize(0.7);
    //hThetaCD_MC_pd->SetFillColor(kSpring-6);
    //hThetaCD_MC_pd->SetFillStyle(3002);
    hThetaCD_MC_pd->Scale(scaleMC01);
    hThetaCD_MC_pd->Draw("same C");

    TLine* MyLine010 = new TLine(40.,0.,40.,Ymax01);
    MyLine010->SetLineColor(2);
    MyLine010->SetLineStyle(1);
    MyLine010->SetLineWidth(1);
    MyLine010->Draw("same");

    TLine* MyLine011 = new TLine(100.,0.,100.,Ymax01);
    MyLine011->SetLineColor(2);
    MyLine011->SetLineStyle(1);
    MyLine011->SetLineWidth(1);
    MyLine011->Draw("same");

    //legend
    TLegend *MyLegend01= new TLegend(0.575, 0.670, 0.895, 0.885);
    MyLegend01->SetFillStyle(1); MyLegend01->SetFillColor(0); MyLegend01->SetLineColor(0); MyLegend01->SetTextSize(0.04);
    MyLegend01->AddEntry(hThetaCD_MC_qf, "WMC: pd #rightarrow ppn_{spec}", "l");
    MyLegend01->AddEntry(hThetaCD_MC_pd, "WMC: pd #rightarrow pd", "l");
    MyLegend01->AddEntry(MyLine010, "cut", "l");
    MyLegend01->Draw("same");

    MyCanvas01->Print("plots/hThetaCD_MC.png","png");

    //
    TCanvas* MyCanvas02 = new TCanvas;

    Double_t Ymax02 = 1.1*hThetaCD_DATA->GetMaximum();
    Double_t scaleMC020 = (hThetaCD_DATA->GetMaximum())/(hThetaCD_MC_qf->GetMaximum());
    Double_t scaleMC021 = (hThetaCD_DATA->GetMaximum())/(hThetaCD_MC_pd->GetMaximum());

    Double_t scaleExp = (hThetaCD_DATA_exp->GetMaximum())/(hThetaCD_DATA->GetMaximum());

    hThetaCD_DATA->SetTitle("");
    hThetaCD_DATA->GetXaxis()->SetTitle("#theta_{CD} [deg]");
    hThetaCD_DATA->GetXaxis()->SetTitleOffset(1.);
    hThetaCD_DATA->GetXaxis()->SetTitleSize(0.06);
    hThetaCD_DATA->GetXaxis()->SetLabelSize(0.05);
    hThetaCD_DATA->GetXaxis()->SetRangeUser(20.,170.);
    hThetaCD_DATA->GetYaxis()->SetTitle("counts");
    hThetaCD_DATA->GetYaxis()->SetTitleOffset(1.);
    hThetaCD_DATA->GetYaxis()->SetTitleSize(0.06);
    hThetaCD_DATA->GetYaxis()->SetLabelSize(0.05);
    hThetaCD_DATA->GetYaxis()->SetRangeUser(0.,Ymax02);

    hThetaCD_DATA->SetLineColor(0);
    hThetaCD_DATA->SetLineStyle(1);
    hThetaCD_DATA->SetLineWidth(1);
    hThetaCD_DATA->SetMarkerColor(0);
    hThetaCD_DATA->SetMarkerStyle(1);
    hThetaCD_DATA->SetMarkerSize(0.);
    hThetaCD_DATA->DrawCopy("p");

    hThetaCD_DATA_sign_cut->SetLineColor(kOrange+1);
    hThetaCD_DATA_sign_cut->SetLineStyle(1);
    hThetaCD_DATA_sign_cut->SetLineWidth(1);
    //hThetaCD_DATA_sign_cut->SetMarkerColor(kOrange+1);
    //hThetaCD_DATA_sign_cut->SetMarkerStyle(5);
    //hThetaCD_DATA_sign_cut->SetMarkerSize(0.7);
    hThetaCD_DATA_sign_cut->SetFillColor(kOrange+1);
    hThetaCD_DATA_sign_cut->SetFillStyle(3002);
    //gStyle->SetHatchesSpacing(0.7); //to define the spacing between hatches.
    //gStyle->SetHatchesLineWidth(1);
    hThetaCD_DATA_sign_cut->Draw("same C");

    hThetaCD_DATA_bkgnd->SetLineColor(kMagenta+2);
    hThetaCD_DATA_bkgnd->SetLineStyle(1);
    hThetaCD_DATA_bkgnd->SetLineWidth(1);
    //hThetaCD_DATA_bkgnd->SetMarkerColor(kMagenta+2);
    //hThetaCD_DATA_bkgnd->SetMarkerStyle(5);
    //hThetaCD_DATA_bkgnd->SetMarkerSize(0.7);
    hThetaCD_DATA_bkgnd->SetFillColor(kMagenta+2);
    hThetaCD_DATA_bkgnd->SetFillStyle(3002);
    //gStyle->SetHatchesSpacing(0.7); //to define the spacing between hatches.
    //gStyle->SetHatchesLineWidth(1);
    hThetaCD_DATA_bkgnd->Draw("same C");

    hThetaCD_MC->SetLineColor(kPink+3);
    hThetaCD_MC->SetLineStyle(1);
    hThetaCD_MC->SetLineWidth(2);
    hThetaCD_MC->SetMarkerColor(kPink+3);
    hThetaCD_MC->SetMarkerStyle(3);
    hThetaCD_MC->SetMarkerSize(0.7);
    hThetaCD_MC->Scale(scaleMC020);
    //hThetaCD_MC->Draw("same C");

    hThetaCD_MC_qf->SetLineColor(kSpring+4);
    hThetaCD_MC_qf->SetLineStyle(1);
    hThetaCD_MC_qf->SetLineWidth(2);
    hThetaCD_MC_qf->SetMarkerColor(kSpring+4);
    hThetaCD_MC_qf->SetMarkerStyle(23);
    hThetaCD_MC_qf->SetMarkerSize(0.7);
    hThetaCD_MC_qf->Scale(scaleMC020);
    //hThetaCD_MC_qf->Draw("same p");

    hThetaCD_MC_pd->SetLineColor(kYellow+1);
    hThetaCD_MC_pd->SetLineStyle(1);
    hThetaCD_MC_pd->SetLineWidth(2);
    hThetaCD_MC_pd->SetMarkerColor(kYellow+1);
    hThetaCD_MC_pd->SetMarkerStyle(3);
    hThetaCD_MC_pd->SetMarkerSize(0.7);
    hThetaCD_MC_pd->Scale(scaleMC021);
    //hThetaCD_MC_pd->Draw("same C");

    hThetaCD_DATA_sign->SetLineColor(94);
    hThetaCD_DATA_sign->SetLineStyle(1);
    hThetaCD_DATA_sign->SetLineWidth(2);
    //hThetaCD_DATA_sign->SetMarkerColor(94);
    //hThetaCD_DATA_sign->SetMarkerStyle(5);
    //hThetaCD_DATA_sign->SetMarkerSize(0.7);
    hThetaCD_DATA_sign->Draw("same C");

    hThetaCD_DATA_exp->SetLineColor(1);
    hThetaCD_DATA_exp->SetLineStyle(1);
    hThetaCD_DATA_exp->SetLineWidth(1);
    hThetaCD_DATA_exp->SetMarkerColor(1);
    hThetaCD_DATA_exp->SetMarkerStyle(23);
    hThetaCD_DATA_exp->SetMarkerSize(0.7);
    hThetaCD_DATA_exp->Scale(1/scaleExp);
    hThetaCD_DATA_exp->Draw("same E1");

    TLine* MyLine020 = new TLine(40.,0.,40.,Ymax02);
    MyLine020->SetLineColor(2);
    MyLine020->SetLineStyle(1);
    MyLine020->SetLineWidth(1);
    MyLine020->Draw("same");

    TLine* MyLine021 = new TLine(100.,0.,100.,Ymax02);
    MyLine021->SetLineColor(2);
    MyLine021->SetLineStyle(1);
    MyLine021->SetLineWidth(1);
    MyLine021->Draw("same");

    //legend
    TLegend *MyLegend02= new TLegend(0.565, 0.640, 0.895, 0.885);
    MyLegend02->SetFillStyle(1); MyLegend02->SetFillColor(0); MyLegend02->SetLineColor(0); MyLegend02->SetTextSize(0.04);
    MyLegend02->AddEntry(hThetaCD_DATA_exp, "experimental points", "pe");
    MyLegend02->AddEntry(hThetaCD_DATA_sign_cut, "signal", "f");
    MyLegend02->AddEntry(hThetaCD_DATA_bkgnd, "background", "f");
    MyLegend02->AddEntry(MyLine020, "cut", "l");
    MyLegend02->Draw("same");

    MyCanvas02->Print("plots/hThetaCD_DATA.png","png");

    //
    TCanvas* MyCanvas03 = new TCanvas;

    Double_t maxY03 = 1.2*(hExcessEnergy_DATA->GetMaximum());

    hExcessEnergy_DATA->SetTitle("");
    hExcessEnergy_DATA->GetXaxis()->SetTitle("excess energy [MeV]");
    hExcessEnergy_DATA->GetXaxis()->SetTitleOffset(1.);
    hExcessEnergy_DATA->GetXaxis()->SetTitleSize(0.06);
    hExcessEnergy_DATA->GetXaxis()->SetLabelSize(0.05);
    hExcessEnergy_DATA->GetYaxis()->SetTitle("counts");
    hExcessEnergy_DATA->GetYaxis()->SetTitleOffset(0.9);
    hExcessEnergy_DATA->GetYaxis()->SetTitleSize(0.06);
    hExcessEnergy_DATA->GetYaxis()->SetLabelSize(0.05);
    //hExcessEnergy_DATA->GetXaxis()->SetRangeUser(-70.,30.);
    hExcessEnergy_DATA->GetYaxis()->SetRangeUser(0.,maxY03);

    hExcessEnergy_DATA->SetLineColor(1);
    hExcessEnergy_DATA->SetLineStyle(1);
    hExcessEnergy_DATA->SetLineWidth(1);
    hExcessEnergy_DATA->SetMarkerColor(1);
    hExcessEnergy_DATA->SetMarkerStyle(23);
    hExcessEnergy_DATA->SetMarkerSize(1);
    hExcessEnergy_DATA->Draw("E1");

    hExcessEnergy_exp->SetLineColor(kGreen+2);
    hExcessEnergy_exp->SetLineStyle(1);
    hExcessEnergy_exp->SetLineWidth(1);
    hExcessEnergy_exp->SetMarkerColor(kGreen+2);
    hExcessEnergy_exp->SetMarkerStyle(4);
    hExcessEnergy_exp->SetMarkerSize(1);
    //hExcessEnergy_exp->Draw("same p");

    hExcessEnergy_sign->SetLineColor(kOrange+1);
    hExcessEnergy_sign->SetLineStyle(1);
    hExcessEnergy_sign->SetLineWidth(1);
    hExcessEnergy_sign->SetMarkerColor(kOrange+1);
    hExcessEnergy_sign->SetMarkerStyle(24);
    hExcessEnergy_sign->SetMarkerSize(0.8);
    hExcessEnergy_sign->Draw("same E1");

    hExcessEnergy_bkgnd->SetLineColor(kMagenta+2);
    hExcessEnergy_bkgnd->SetLineStyle(1);
    hExcessEnergy_bkgnd->SetLineWidth(1);
    hExcessEnergy_bkgnd->SetMarkerColor(kMagenta+2);
    hExcessEnergy_bkgnd->SetMarkerStyle(21);
    hExcessEnergy_bkgnd->SetMarkerSize(0.8);
    hExcessEnergy_bkgnd->Draw("same E1");

    hExcessEnergy_diff->SetLineColor(kGreen);
    hExcessEnergy_diff->SetLineStyle(1);
    hExcessEnergy_diff->SetLineWidth(1);
    hExcessEnergy_diff->SetMarkerColor(kGreen);
    hExcessEnergy_diff->SetMarkerStyle(4);
    hExcessEnergy_diff->SetMarkerSize(1);
    //hExcessEnergy_diff->Draw("same p");

    //legend
    TLegend *MyLegend03= new TLegend(0.550, 0.430, 0.890, 0.620);
    MyLegend03->SetFillStyle(1); MyLegend03->SetFillColor(0); MyLegend03->SetLineColor(0); MyLegend03->SetTextSize(0.04);
    MyLegend03->AddEntry(hExcessEnergy_DATA, "experimental points", "pe");
    MyLegend03->AddEntry(hExcessEnergy_sign, "signal", "pe");
    MyLegend03->AddEntry(hExcessEnergy_bkgnd, "background pd #rightarrow pd", "pe");
    MyLegend03->Draw("same");

    MyCanvas03->Print("plots/hExcessEnergy.png","png");

}


