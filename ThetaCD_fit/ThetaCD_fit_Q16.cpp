/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2017-07
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-05
***********************************************/

//Macro for fitting ThetaCD distributions

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

/////////////////////////////////////////////////////////////////////////

//Gauss function
Double_t Gauss_sign(Double_t *x, Double_t *p) {
    Double_t out_gauss = 0.;
    out_gauss = p[0]*(exp(-0.5*(x[0] - p[1])*(x[0] - p[1])/(p[2]*p[2])));
    return out_gauss;
}

//Total function = Gauss_1 + Gauss_2
Double_t Gauss_sign_bkgnd(Double_t *y, Double_t *r) {
    Double_t out_total = 0.;
    Double_t Gauss_1 = 0.;
    Double_t Gauss_2 = 0.;
    Gauss_1 = r[0]*(exp(-0.5*(y[0] - r[1])*(y[0] - r[1])/(r[2]*r[2])));
    Gauss_2 = r[3]*(exp(-0.5*(y[0] - r[4])*(y[0] - r[4])/(r[5]*r[5])));
    out_total = Gauss_1 + Gauss_2;
    return out_total;
}

/////////////////////////////////////////////////////////////////////////

void ThetaCD_fit_Q16() {

    Int_t q = 16;

    cout<<fixed;
    cout.precision(2);

    //number of parameters//
    sign_param = 3;
    bkgnd_param = 3;
    total_param = sign_param + bkgnd_param;

    //draw options//
    begin_draw_sign = 40.;
    end_draw_sign = 100.;
    begin_draw_bkgnd = 40.;
    end_draw_bkgnd = 100.;
    begin_draw_total = 40.;
    end_draw_total = 100.;

    //open the ROOT file with data
    f[0] = new TFile("../input/DATA_ppn_qf_offset.root","READ");
    f[1] = new TFile("../input/MC-ppn_qf-PARIS-x6.root","READ");
    //f[1] = new TFile("../input/MC-ppn_qf-CDBONN-x6.root","READ");

    //create new file.root
    MyFile = new TFile(Form("files/ThetaCD_Q%d.root",q),"RECREATE");

    f[0]->cd("Histograms/DATA_lev4/ThetaCD_lev4/ThetaCD_Q");
    hThetaCD_DATA_Q = (TH1D*)gDirectory->Get(Form("hThetaCD_lev4_Q%d",q));
    //hThetaCD_DATA_Q->Rebin(2);

    Double_t bin_content = 0.;

    //main loop
    for (Int_t i = 1; i < 16; i++) {

        f[0]->cd("Histograms/DATA_lev5/ThetaCD_lev5");

        hThetaCD_DATA_i[i] = (TH1D*)gDirectory->Get(Form("ThetaCD_Q%d/hThetaCD_lev5_Q%d_Th%d",q,q,i));
        hThetaCD_DATA_i[i]->Rebin(2);

        hThetaCD_DATA_sign_i[i] = new TH1D(Form("hThetaCD_DATA_sign_Q%d_Th%d",q,i),"",180,0.,180.);
        hThetaCD_DATA_bkgnd_i[i] = new TH1D(Form("hThetaCD_DATA_bkgnd_Q%d_Th%d",q,i),"",180,0.,180.);

        hThetaCD_DATA_sign_i[i] = (TH1D*)hThetaCD_DATA_i[i]->Clone(Form("hThetaCD_DATA_sign_Q%d_Th%d",q,i));
        //hThetaCD_DATA_bkgnd_i[i] = (TH1D*)hThetaCD_DATA_i[i]->Clone(Form("hThetaCD_DATA_bkgnd_Q%d_Th%d",q,i));

        for(Int_t j = 1; j < 181; j++) {
            bin_content = hThetaCD_DATA_i[i]->GetBinContent(j);

            if(i >= 12) {
                hThetaCD_DATA_bkgnd_i[i]->SetBinContent(j,bin_content);
            }
            else {
                hThetaCD_DATA_bkgnd_i[i]->SetBinContent(j,0.);
            }
        }

        if (i >= 12) {

            if (i == 12) {
                begin_fit_sign = 62.;
                end_fit_sign = 75.;
                begin_fit_bkgnd = 74.;
                end_fit_bkgnd = 79.;
                begin_fit_total = 62.;
                end_fit_total = 79.;
            }

            if (i == 13) {
                begin_fit_sign = 59.;
                end_fit_sign = 71.5;
                begin_fit_bkgnd = 72.5;
                end_fit_bkgnd = 77.;
                begin_fit_total = 59.;
                end_fit_total = 77.;
            }

            if (i == 14) {
                begin_fit_sign = 60.;
                end_fit_sign = 71.;
                begin_fit_bkgnd = 71.5;
                end_fit_bkgnd = 77.;
                begin_fit_total = 60.;
                end_fit_total = 77.;
            }

            if (i == 15) {
                begin_fit_sign = 59.;
                end_fit_sign = 70.;
                begin_fit_bkgnd = 71.5;
                end_fit_bkgnd = 75.;
                begin_fit_total = 59.;
                end_fit_total = 75.;
            }

            ////SIGNAL FIT////
            signalFunct[i] = new TF1(Form("signalFunct_%d",i),Gauss_sign,begin_draw_sign,end_draw_sign,sign_param);

            //set some parameters to start
            signalFunct[i]->SetParameter(0,900.);
            signalFunct[i]->SetParameter(1,70.1);
            signalFunct[i]->SetParameter(2,7.5);

            hThetaCD_DATA_i[i]->Fit(signalFunct[i],"","",begin_fit_sign,end_fit_sign);

            ////BACKGROUND FIT////
            backgroundFunct[i] = new TF1(Form("backgroundFunct_%d",1),Gauss_sign,begin_draw_bkgnd,end_draw_bkgnd,bkgnd_param);

            //set some parameters to start
            backgroundFunct[i]->SetParameter(0,850.);
            backgroundFunct[i]->SetParameter(1,75.);
            backgroundFunct[i]->SetParameter(2,3.);

            hThetaCD_DATA_i[i]->Fit(backgroundFunct[i],"","",begin_fit_bkgnd,end_fit_bkgnd);

            ////TOTAL////
            totalFunct[i] = new TF1(Form("totalFunct_%d",i),Gauss_sign_bkgnd,begin_draw_total,end_draw_total,total_param);

            //take parameters from previous fits
            for (Int_t k = 0; k < sign_param; k++) {
                totalFunct[i]->SetParameter(k,signalFunct[i]->GetParameter(k));
            }

            for(Int_t k = sign_param; k < total_param; k++) {
                totalFunct[i]->SetParameter(k,backgroundFunct[i]->GetParameter(k-sign_param));
            }

            //FIT
            hThetaCD_DATA_i[i]->Fit(totalFunct[i],"","",begin_fit_total,end_fit_total);

            //set background for whole spectrum taking parameters from total function
            signalFunct_end[i] = new TF1(Form("signalFunct_end_%d",i),Gauss_sign,begin_draw_sign,end_draw_sign,sign_param);
            backgroundFunct_end[i] = new TF1(Form("backgroundFunct_end_%d",i),Gauss_sign,begin_draw_bkgnd,end_draw_bkgnd,bkgnd_param);

            for(Int_t k = 0; k < sign_param; k++) {
                signalFunct_end[i]->SetParameter(k,totalFunct[i]->GetParameter(k));
            }

            for(Int_t k = 0; k < bkgnd_param; k++) {
                backgroundFunct_end[i]->SetParameter(k,totalFunct[i]->GetParameter(k+sign_param));
            }

            hThetaCD_DATA_sign_i[i]->Add(backgroundFunct_end[i],-1);
            //hThetaCD_DATA_bkgnd_i[i]->Add(signalFunct_end[i],-1);
            hThetaCD_DATA_bkgnd_i[i]->Add(hThetaCD_DATA_sign_i[i],-1);

        }

        Integral_DATA = hThetaCD_DATA_i[i]->Integral();
        Integral_signal = hThetaCD_DATA_sign_i[i]->Integral();
        Integral_bkgnd = hThetaCD_DATA_bkgnd_i[i]->Integral();

        //cout<<"DATA = "<<Integral_DATA<<endl;
        //cout<<"signal = "<<Integral_signal<<endl;
        //cout<<"background = "<<Integral_bkgnd<<endl;

        //all ThetaFD bins together
        hThetaCD_DATA->Add(hThetaCD_DATA_i[i]);
        hThetaCD_DATA_sign->Add(hThetaCD_DATA_sign_i[i]);
        hThetaCD_DATA_bkgnd->Add(hThetaCD_DATA_bkgnd_i[i]);

        //Monte Carlo simulations
        f[1]->cd("Histograms/DATA_lev4/ThetaCD_lev4");

        hThetaCD_MC_i[i] = (TH1D*)gDirectory->Get(Form("ThetaCD_Q%d/hThetaCD_lev4_Q%d_Th%d",q,q,i));
        hThetaCD_MC_i[i]->Rebin(2);

        hThetaCD_MC->Add(hThetaCD_MC_i[i]);

        //new file with results
        ofstream results;
        results.open(Form("files/ThetaCD_Q%d.dat",q), ios::app);
        results<<Form("%i %g %g %g %g",i,Integral_DATA,Integral_signal,Integral_bkgnd,(Integral_bkgnd/Integral_DATA)*100)<<endl;
        results.close();

    }

    cout<<"background for Q ("<<-70.+(q-1)*2.5<<","<<-67.5+(q-1)*2.5<<") MeV is equal to "<<(hThetaCD_DATA_bkgnd->Integral()/hThetaCD_DATA->Integral())*100<<"%"<<endl;

    //create new output file.root
    MyFile->cd();
    hThetaCD_MC->Write(Form("hThetaCD_MC_Q%d",q));
    hThetaCD_DATA->Write(Form("hThetaCD_DATA_Q%d",q));
    hThetaCD_DATA_sign->Write(Form("hThetaCD_DATA_sign_Q%d",q));
    hThetaCD_DATA_bkgnd->Write(Form("hThetaCD_DATA_bckr_Q%d",q));
    for (Int_t k = 1; k < 16; k++) {
        hThetaCD_MC_i[k]->Write(Form("hThetaCD_MC_Q%d_Th%d",q,k));
        hThetaCD_DATA_i[k]->Write(Form("hThetaCD_DATA_Q%d_Th%d",q,k));
        hThetaCD_DATA_sign_i[k]->Write(Form("hThetaCD_DATA_sign_Q%d_Th%d",q,k));
        hThetaCD_DATA_bkgnd_i[k]->Write(Form("hThetaCD_DATA_bkgnd_Q%d_Th%d",q,k));
    }
    MyFile->Close();

//////////////////////////////////////////////////////////////////////////////////////////////////////////

    //set no statistics on histograms and palette
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.10);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");

    //draw histograms on separated canvas

    TCanvas* MyCanvas01 = new TCanvas; //create Canvas

    for (Int_t j = 12; j < 16; j++) {

        hThetaCD_MC_i[j]->Rebin(2);

        Double_t Ymax01 = 1.1*hThetaCD_DATA_i[j]->GetMaximum();
        Double_t scaleMC01 = (hThetaCD_DATA_i[j]->GetMaximum())/(hThetaCD_MC_i[j]->GetMaximum());

        TString descr01 = Form("Q #in (%G,%G) MeV, #theta_{FD} #in (%G#circ,%G#circ)",-70.+(q-1)*2.5,-67.5+(q-1)*2.5,3.+(j-1),4.+(j-1));

        hThetaCD_DATA_i[j]->SetTitle(descr01);
        hThetaCD_DATA_i[j]->GetXaxis()->SetTitle("#theta_{CD},#circ");
        hThetaCD_DATA_i[j]->GetXaxis()->SetTitleOffset(1.);
        hThetaCD_DATA_i[j]->GetXaxis()->SetTitleSize(0.06);
        hThetaCD_DATA_i[j]->GetXaxis()->SetLabelSize(0.05);
        hThetaCD_DATA_i[j]->GetXaxis()->SetRangeUser(20.,170.);
        hThetaCD_DATA_i[j]->GetYaxis()->SetTitle("counts");
        hThetaCD_DATA_i[j]->GetYaxis()->SetTitleOffset(1.);
        hThetaCD_DATA_i[j]->GetYaxis()->SetTitleSize(0.06);
        hThetaCD_DATA_i[j]->GetYaxis()->SetLabelSize(0.05);
        hThetaCD_DATA_i[j]->GetYaxis()->SetRangeUser(0.,Ymax01);

        hThetaCD_DATA_i[j]->SetLineColor(kBlack);
        hThetaCD_DATA_i[j]->SetLineStyle(1);
        hThetaCD_DATA_i[j]->SetLineWidth(1);
        hThetaCD_DATA_i[j]->SetMarkerColor(kBlack);
        hThetaCD_DATA_i[j]->SetMarkerStyle(23);
        hThetaCD_DATA_i[j]->SetMarkerSize(0.7);
        hThetaCD_DATA_i[j]->Draw("pE1");

        hThetaCD_MC_i[j]->SetLineColor(kAzure-3);
        hThetaCD_MC_i[j]->SetLineStyle(1);
        hThetaCD_MC_i[j]->SetLineWidth(1);
        hThetaCD_MC_i[j]->SetMarkerColor(kAzure-3);
        hThetaCD_MC_i[j]->SetMarkerStyle(3);
        hThetaCD_MC_i[j]->SetMarkerSize(0.7);
        hThetaCD_MC_i[j]->Scale(scaleMC01);
        //hThetaCD_MC_i[j]->Draw("same E1");

        hThetaCD_DATA_sign_i[j]->SetLineColor(kOrange+1);
        hThetaCD_DATA_sign_i[j]->SetLineStyle(1);
        hThetaCD_DATA_sign_i[j]->SetLineWidth(2);
        hThetaCD_DATA_sign_i[j]->SetMarkerColor(kOrange+1);
        hThetaCD_DATA_sign_i[j]->SetMarkerStyle(5);
        hThetaCD_DATA_sign_i[j]->SetMarkerSize(0.7);
        hThetaCD_DATA_sign_i[j]->Draw("same C");

        hThetaCD_DATA_bkgnd_i[j]->SetLineColor(kMagenta+2);
        hThetaCD_DATA_bkgnd_i[j]->SetLineStyle(1);
        hThetaCD_DATA_bkgnd_i[j]->SetLineWidth(2);
        hThetaCD_DATA_bkgnd_i[j]->SetMarkerColor(kMagenta+2);
        hThetaCD_DATA_bkgnd_i[j]->SetMarkerStyle(5);
        hThetaCD_DATA_bkgnd_i[j]->SetMarkerSize(0.7);
        hThetaCD_DATA_bkgnd_i[j]->Draw("same C");

        totalFunct[j]->SetLineColor(kCyan-3);
        totalFunct[j]->SetLineStyle(1);
        totalFunct[j]->SetLineWidth(2);
        totalFunct[j]->Draw("same C");

        signalFunct[j]->SetLineColor(kTeal+4);
        signalFunct[j]->SetLineStyle(1);
        signalFunct[j]->SetLineWidth(2);
        //signalFunct[j]->Draw("same C");

        signalFunct_end[j]->SetLineColor(kRed+2);
        signalFunct_end[j]->SetLineStyle(1);
        signalFunct_end[j]->SetLineWidth(2);
        //signalFunct_end[j]->Draw("same C");

        backgroundFunct[j]->SetLineColor(kSpring+4);
        backgroundFunct[j]->SetLineStyle(1);
        backgroundFunct[j]->SetLineWidth(2);
        //backgroundFunct[j]->Draw("same C");

        backgroundFunct_end[j]->SetLineColor(kBlue+1);
        backgroundFunct_end[j]->SetLineStyle(1);
        backgroundFunct_end[j]->SetLineWidth(2);
        //backgroundFunct_end[j]->Draw("same C");

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
        TLegend *MyLegend01= new TLegend(0.540, 0.560, 0.880, 0.885);
        MyLegend01->SetFillStyle(0); MyLegend01->SetFillColor(0); MyLegend01->SetLineColor(0); MyLegend01->SetTextSize(0.04);
        MyLegend01->AddEntry(hThetaCD_DATA_i[j], "experimental points", "ep");
        MyLegend01->AddEntry(hThetaCD_DATA_sign_i[j], "signal", "l");
        MyLegend01->AddEntry(totalFunct[j], "fitting function", "l");
        MyLegend01->AddEntry(hThetaCD_DATA_bkgnd_i[j], "background: pd #rightarrow pd", "l");
        //MyLegend01->AddEntry(hThetaCD_MC_i[j], "WMC: pd #rightarrow ppn_{spectator}", "ep");
        MyLegend01->AddEntry(MyLine010, "cut: #theta_{CD}#in(40,100)#circ", "l");
        MyLegend01->Draw("same");

        MyCanvas01->Print(Form("plots/hThetaCD_DATA_Q%d_Th%d.png",q,j),"png");

    }

    TCanvas* MyCanvas02 = new TCanvas;

    Double_t Ymax02 = 1.1*hThetaCD_DATA->GetMaximum();
    Double_t scaleMC020 = (hThetaCD_DATA->GetMaximum())/(hThetaCD_MC->GetMaximum());

    TString descr02 = Form("Q #in (%G,%G) MeV",-70.+(q-1)*2.5,-67.5+(q-1)*2.5);

    hThetaCD_DATA->SetTitle(descr02);
    hThetaCD_DATA->GetXaxis()->SetTitle("#theta_{CD},#circ");
    hThetaCD_DATA->GetXaxis()->SetTitleOffset(1.);
    hThetaCD_DATA->GetXaxis()->SetTitleSize(0.06);
    hThetaCD_DATA->GetXaxis()->SetLabelSize(0.05);
    hThetaCD_DATA->GetXaxis()->SetRangeUser(20.,170.);
    hThetaCD_DATA->GetYaxis()->SetTitle("counts");
    hThetaCD_DATA->GetYaxis()->SetTitleOffset(1.);
    hThetaCD_DATA->GetYaxis()->SetTitleSize(0.06);
    hThetaCD_DATA->GetYaxis()->SetLabelSize(0.05);
    hThetaCD_DATA->GetYaxis()->SetRangeUser(0.,Ymax02);

    hThetaCD_DATA->SetLineColor(kBlack);
    hThetaCD_DATA->SetLineStyle(1);
    hThetaCD_DATA->SetLineWidth(1);
    hThetaCD_DATA->SetMarkerColor(kBlack);
    hThetaCD_DATA->SetMarkerStyle(1);
    hThetaCD_DATA->SetMarkerSize(0.);
    hThetaCD_DATA->Draw("p");

    hThetaCD_MC->SetLineColor(kAzure-3);
    hThetaCD_MC->SetLineStyle(1);
    hThetaCD_MC->SetLineWidth(1);
    hThetaCD_MC->SetMarkerColor(kAzure-3);
    hThetaCD_MC->SetMarkerStyle(3);
    hThetaCD_MC->SetMarkerSize(0.7);
    hThetaCD_MC->Scale(scaleMC020);
    hThetaCD_MC->Draw("same E1");

    hThetaCD_DATA_Q->SetLineColor(kBlack);
    hThetaCD_DATA_Q->SetLineStyle(1);
    hThetaCD_DATA_Q->SetLineWidth(1);
    hThetaCD_DATA_Q->SetMarkerColor(kBlack);
    hThetaCD_DATA_Q->SetMarkerStyle(23);
    hThetaCD_DATA_Q->SetMarkerSize(0.7);
    hThetaCD_DATA_Q->Scale(2);
    hThetaCD_DATA_Q->Draw("same E1");

    hThetaCD_DATA_sign->SetLineColor(kOrange+1);
    hThetaCD_DATA_sign->SetLineStyle(1);
    hThetaCD_DATA_sign->SetLineWidth(2);
    hThetaCD_DATA_sign->SetMarkerColor(kOrange+1);
    hThetaCD_DATA_sign->SetMarkerStyle(5);
    hThetaCD_DATA_sign->SetMarkerSize(0.7);
    hThetaCD_DATA_sign->Draw("same C");

    hThetaCD_DATA_bkgnd->SetLineColor(kMagenta+2);
    hThetaCD_DATA_bkgnd->SetLineStyle(1);
    hThetaCD_DATA_bkgnd->SetLineWidth(2);
    hThetaCD_DATA_bkgnd->SetMarkerColor(kMagenta+2);
    hThetaCD_DATA_bkgnd->SetMarkerStyle(5);
    hThetaCD_DATA_bkgnd->SetMarkerSize(0.7);
    hThetaCD_DATA_bkgnd->Draw("same C");

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
    TLegend *MyLegend02= new TLegend(0.540, 0.600, 0.880, 0.885);
    MyLegend02->SetFillStyle(0); MyLegend02->SetFillColor(0); MyLegend02->SetLineColor(0); MyLegend02->SetTextSize(0.04);
    MyLegend02->AddEntry(hThetaCD_DATA_Q, "experimental points", "pe");
    MyLegend02->AddEntry(hThetaCD_DATA_sign, "signal", "l");
    MyLegend02->AddEntry(hThetaCD_DATA_bkgnd, "background: pd #rightarrow pd", "l");
    MyLegend02->AddEntry(hThetaCD_MC, "WMC: pd #rightarrow ppn_{spectator}", "ep");
    MyLegend02->AddEntry(MyLine020, "cut: #theta_{CD}#in(40,100)#circ", "l");
    MyLegend02->Draw("same");

    MyCanvas02->Print(Form("plots/hThetaCD_DATA_Q%d.png",q),"png");

}
