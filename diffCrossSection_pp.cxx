/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2016-07
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-06
***********************************************/

//Comparison of the differential cross sections for the proton-proton elastic scattering from the EDDA Collaboration measurents and SAID programe calculations

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
#include "TString.h"
#include <iostream>
#include <fstream>

void diffCrossSection_pp() {

    //DATA from EDDA measurements

    TGraphErrors* XS_graph[181];

    //Theta_scatt^{*} = 35 deg
    double x_35[] = {1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125,1.8375,1.8625,1.8875,1.9125,1.9375,1.9625,1.9875,2.0125};
    double y_35[] = {7.235,7.212,7.171,7.235,7.124,7.089,7.003,7.028,6.983,6.946,7.002,6.729,6.772,6.667,6.573,6.496,6.404,6.343,6.162};
    //total error (stat+syst)
    double y_err_35[] = {0.206,0.205,0.201,0.203,0.2,0.199,0.196,0.197,0.196,0.195,0.196,0.189,0.190,0.189,0.187,0.183,0.180,0.179,0.174};
    XS_graph[35] = new TGraphErrors(20,x_35,y_35,0,y_err_35);

    XS_graph[35]->SetMarkerStyle(23);
    XS_graph[35]->SetMarkerSize(0.8);
    XS_graph[35]->SetMarkerColor(1);

    //Theta_scatt^{*} = 41 deg
    double x_41[] = {0.9875,1.0125,1.0375,1.0625,1.0875,1.1125,1.1375,1.1625,1.1875,1.2125,1.2375,1.2625,1.2875,1.3125,1.3375,1.3625,1.3875,1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};
    double y_41[] = {4.066,4.190,4.233,4.364,4.478,4.601,4.747,4.881,4.939,5.161,5.195,5.305,5.374,5.518,5.523,5.580,5.626,5.600,5.573,5.500,5.566,5.533,5.460,5.369,5.237,5.193,5.141,5.016,4.952,4.852,4.828,4.732,4.671};
    //total error (stat+syst)
    double y_err_41[] = {0.117,0.121,0.122,0.125,0.128,0.128,0.132,0.136,0.137,0.143,0.144,0.147,0.149,0.153,0.153,0.155,0.156,0.155,0.155,0.153,0.155,0.154,0.153,0.150,0.146,0.145,0.143,0.141,0.139,0.136,0.135,0.132,0.131};
    XS_graph[41] = new TGraphErrors(33,x_41,y_41,0,y_err_41);

    XS_graph[41]->SetMarkerStyle(23);
    XS_graph[41]->SetMarkerSize(0.8);
    XS_graph[41]->SetMarkerColor(1);

    //Theta_scatt^{*} = 55 deg
    double x_55[] = {0.7375,0.7625,0.7875,0.8125,0.8375,0.8625,0.8875,0.9125,0.9375,0.9625,0.9875,1.0125,1.0375,1.0625,1.0875,1.1125,1.1375,1.1625,1.1875,1.2125,1.2375,1.2625,1.2875,1.3125,1.3375,1.3625,1.3875,1.4125,1.4375,1462.5,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125,1.8375,1.8625,1.8875,1.9125,1.9375,1.9625,1.9875,2.0125};
    double y_55[] = {3.772,3.835,3.870,3.824,3.861,3.921,3.952,3.905,3.900,3.909,3.908,3.964,3.940,3.948,3.987,3.941,3.978,3.927,3.866,3.896,3.826,3.720,3.621,3.553,3.513,3.385,3.372,3.229,3.063,3.055,2.953,2.855,2.756,2.652,2.576,2.552,2.422,2.358,2.274,2.205,2.148,2.073,2.012,1.946,1.858,1.857,1.712,1.677,1.637,1.557,1.486,1.442};
    //total error (stat+syst)
    double y_err_55[] = {0.106,0.104,0.109,0.108,0.109,0.111,0.112,0.111,0.111,0.111,0.111,0.114,0.113,0.112,0.113,0.112,0.111,0.109,0.108,0.109,0.107,0.104,0.101,0.099,0.098,0.095,0.094,0.091,0.086,0.086,0.084,0.081,0.078,0.076,0.074,0.073,0.069,0.068,0.066,0.064,0.062,0.06,0.058,0.057,0.054,0.054,0.052,0.051,0.048,0.046,0.044,0.043};
    XS_graph[55] = new TGraphErrors(53,x_55,y_55,0,y_err_55);

    XS_graph[55]->SetMarkerStyle(23);
    XS_graph[55]->SetMarkerSize(0.8);
    XS_graph[55]->SetMarkerColor(1);

    //Theta_scatt^{*} = 75 deg
    double x_75[] = {0.7125,0.7375,0.7625,0.7875,0.8125,0.8375,0.8625,0.8875,0.9125,0.9375,0.9625,0.9875,1.0125,1.0375,1.0625,1.0875,1.1125,1.1375,1.1625,1.1875,1.2125,1.2375,1.2625,1.2875,1.3125,1.3375,1.3625,1.3875,1.4125,1.4375,1462.5,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125,1.8375,1.8625,1.8875,1.9125,1.9375,1.9625,1.9875,2.0125};
    double y_75[] = {3.703,3.765,3.859,3.821,3.844,3.866,3.888,3.872,3.842,3.829,3.725,3.718,3.563,3.543,3.476,3.386,3.219,3.270,3.111,3.006,2.813,2.683,2.505,2.385,2.151,2.042,1.877,1.735,1.621,1.486,1.437,1.304,1.239,1.122,1.097,1.043,0.967,0.912,0.884,0.838,0.809,0.778,0.752,0.723,0.700,0.674,0.660,0.613,0.617,0.594,0.571,0.545,0.540};
    //total error (stat+syst)
    double y_err_75[] = {0.104,0.106,0.108,0.107,0.108,0.109,0.109,0.109,0.108,0.108,0.105,0.105,0.102,0.102,0.099,0.097,0.092,0.091,0.087,0.084,0.079,0.075,0.071,0.067,0.061,0.058,0.054,0.050,0.047,0.043,0.042,0.039,0.037,0.033,0.033,0.032,0.029,0.028,0.027,0.026,0.025,0.024,0.023,0.022,0.022,0.021,0.021,0.02,0.02,0.019,0.018,0.018,0.018};
    XS_graph[75] = new TGraphErrors(54,x_75,y_75,0,y_err_75);

    XS_graph[75]->SetMarkerStyle(23);
    XS_graph[75]->SetMarkerSize(0.8);
    XS_graph[75]->SetMarkerColor(1);

    //Theta_scatt^{*} = 89 deg
    double x_89[] = {0.7125,0.7375,0.7625,0.7875,0.8125,0.8375,0.8625,0.8875,0.9125,0.9375,0.9625,0.9875,1.0125,1.0375,1.0625,1.0875,1.1125,1.1375,1.1625,1.1875,1.2125,1.2375,1.2625,1.2875,1.3125,1.3375,1.3625,1.3875,1.4125,1.4375,1462.5,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125,1.8375,1.8625,1.8875,1.9125,1.9375,1.9625,1.9875,2.0125};
    double y_89[] = {3.886,3.822,3.840,3.836,3.782,3.771,3.829,3.773,3.797,3.742,3.774,3.754,3.751,3.667,3.640,3.506,3.353,3.320,3.159,3.019,2.768,2.575,2.306,2.118,1.924,1.721,1.538,1.402,1.293,1.169,1.038,0.977,0.914,0.841,0.776,0.740,0.695,0.675,0.623,0.597,0.601,0.572,0.553,0.537,0.526,0.507,0.489,0.491,0.485,0.454,0.453,0.437,0.440};
    //total error (stat+syst)
    double y_err_89[] = {0.109,0.107,0.108,0.108,0.106,0.106,0.108,0.106,0.107,0.106,0.107,0.106,0.107,0.105,0.103,0.100,0.096,0.093,0.088,0.084,0.078,0.072,0.065,0.060,0.055,0.049,0.044,0.041,0.038,0.035,0.031,0.030,0.028,0.026,0.025,0.024,0.022,0.021,0.020,0.019,0.019,0.018,0.018,0.017,0.017,0.017,0.016,0.017,0.017,0.015,0.015,0.015,0.015};
    XS_graph[89] = new TGraphErrors(54,x_89,y_89,0,y_err_89);

    XS_graph[89]->SetMarkerStyle(23);
    XS_graph[89]->SetMarkerSize(0.8);
    XS_graph[89]->SetMarkerColor(1);


///////////////////////////////////////////////////

    //differential cross sections from SAID for Theta_scatt^{*} = (0,180) deg and p_pbeam = (0,2) GeV/c (beam momentum range taken from the histogram (simulation result))

    TGraph *gXS[181];
    FILE *said_file;

    for (int th = 0; th < 181; th++) {

        said_file = fopen(Form("input/SAID_data/Theta_%i.txt", th), "r");

        gXS[th] = new TGraph();

        float p_pbeam = 0.;
        float XS = 0.;
        int i = -1;

        while (!feof(said_file)) {
            i++;
            fscanf(said_file, "%f %f\n", &p_pbeam, &XS);
            p_pbeam = p_pbeam/1000.;
            gXS[th]->SetPoint(i,p_pbeam,XS);
        }
    }


///////////////////////////////////////////////////

    //effective beam momentum (simulation result)

    TH1F *hBeam_momentum_eff = new TH1F("hBeam_momentum_eff","",100,0.,3.);

    ifstream Mom_file;
    Mom_file.open("input/ProtonVariables-PARIS-1.dat");
    //Mom_file.open("input/ProtonVariables-CDBONN-1.dat");

    double p_beam_eff;
    double Theta_scatt_cm;

    for (int j = 1; j < 1000001; j++) {

        Mom_file>>p_beam_eff>>Theta_scatt_cm;

        hBeam_momentum_eff->Fill(p_beam_eff);

    }

///////////////////////////////////////////////////

    //set no statistics on histograms and pallette
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.15);


    TCanvas* MyCanvas01 = new TCanvas;

    //hBeam_momentum_eff->SetTitle("Effective beam momentum");
    hBeam_momentum_eff->GetXaxis()->SetTitle("effective beam momentum, GeV/c");
    hBeam_momentum_eff->GetYaxis()->SetTitle("d#sigma/d#Omega, mb/sr");
    hBeam_momentum_eff->GetXaxis()->SetTitleSize(0.06);
    hBeam_momentum_eff->GetXaxis()->SetTitleOffset(1.1);
    hBeam_momentum_eff->GetXaxis()->SetLabelSize(0.05);
    hBeam_momentum_eff->GetYaxis()->SetTitleSize(0.06);
    hBeam_momentum_eff->GetYaxis()->SetTitleOffset(0.7);
    hBeam_momentum_eff->GetYaxis()->SetLabelSize(0.05);
    hBeam_momentum_eff->GetXaxis()->SetRangeUser(0.7,2.5);
    hBeam_momentum_eff->GetYaxis()->SetRangeUser(0.0,20.0);

    hBeam_momentum_eff->SetLineWidth(2);
    hBeam_momentum_eff->SetLineColor(kPink);
    hBeam_momentum_eff->SetLineStyle(2);
    hBeam_momentum_eff->Scale(0.00012);
    hBeam_momentum_eff->DrawCopy("C");

    gXS[35]->SetLineColor(51);
    gXS[35]->SetLineWidth(1);
    gXS[35]->Draw("same");

    gXS[41]->SetLineColor(60);
    gXS[41]->SetLineWidth(1);
    gXS[41]->Draw("same");

    gXS[55]->SetLineColor(66);
    gXS[55]->SetLineWidth(1);
    gXS[55]->Draw("same");

    gXS[75]->SetLineColor(75);
    gXS[75]->SetLineWidth(1);
    gXS[75]->Draw("same");

    gXS[89]->SetLineColor(92);
    gXS[89]->SetLineWidth(1);
    gXS[89]->Draw("same");

    gXS[150]->SetLineColor(1);
    gXS[150]->SetLineWidth(1);
    gXS[150]->Draw("same");

    XS_graph[35]->Draw("same p");
    XS_graph[41]->Draw("same p");
    XS_graph[55]->Draw("same p");
    XS_graph[75]->Draw("same p");
    XS_graph[89]->Draw("same p");

    TLegend *MyLegend01 = new TLegend(0.145, 0.475, 0.345, 0.885);
    MyLegend01->SetFillStyle(1); MyLegend01->SetFillColor(0); MyLegend01->SetLineColor(0); MyLegend01->SetTextSize(0.04);
    MyLegend01->SetBorderSize(5);
    MyLegend01->AddEntry( XS_graph[35], "EDDA" , "p");
    MyLegend01->AddEntry((TObject*)0, "SAID:", "");
    MyLegend01->AddEntry( gXS[35], "#theta^{cm} = 35#circ" , "l");
    MyLegend01->AddEntry( gXS[41], "#theta^{cm} = 41#circ" , "l");
    MyLegend01->AddEntry( gXS[55], "#theta^{cm} = 55#circ" , "l");
    MyLegend01->AddEntry( gXS[75], "#theta^{cm} = 75#circ" , "l");
    MyLegend01->AddEntry( gXS[89], "#theta^{cm} = 89#circ" , "l");
    MyLegend01->AddEntry( gXS[150], "#theta^{cm} = 150#circ" , "l");
    MyLegend01->Draw("same");

    TLegend *MyLegend01_pb = new TLegend(0.600, 0.820, 0.885, 0.885);
    MyLegend01_pb->SetFillStyle(1); MyLegend01_pb->SetFillColor(0); MyLegend01_pb->SetLineColor(0); MyLegend01_pb->SetTextSize(0.04);
    MyLegend01_pb->AddEntry( hBeam_momentum_eff, "beam momentum" , "l");
    MyLegend01_pb->Draw("same");

    MyCanvas01->Print("output/plots/hCrossSections.png","png");
    MyCanvas01->Print("output/plots/hCrossSections.eps","eps");

    //
    TCanvas* MyCanvas01a = new TCanvas;

    hBeam_momentum_eff->GetXaxis()->SetTitle("\\hbox{efektywny pęd wiązki, GeV/c}");
    hBeam_momentum_eff->GetXaxis()->SetTitle("\\hbox{p}_{\\hbox{wiązki}}^{\\hbox{ef}}, \\hbox{GeV/c}");
    hBeam_momentum_eff->DrawCopy("C");

    gXS[35]->Draw("same");
    gXS[41]->Draw("same");
    gXS[55]->Draw("same");
    gXS[75]->Draw("same");
    gXS[89]->Draw("same");
    gXS[150]->Draw("same");

    XS_graph[35]->Draw("same p");
    XS_graph[41]->Draw("same p");
    XS_graph[55]->Draw("same p");
    XS_graph[75]->Draw("same p");
    XS_graph[89]->Draw("same p");

    TLegend *MyLegend01a = new TLegend(0.155, 0.475, 0.345, 0.885);
    MyLegend01a->SetFillStyle(1001); MyLegend01a->SetFillColor(19); MyLegend01a->SetLineColor(0); MyLegend01a->SetTextSize(0.04);
    MyLegend01a->SetBorderSize(5);
    MyLegend01a->AddEntry( XS_graph[35], "EDDA" , "p");
    MyLegend01a->AddEntry((TObject*)0, "SAID:", "");
    MyLegend01a->AddEntry( gXS[35], "#theta^{cm} = 35#circ" , "l");
    MyLegend01a->AddEntry( gXS[41], "#theta^{cm} = 41#circ" , "l");
    MyLegend01a->AddEntry( gXS[55], "#theta^{cm} = 55#circ" , "l");
    MyLegend01a->AddEntry( gXS[75], "#theta^{cm} = 75#circ" , "l");
    MyLegend01a->AddEntry( gXS[89], "#theta^{cm} = 89#circ" , "l");
    MyLegend01a->AddEntry( gXS[150], "#theta^{cm} = 150#circ" , "l");
    MyLegend01a->Draw();

    TLegend *MyLegend01a_pb = new TLegend(0.650, 0.820, 0.885, 0.885);
    MyLegend01a_pb->SetFillStyle(1001); MyLegend01a_pb->SetFillColor(19); MyLegend01a_pb->SetLineColor(0); MyLegend01a_pb->SetTextSize(0.04);
    MyLegend01a_pb->SetBorderSize(5);
    MyLegend01a_pb->AddEntry( hBeam_momentum_eff, "\\hbox{pęd wiązki}" , "l");
    MyLegend01a_pb->Draw();

    MyCanvas01a->Print("output/plots/hCrossSections_pl.png","png");
    MyCanvas01a->Print("output/plots/hCrossSections_pl.eps","eps");

}

