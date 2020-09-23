/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2017-07
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-06
***********************************************/

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
#include <TMathText.h>
#include <TLatex.h>

void drawHistograms_pl() {

    TFile* myFile[5];

    myFile[0] = new TFile("input/DATA_ppn_qf_offset.root","READ");
    myFile[1] = new TFile("input/MC-ppn_qf-PARIS-x6.root","READ");
    //myFile[1] = new TFile("input/MC-ppn_qf-CDBONN-x6.root","READ");
    myFile[2] = new TFile("input/MC-dnpi_plus_thetacut-x6.root","READ");
    myFile[3] = new TFile("input/MC-pd-momcut.root","READ");

    ///////////////////////////////////////////////

    TH2F* hEdepPSBvsSEC[5];
    TH2F* hThetaFDvsThetaCD[5];

    for (int i = 0; i < 2; ++i) {

        myFile[i]->cd("Histograms");

        hEdepPSBvsSEC[i] = (TH2F*)gDirectory->Get("DATA_lev3/hEdepPSBvsSEC_lev3");
        hThetaFDvsThetaCD[i] = (TH2F*)gDirectory->Get("DATA_lev4/hThetaFDvsThetaCD_lev4");

    }

    myFile[2]->cd("Histograms");
    hEdepPSBvsSEC[2] = (TH2F*)gDirectory->Get("DATA_lev4/hEdepPSBvsSEC_lev4");

    myFile[3]->cd("Histograms");
    hThetaFDvsThetaCD[2] = (TH2F*)gDirectory->Get("DATA_lev4/hThetaFDvsThetaCD_lev4");

    //
    double x[2];
    double y[2];
    x[0] = 0.0010;
    y[1] = 0.00015;

    y[0] = -0.00618182*x[0] + 0.003400;
    x[1] = (0.00015 - 0.003400)/(-0.00618182);

    ////
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetPadRightMargin(0.133);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPalette(55);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTextFont(62);    

    TCanvas* MyCanvas01 = new TCanvas("MyCanvas01","",465,500);

    hEdepPSBvsSEC[0]->GetXaxis()->SetTitle("#DeltaE(SEC), GeV");
    hEdepPSBvsSEC[0]->GetYaxis()->SetTitle("#DeltaE(PSB), GeV");
    hEdepPSBvsSEC[0]->GetXaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[0]->GetXaxis()->SetTitleOffset(1.0);
    hEdepPSBvsSEC[0]->GetXaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[0]->GetYaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[0]->GetYaxis()->SetTitleOffset(1.55);
    hEdepPSBvsSEC[0]->GetYaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[0]->GetXaxis()->SetRangeUser(0.0010,0.7);
    hEdepPSBvsSEC[0]->GetYaxis()->SetRangeUser(0.00015,0.012);
    hEdepPSBvsSEC[0]->GetZaxis()->SetLabelSize(0.05);
    //hEdepPSBvsSEC[0]->SetMaximum(8000.);
    //hEdepPSBvsSEC[0]->RebinX(2);
    //hEdepPSBvsSEC[0]->RebinY(2);
    //gPad->SetLogz();
    hEdepPSBvsSEC[0]->Draw("colz");

    TPaveText *pt01 = new TPaveText(0., 0.010, 0.602, 0.012);
    TText *t01 = pt01->AddText("Dane eksperymentalne");
    pt01->SetFillColor(19);
    t01->SetTextSize(0.06);
    //pt01->Draw();

    TPaveText *pt010 = new TPaveText(0.05, 0.005, 0.05, 0.005,"PROTON");
    pt010->SetTextSize(0.06);
    pt010->SetFillColor(0);
    pt010->SetTextColor(0);
    pt010->SetTextAlign(22);
    pt010->AddText("p");

    TPaveText *pt011 = new TPaveText(0.05, 0.0025, 0.05, 0.0025,"PION");
    pt011->SetTextSize(0.07);
    pt011->SetFillColor(0);
    pt011->SetTextColor(0);
    pt011->SetTextAlign(22);
    pt011->AddText("#pi^{+}");

    TLine* line010 = new TLine(x[0],y[0],x[1],y[1]);
    line010->SetLineColor(kRed);
    line010->SetLineWidth(2);
    line010->Draw("same");

    pt010->Draw("same");
    pt011->Draw("same");

    TLegend *MyLegend00 = new TLegend(0.385, 0.770, 0.855, 0.885);
    MyLegend00->SetFillStyle(1001); MyLegend00->SetFillColor(19); MyLegend00->SetLineColor(1); MyLegend00->SetBorderSize(5);
    MyLegend00->SetTextSize(0.04);
    MyLegend00->AddEntry((TObject*)0, "Dane eksperyment.", "");
    MyLegend00->AddEntry(line010, "\\hbox{cięcie}", "l");
    MyLegend00->Draw();

    MyCanvas01->Print("output/plots/hEdepPSBvsSEC_DATA_pl.png","png");
    MyCanvas01->Print("output/plots/hEdepPSBvsSEC_DATA_pl.eps","eps");

    //
    TCanvas* MyCanvas02 = new TCanvas("MyCanvas02","",465,500);

    hEdepPSBvsSEC[1]->GetXaxis()->SetTitle("#DeltaE(SEC), GeV");
    hEdepPSBvsSEC[1]->GetYaxis()->SetTitle("#DeltaE(PSB), GeV");
    hEdepPSBvsSEC[1]->GetXaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[1]->GetXaxis()->SetTitleOffset(1.0);
    hEdepPSBvsSEC[1]->GetXaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[1]->GetYaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[1]->GetYaxis()->SetTitleOffset(1.55);
    hEdepPSBvsSEC[1]->GetYaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[1]->GetXaxis()->SetRangeUser(0.0010,0.7);
    hEdepPSBvsSEC[1]->GetYaxis()->SetRangeUser(0.00015,0.012);
    hEdepPSBvsSEC[1]->GetZaxis()->SetLabelSize(0.05);
    //hEdepPSBvsSEC[1]->SetMaximum(8000.);
    //hEdepPSBvsSEC[1]->RebinX(2);
    //hEdepPSBvsSEC[1]->RebinY(2);
    //gPad->SetLogz();
    hEdepPSBvsSEC[1]->Draw("colz");

    TPaveText *pt02 = new TPaveText(0., 0.010, 0.474, 0.012);
    TText *t02 = pt02->AddText("WMC: pd #rightarrow ppn_{sp}");
    pt02->SetFillColor(19);
    t02->SetTextSize(0.06);
    //pt02->Draw();

    TPaveText *pt020 = new TPaveText(0.05, 0.005, 0.05, 0.005,"PROTON");
    pt020->SetTextSize(0.06);
    pt020->SetFillColor(0);
    pt020->SetTextColor(0);
    pt020->SetTextAlign(22);
    pt020->AddText("p");

    TPaveText *pt021 = new TPaveText(0.05, 0.0025, 0.05, 0.0025,"PION");
    pt021->SetTextSize(0.06);
    pt021->SetFillColor(0);
    pt021->SetTextColor(0);
    pt021->SetTextAlign(22);
    pt021->AddText("#pi^{+}");

    TLine* line020 = new TLine(x[0],y[0],x[1],y[1]);
    line020->SetLineColor(kRed);
    line020->SetLineWidth(2);
    line020->Draw("same");

    pt020->Draw("same");
    pt021->Draw("same");

    TLegend *MyLegend01 = new TLegend(0.385, 0.770, 0.855, 0.885);
    MyLegend01->SetFillStyle(1001); MyLegend01->SetFillColor(19); MyLegend01->SetLineColor(1); MyLegend01->SetBorderSize(5);
    MyLegend01->SetTextSize(0.04);
    MyLegend01->AddEntry((TObject*)0, "WMC: pd #rightarrow ppn_{sp}", "");
    MyLegend01->AddEntry(line020, "\\hbox{cięcie}", "l");
    MyLegend01->Draw();

    MyCanvas02->Print("output/plots/hEdepPSBvsSEC_MC_pl.png","png");
    MyCanvas02->Print("output/plots/hEdepPSBvsSEC_MC_pl.eps","eps");

    //
    TCanvas* MyCanvas03 = new TCanvas("MyCanvas03","",465,500);

    hEdepPSBvsSEC[2]->GetXaxis()->SetTitle("#DeltaE(SEC), GeV");
    hEdepPSBvsSEC[2]->GetYaxis()->SetTitle("#DeltaE(PSB), GeV");
    hEdepPSBvsSEC[2]->GetXaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[2]->GetXaxis()->SetTitleOffset(1.0);
    hEdepPSBvsSEC[2]->GetXaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[2]->GetYaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[2]->GetYaxis()->SetTitleOffset(1.55);
    hEdepPSBvsSEC[2]->GetYaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[2]->GetXaxis()->SetRangeUser(0.0010,0.5);
    hEdepPSBvsSEC[2]->GetYaxis()->SetRangeUser(0.00015,0.012);
    hEdepPSBvsSEC[2]->GetZaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[2]->GetXaxis()->SetNdivisions(5,5,0, kTRUE);
    //hEdepPSBvsSEC[2]->SetMaximum(8000.);
    //hEdepPSBvsSEC[2]->RebinX(2);
    //hEdepPSBvsSEC[2]->RebinY(2);
    //gPad->SetLogz();
    hEdepPSBvsSEC[2]->Draw("colz");

    TPaveText *pt03 = new TPaveText(0., 0.010, 0.315, 0.012);
    TText *t03 = pt03->AddText("WMC: pd #rightarrow d#pi^{+}n_{sp}");
    pt03->SetFillColor(19);
    t03->SetTextSize(0.06);
    t03->SetTextFont(42);
    //pt03->Draw();

    TPaveText *pt030 = new TPaveText(0.1, 0.008, 0.1, 0.008,"DEUTERON");
    pt030->SetTextSize(0.06);
    pt030->SetFillColor(0);
    pt030->SetTextColor(0);
    pt030->SetTextAlign(22);
    pt030->AddText("d");

    TPaveText *pt031 = new TPaveText(0.05, 0.0025, 0.05, 0.0025,"PION");
    pt031->SetTextSize(0.07);
    pt031->SetFillColor(0);
    pt031->SetTextColor(0);
    pt031->SetTextAlign(22);
    pt031->AddText("#pi^{+}");

    TLine* line030 = new TLine(x[0],y[0],x[1],y[1]);
    line030->SetLineColor(kRed);
    line030->SetLineWidth(2);
    //line030->Draw("same");

    //pt030->Draw("same");
    pt031->Draw("same");

    TLegend *MyLegend02 = new TLegend(0.385, 0.815, 0.855, 0.885);
    MyLegend02->SetFillStyle(1001); MyLegend02->SetFillColor(19); MyLegend02->SetLineColor(1); MyLegend02->SetBorderSize(5);
    MyLegend02->SetTextSize(0.04);
    MyLegend02->AddEntry((TObject*)0, "WMC: pd #rightarrow d#pi^{+}n_{sp}", "");
    //MyLegend02->AddEntry(line020, "\\hbox{cięcie}", "l");
    MyLegend02->Draw();

    MyCanvas03->Print("output/plots/hEdepPSBvsSEC_pi_pl.png","png");
    MyCanvas03->Print("output/plots/hEdepPSBvsSEC_pi_pl.eps","eps");

    //
    TCanvas* MyCanvas04 = new TCanvas("MyCanvas04","",465,500);

    hThetaFDvsThetaCD[0]->GetXaxis()->SetTitle("#theta_{FD},#circ");
    hThetaFDvsThetaCD[0]->GetYaxis()->SetTitle("#theta_{CD},#circ");
    hThetaFDvsThetaCD[0]->GetXaxis()->SetTitleSize(0.06);
    hThetaFDvsThetaCD[0]->GetXaxis()->SetTitleOffset(1.);
    hThetaFDvsThetaCD[0]->GetXaxis()->SetLabelSize(0.05);
    hThetaFDvsThetaCD[0]->GetYaxis()->SetTitleSize(0.06);
    hThetaFDvsThetaCD[0]->GetYaxis()->SetTitleOffset(1.3);
    hThetaFDvsThetaCD[0]->GetYaxis()->SetLabelSize(0.05);
    hThetaFDvsThetaCD[0]->GetXaxis()->SetRangeUser(3.,18.);
    hThetaFDvsThetaCD[0]->GetYaxis()->SetRangeUser(20.,170.);
    hThetaFDvsThetaCD[0]->GetZaxis()->SetLabelSize(0.05);
    //hThetaFDvsThetaCD[0]->SetMaximum(8000.);
    //hThetaFDvsThetaCD[0]->RebinX(2);
    //hThetaFDvsThetaCD[0]->RebinY(2);
    gPad->SetLogz();
    hThetaFDvsThetaCD[0]->Draw("colz");

    TLine* line040 = new TLine(3.,40.,18.,40.);
    line040->SetLineColor(kRed);
    line040->SetLineWidth(2);
    line040->Draw("same");

    TLine* line041 = new TLine(3.,100.,18.,100.);
    line041->SetLineColor(kRed);
    line041->SetLineWidth(2);
    line041->Draw("same");

    TPaveText *pt04 = new TPaveText(3.,150.,15.7,170.);
    TText *t04 = pt04->AddText("Dane eksperymentalne");
    pt04->SetFillColor(19);
    t04->SetTextSize(0.06);
    //pt04->Draw();

    TPaveText *pt040 = new TPaveText(16., 65., 16., 65.,"A");
    pt040->SetTextSize(0.05);
    pt040->SetFillColor(0);
    pt040->SetTextColor(0);
    pt040->SetTextAlign(22);
    pt040->AddText("a");
    pt040->Draw("same");

    TPaveText *pt041 = new TPaveText(15., 82., 15., 82.,"B");
    pt041->SetTextSize(0.06);
    pt041->SetFillColor(0);
    pt041->SetTextColor(0);
    pt041->SetTextAlign(22);
    pt041->AddText("b");
    pt041->Draw("same");

    TPaveText *pt042 = new TPaveText(14., 125., 14., 125.,"C");
    pt042->SetTextSize(0.06);
    pt042->SetFillColor(0);
    pt042->SetTextColor(0);
    pt042->SetTextAlign(22);
    pt042->AddText("c");
    pt042->Draw("same");

    TLegend *MyLegend04 = new TLegend(0.385, 0.770, 0.855, 0.885);
    MyLegend04->SetFillStyle(1001); MyLegend04->SetFillColor(19); MyLegend04->SetLineColor(1); MyLegend04->SetBorderSize(5);
    MyLegend04->SetTextSize(0.04);
    MyLegend04->AddEntry((TObject*)0, "Dane eksperyment.", "");
    MyLegend04->AddEntry(line040, "\\hbox{cięcie}", "l");
    MyLegend04->Draw();

    MyCanvas04->Print("output/plots/hThetaFDvsThetaCD_DATA_pl.png","png");
    MyCanvas04->Print("output/plots/hThetaFDvsThetaCD_DATA_pl.eps","eps");

    //
    TCanvas* MyCanvas05 = new TCanvas("MyCanvas05","",465,500);

    hThetaFDvsThetaCD[1]->GetXaxis()->SetTitle("#theta_{FD},#circ");
    hThetaFDvsThetaCD[1]->GetYaxis()->SetTitle("#theta_{CD},#circ");
    hThetaFDvsThetaCD[1]->GetXaxis()->SetTitleSize(0.06);
    hThetaFDvsThetaCD[1]->GetXaxis()->SetTitleOffset(1.);
    hThetaFDvsThetaCD[1]->GetXaxis()->SetLabelSize(0.05);
    hThetaFDvsThetaCD[1]->GetYaxis()->SetTitleSize(0.06);
    hThetaFDvsThetaCD[1]->GetYaxis()->SetTitleOffset(1.3);
    hThetaFDvsThetaCD[1]->GetYaxis()->SetLabelSize(0.05);
    hThetaFDvsThetaCD[1]->GetXaxis()->SetRangeUser(3.,18.);
    hThetaFDvsThetaCD[1]->GetYaxis()->SetRangeUser(20.,170.);
    hThetaFDvsThetaCD[1]->GetZaxis()->SetLabelSize(0.05);
    //hThetaFDvsThetaCD[1]->SetMaximum(8000.);
    //hThetaFDvsThetaCD[1]->RebinX(2);
    //hThetaFDvsThetaCD[1]->RebinY(2);
    gPad->SetLogz();
    hThetaFDvsThetaCD[1]->Draw("colz");

    TLine* line050 = new TLine(3.,40.,18.,40.);
    line050->SetLineColor(kRed);
    line050->SetLineWidth(2);
    line050->Draw("same");

    TLine* line051 = new TLine(3.,100.,18.,100.);
    line051->SetLineColor(kRed);
    line051->SetLineWidth(2);
    line051->Draw("same");

    TPaveText *pt05 = new TPaveText(3.,150.,13.4,170.);
    TText *t05 = pt05->AddText("WMC: pd #rightarrow ppn_{sp}");
    pt05->SetFillColor(19);
    t05->SetTextSize(0.06);
    //pt05->Draw();

    TLegend *MyLegend05 = new TLegend(0.385, 0.770, 0.855, 0.885);
    MyLegend05->SetFillStyle(1001); MyLegend05->SetFillColor(19); MyLegend05->SetLineColor(1); MyLegend05->SetBorderSize(5);
    MyLegend05->SetTextSize(0.04);
    MyLegend05->AddEntry((TObject*)0, "WMC: pd #rightarrow ppn_{sp}", "");
    MyLegend05->AddEntry(line050, "\\hbox{cięcie}", "l");
    MyLegend05->Draw();

    MyCanvas05->Print("output/plots/hThetaFDvsThetaCD_MC_pl.png","png");
    MyCanvas05->Print("output/plots/hThetaFDvsThetaCD_MC_pl.eps","eps");

    //
    TCanvas* MyCanvas06 = new TCanvas("MyCanvas06","",465,500);

    hThetaFDvsThetaCD[2]->GetXaxis()->SetTitle("#theta_{FD},#circ");
    hThetaFDvsThetaCD[2]->GetYaxis()->SetTitle("#theta_{CD},#circ");
    hThetaFDvsThetaCD[2]->GetXaxis()->SetTitleSize(0.06);
    hThetaFDvsThetaCD[2]->GetXaxis()->SetTitleOffset(1.);
    hThetaFDvsThetaCD[2]->GetXaxis()->SetLabelSize(0.05);
    hThetaFDvsThetaCD[2]->GetYaxis()->SetTitleSize(0.06);
    hThetaFDvsThetaCD[2]->GetYaxis()->SetTitleOffset(1.3);
    hThetaFDvsThetaCD[2]->GetYaxis()->SetLabelSize(0.05);
    hThetaFDvsThetaCD[2]->GetXaxis()->SetRangeUser(3.,18.);
    hThetaFDvsThetaCD[2]->GetYaxis()->SetRangeUser(20.,170.);
    hThetaFDvsThetaCD[2]->GetZaxis()->SetLabelSize(0.05);
    //hThetaFDvsThetaCD[2]->SetMaximum(8000.);
    //hThetaFDvsThetaCD[2]->RebinX(2);1111
    //hThetaFDvsThetaCD[2]->RebinY(2);
    gPad->SetLogz();
    hThetaFDvsThetaCD[2]->Draw("colz");

    TPaveText *pt060 = new TPaveText(15., 82., 15., 82.,"PROTON");
    pt060->SetTextSize(0.06);
    pt060->SetFillColor(0);
    pt060->SetTextColor(1);
    pt060->SetTextAlign(22);
    pt060->AddText("p");
    pt060->Draw("same");

    TPaveText *pt061 = new TPaveText(14., 125., 14., 125.,"DEUTERON");
    pt061->SetTextSize(0.06);
    pt061->SetFillColor(0);
    pt061->SetTextColor(1);
    pt061->SetTextAlign(22);
    pt061->AddText("d");
    pt061->Draw("same");

    TLine* line060 = new TLine(3.,40.,18.,40.);
    line060->SetLineColor(kRed);
    line060->SetLineWidth(2);
    //line060->Draw("same");

    TLine* line061 = new TLine(3.,100.,18.,100.);
    line061->SetLineColor(kRed);
    line061->SetLineWidth(2);
    //line061->Draw("same");

    TPaveText *pt06 = new TPaveText(3.,150.,12.0,170.);
    TText *t06 = pt06->AddText("WMC: pd #rightarrow pd");
    pt06->SetFillColor(19);
    t06->SetTextSize(0.06);
    //pt06->Draw();

    TLegend *MyLegend06 = new TLegend(0.385, 0.815, 0.855, 0.885);
    MyLegend06->SetFillStyle(1001); MyLegend06->SetFillColor(19); MyLegend06->SetLineColor(1); MyLegend06->SetBorderSize(5);
    MyLegend06->SetTextSize(0.04);
    MyLegend06->AddEntry((TObject*)0, "WMC: pd #rightarrow d#pi^{+}n_{sp}", "");
    //MyLegend06->AddEntry(line020, "\\hbox{cięcie}", "l");
    MyLegend06->Draw();

    MyCanvas06->Print("output/plots/hThetaFDvsThetaCD_pd_pl.png","png");
    MyCanvas06->Print("output/plots/hThetaFDvsThetaCD_pd_pl.eps","eps");

}
