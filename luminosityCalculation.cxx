/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2017-07
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-09
***********************************************/

//Macro to determine the luminosity for the proton-deutron collisions using quasi-free proton-proton scattering with WASA-at-COSY facility

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

//Quadratic background function (for deltaPhi distributions)
const Int_t degPolPhi = 2;
const Int_t nBkgndParam = degPolPhi + 1;
Double_t background(Double_t *x, Double_t *par) {
    Double_t fpolb = 0.;
    for(Int_t i = 0; i < nBkgndParam; i++) {
        fpolb += par[0+i]*pow(x[0],degPolPhi-i);
    };
    return fpolb;
}

//Third degree polynomial function (for luminosity curve)
const Int_t degPol = 3;
const Int_t nParam = degPol + 1;
Double_t polynomialFunc(Double_t *x, Double_t *par) {
    Double_t fpol = 0.;
    for(Int_t i = 0; i < nParam; i++) {
        fpol += par[0+i]*pow(x[0],degPol-i);
    };
    return fpol;
}

//
void luminosityCalculation(){

    cout<<fixed;
    cout.precision(2);

    TFile* myFile[5];

    myFile[0] = new TFile("output/files/ThetaCD_merged.root","READ");
    myFile[1] = new TFile("input/MC-ppn_qf-PARIS-x6.root","READ");
    myFile[2] = new TFile("input/MC-ppn_qf-CDBONN-x6.root","READ");
    myFile[3] = new TFile("input/DATA_ppn_qf_offset.root","READ");

    TH1D *hLuminosity_syst = new TH1D("hLuminosity_syst","",40,-70.,30.);
    TH1D *hLuminosity = new TH1D("hLuminosity","",40,-70.,30.);

    TH1D *hSignalThetaCD;
    TH1D *hBkgndThetaCD;

    TH1D *hGenerated_WMC[12];
    TH1D *hExcessEnergy_WMC_XS[12];
    TH1D *hExcessEnergy_WMC_XS_sq[1];

    TH1D *hExcessEnergy_DATA[12];
    TH1D *hExcessEnergy_DATA_clear[12];
    TH1D *hExcessEnergy_DATA_accept[12];

    TH1D *hDeltaPhi[12];
    TH1D *hDeltaPhi_clear[12];
    TH1D *hDeltaPhi_MC;

    TH1D *hDeltaPhi_Q[12][41];
    TH1D *hDeltaPhi_cut_Q[12][41];
    TH1D *hDeltaPhi_clear_Q[12][41];

    TF1 *fitBkgnd[12][41];

    Double_t integralDATA[12][41];

    Double_t beginFitBkgnd = 0.;
    Double_t endFitBkgnd = 360.;
    Double_t beginDrawBkgnd = 0.;
    Double_t endDrawBkgnd = 360.;

///////////////////////////////////////////

    //Signal and background in thetaCD spectrum
    myFile[0]->cd();
    hSignalThetaCD = (TH1D*)gDirectory->Get("hExcessEnergy_sign");
    hBkgndThetaCD = (TH1D*)gDirectory->Get("hExcessEnergy_bkgnd");

    //WMC
    myFile[1]->cd("Histograms");    //PARIS
    for (Int_t j = 0; j < 11; j++) {
        hGenerated_WMC[j] = (TH1D*)gDirectory->Get("WMC/hGenerated_Q");
    }

    hDeltaPhi_MC = (TH1D*)gDirectory->Get("DATA_lev5/hDeltaPhi_sym_lev5");

    //main cuts
    hExcessEnergy_WMC_XS_sq[0] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_wght_sq_lev5");
    hExcessEnergy_WMC_XS[0] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_wght_bilin_lev5");
    //deltaPhi
    hExcessEnergy_WMC_XS[1] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_wght_bilin_lev5");
    hExcessEnergy_WMC_XS[2] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_wght_bilin_lev5");
    //interpol
    hExcessEnergy_WMC_XS[3] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_wght_lin_lev5");
    hExcessEnergy_WMC_XS[4] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_wght_prx_lev5");
    //new cuts
    hExcessEnergy_WMC_XS[5] = (TH1D*)gDirectory->Get("DATA_syst_lev5_E1/hQ_wght_bilin_syst_lev5_cut1");
    hExcessEnergy_WMC_XS[6] = (TH1D*)gDirectory->Get("DATA_syst_lev5_E2/hQ_wght_bilin_syst_lev5_cut2");
    hExcessEnergy_WMC_XS[7] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T3/hQ_wght_bilin_syst_lev5_cut5");
    hExcessEnergy_WMC_XS[8] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T4/hQ_wght_bilin_syst_lev5_cut6");
    hExcessEnergy_WMC_XS[9] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T5/hQ_wght_bilin_syst_lev5_cut7");
    hExcessEnergy_WMC_XS[10] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T6/hQ_wght_bilin_syst_lev5_cut8");

    //Fermi
    myFile[2]->cd("Histograms");    //CDBONN
    hGenerated_WMC[11] = (TH1D*)gDirectory->Get("WMC/hGenerated_Q");
    hExcessEnergy_WMC_XS[11] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_wght_bilin_lev5");

    //DATA
    myFile[3]->cd("Histograms");

    //main cuts
    hExcessEnergy_DATA[0] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_lev5");
    //deltaPhi
    hExcessEnergy_DATA[1] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_lev5");
    hExcessEnergy_DATA[2] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_lev5");
    //interpol
    hExcessEnergy_DATA[3] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_lev5");
    hExcessEnergy_DATA[4] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_lev5");
    //new cuts
    hExcessEnergy_DATA[5] = (TH1D*)gDirectory->Get("DATA_syst_lev5_E1/hQ_syst_lev5_cut1");
    hExcessEnergy_DATA[6] = (TH1D*)gDirectory->Get("DATA_syst_lev5_E2/hQ_syst_lev5_cut2");
    hExcessEnergy_DATA[7] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T3/hQ_syst_lev5_cut5");
    hExcessEnergy_DATA[8] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T4/hQ_syst_lev5_cut6");
    hExcessEnergy_DATA[9] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T5/hQ_syst_lev5_cut7");
    hExcessEnergy_DATA[10] = (TH1D*)gDirectory->Get("DATA_syst_lev5_T6/hQ_syst_lev5_cut8");
    //Fermi
    hExcessEnergy_DATA[11] = (TH1D*)gDirectory->Get("DATA_lev5/hQ_lev5");

    for(Int_t k = 0; k < 12; k++){

        hExcessEnergy_DATA_clear[k] = new TH1D(Form("hExcessEnergy_DATA_clear_%d",k),"",40,-70.,30.);
        hExcessEnergy_DATA_accept[k] = new TH1D(Form("hExcessEnergy_DATA_accept_%d",k),"",40,-70.,30.);

        hDeltaPhi[k] = new TH1D(Form("hDeltaPhi_%d",k),"",360,0.,360.);
        hDeltaPhi_clear[k] = new TH1D(Form("hDeltaPhi_clear_%d",k),"",360,0.,360.);

        for(Int_t i = 1; i < 41; i++) {   //loop over coplanarity histograms for each Q bin

            //main cuts
            hDeltaPhi_Q[0][i]=(TH1D*)gDirectory->Get(Form("DATA_lev5/DeltaPhi_lev5/hDeltaPhi_lev5_Q%d",i));  //Delta Phi for i-th bin before bkgnd rejection
            //deltaPhi
            hDeltaPhi_Q[1][i]=(TH1D*)gDirectory->Get(Form("DATA_lev5/DeltaPhi_lev5/hDeltaPhi_lev5_Q%d",i));
            hDeltaPhi_Q[2][i]=(TH1D*)gDirectory->Get(Form("DATA_lev5/DeltaPhi_lev5/hDeltaPhi_lev5_Q%d",i));
            //interpol
            hDeltaPhi_Q[3][i]=(TH1D*)gDirectory->Get(Form("DATA_lev5/DeltaPhi_lev5/hDeltaPhi_lev5_Q%d",i));
            hDeltaPhi_Q[4][i]=(TH1D*)gDirectory->Get(Form("DATA_lev5/DeltaPhi_lev5/hDeltaPhi_lev5_Q%d",i));
            //new cuts
            hDeltaPhi_Q[5][i]=(TH1D*)gDirectory->Get(Form("DATA_syst_lev5_E1/DeltaPhi_syst_lev5_cut1/hDeltaPhi_syst_lev5_cut1_Q%d",i));
            hDeltaPhi_Q[6][i]=(TH1D*)gDirectory->Get(Form("DATA_syst_lev5_E2/DeltaPhi_syst_lev5_cut2/hDeltaPhi_syst_lev5_cut2_Q%d",i));
            hDeltaPhi_Q[7][i]=(TH1D*)gDirectory->Get(Form("DATA_syst_lev5_T3/DeltaPhi_syst_lev5_cut5/hDeltaPhi_syst_lev5_cut5_Q%d",i));
            hDeltaPhi_Q[8][i]=(TH1D*)gDirectory->Get(Form("DATA_syst_lev5_T4/DeltaPhi_syst_lev5_cut6/hDeltaPhi_syst_lev5_cut6_Q%d",i));
            hDeltaPhi_Q[9][i]=(TH1D*)gDirectory->Get(Form("DATA_syst_lev5_T5/DeltaPhi_syst_lev5_cut7/hDeltaPhi_syst_lev5_cut7_Q%d",i));
            hDeltaPhi_Q[10][i]=(TH1D*)gDirectory->Get(Form("DATA_syst_lev5_T6/DeltaPhi_syst_lev5_cut8/hDeltaPhi_syst_lev5_cut8_Q%d",i));
            //Fermi
            hDeltaPhi_Q[11][i]=(TH1D*)gDirectory->Get(Form("DATA_lev5/DeltaPhi_lev5/hDeltaPhi_lev5_Q%d",i));

            hDeltaPhi_cut_Q[k][i]= new TH1D(Form("hDeltaPhi_cut_Q%d_%d",i,k),"",360,0.,360.);
            //hDelta_Phi_clear_Q[k][i]= (TH1D*)hDeltaPhi_Q[k][i]->Clone(Form("hDeltaPhi_clear_Q%d_%d",i,k));

            Double_t binContent = 0.;

            for(Int_t j = 1; j < 361; j++) {
                binContent = hDeltaPhi_Q[k][i]->GetBinContent(j);
                //hDeltaPhi_cut_Q[k][i] = (TH1D*)hDeltaPhi_Q[k][i]->Clone(Form("hDeltaPhi_cut_Q%d_%d",i,k));

                if (k == 1) {   //a//
                    if (j > 80 && j < 280) {
                        hDeltaPhi_cut_Q[k][i]->SetBinContent(j,0.);
                    }
                    else {
                        hDeltaPhi_cut_Q[k][i]->SetBinContent(j,binContent);
                    }
                }                   //a//
                else {
                    if (k == 2) {   //b//
                        if (j > 100 && j < 260) {
                            hDeltaPhi_cut_Q[k][i]->SetBinContent(j,0.);
                        }
                        else {
                            hDeltaPhi_cut_Q[k][i]->SetBinContent(j,binContent);
                        }
                    }               //b//
                    else {          //c//
                        if (j > 90 && j < 270) {
                            hDeltaPhi_cut_Q[k][i]->SetBinContent(j,0.);
                        }
                        else {
                            hDeltaPhi_cut_Q[k][i]->SetBinContent(j,binContent);
                        }
                    }               //c//
                }
            }

            hDeltaPhi_clear_Q[k][i] = (TH1D*)hDeltaPhi_Q[k][i]->Clone(Form("hDeltaPhi_clear_Q%d_%d",i,k));

            //fit
            fitBkgnd[k][i] = new TF1(Form("fitBkgnd_%d_%d",k,i),background,beginDrawBkgnd,endDrawBkgnd,nBkgndParam);

            hDeltaPhi_cut_Q[k][i]->Fit(Form("fitBkgnd_%d_%d",k,i),"","",beginFitBkgnd,endFitBkgnd);

            hDeltaPhi_clear_Q[k][i]->Add(fitBkgnd[k][i],-1);

            //all Q bins together
            hDeltaPhi[k]->Add(hDeltaPhi_Q[k][i]);
            hDeltaPhi_clear[k]->Add(hDeltaPhi_clear_Q[k][i]);

            //
            integralDATA[k][i] = hDeltaPhi_clear_Q[k][i]->Integral();
            hExcessEnergy_DATA_clear[k]->SetBinContent(i,integralDATA[k][i]);

        }

        hExcessEnergy_DATA_accept[k] = (TH1D*)hExcessEnergy_DATA_clear[k]->Clone(Form("hExcessEnergy_DATA_accept_%d",k));
        hExcessEnergy_DATA_accept[k]->Add(hBkgndThetaCD,-1);

    }

    //////////////////////////////////////////////
/*
    TH1D* hLuminosityPart1 = new TH1D("hLuminosityPart1","",40,-70.,30.);

    hLuminosityPart1->Divide(hGenerated_WMC[0],hExcessEnergy_WMC_XS[0],1/(2*TMath::Pi()),1e6,"");   //1e6 - calculation of XS in [nb]
    hLuminosity->Multiply(hLuminosityPart1,hExcessEnergy_DATA_accept[0],1.,4e3);                    //4000 - prescaling factor for Tr. 21
*/
    Double_t luminosity[12];

    //for statistical error estimation
    Double_t Ngen[12],Nexp[12],Nwmc[12];
    Double_t errNgen, errNexp, errNwmc;
    Double_t errStatLuminPart[3];
    Double_t errStatLumin;
    Double_t errStatLuminSq = 0.;

    //for systematics
    Double_t Ndata,Nbkgr;
    Double_t errSystLuminPart[12];
    Double_t errSystLuminPartTotal[12];
    Double_t errSystLumin;
    Double_t errSystLuminTotal = 0.;

    errSystLuminPartTotal[0] = 0.;
    errSystLuminPartTotal[1] = 0.;
    errSystLuminPartTotal[2] = 0.;
    errSystLuminPartTotal[3] = 0.;
    errSystLuminPartTotal[4] = 0.;
    errSystLuminPartTotal[5] = 0.;
    errSystLuminPartTotal[6] = 0.;
    errSystLuminPartTotal[7] = 0.;
    errSystLuminPartTotal[8] = 0.;

    for (Int_t l = 1; l < 41; l++) {

        for (Int_t m = 0; m < 12; m++) {

            Ngen[m] = hGenerated_WMC[m]->GetBinContent(l);
            Nexp[m] = hExcessEnergy_DATA_accept[m]->GetBinContent(l);
            Nwmc[m] = hExcessEnergy_WMC_XS[m]->GetBinContent(l);

            luminosity[m] = (1/(2*TMath::Pi()))*(Ngen[m]/(Nwmc[m]*1e6))*Nexp[m]*4000;

        }

        //errNgen = hGenerated_WMC[0]->GetBinError(l);
        errNexp = hExcessEnergy_DATA_accept[0]->GetBinError(l);
        errNwmc = hExcessEnergy_WMC_XS_sq[0]->GetBinError(l);

        //errStatLuminPart[0] = (Ngen/Nwmc)*errNgen;
        errStatLuminPart[1] = (Ngen[0]/Nwmc[0])*errNexp;
        errStatLuminPart[2] = ((Ngen[0]*Nexp[0])/(Nwmc[0]*Nwmc[0]))*errNwmc;

        //errStatLumin = (4000./(2*TMath::Pi()*1e6))*TMath::Sqrt(errStatLuminPart[0]*errStatLuminPart[0] + errStatLuminPart[1]*errStatLuminPart[1] + errStatLuminPart[2]*errStatLuminPart[2]);
        errStatLumin = (4000./(2*TMath::Pi()*1e6))*TMath::Sqrt(errStatLuminPart[1]*errStatLuminPart[1] + errStatLuminPart[2]*errStatLuminPart[2]);

        hLuminosity->SetBinContent(l,luminosity[0]);
        hLuminosity->SetBinError(l,errStatLumin);

        //for error of total integrated luminosity
        errStatLuminSq = errStatLuminSq + errStatLumin*errStatLumin;

        //systematics
        Ndata = hExcessEnergy_DATA[0]->GetBinContent(l);
        Nbkgr = hBkgndThetaCD->GetBinContent(l);

        errSystLuminPart[0] = 0.5*(TMath::Abs(luminosity[1]-luminosity[0]) + TMath::Abs(luminosity[2]-luminosity[0]));  //deltaPhi
        errSystLuminPart[1] = 0.5*(TMath::Abs(luminosity[3]-luminosity[0]) + TMath::Abs(luminosity[4]-luminosity[0]));  //interpol
        errSystLuminPart[2] = 0.5*(TMath::Abs(luminosity[5]-luminosity[0]) + TMath::Abs(luminosity[6]-luminosity[0]));  //energy
        errSystLuminPart[3] = 0.5*(TMath::Abs(luminosity[7]-luminosity[0]) + TMath::Abs(luminosity[8]-luminosity[0]));  //angle 1
        errSystLuminPart[4] = 0.5*(TMath::Abs(luminosity[9]-luminosity[0]) + TMath::Abs(luminosity[10]-luminosity[0])); //angle 2
        errSystLuminPart[5] = 0.5*TMath::Abs(luminosity[11]-luminosity[0]); //Fermi
        errSystLuminPart[6] = 0.5*(Nbkgr/Ndata)*luminosity[0];              //background in thetaCD
        errSystLuminPart[7] = 0.0225*luminosity[0];                         //shadow effect
        errSystLuminPart[8] = 0.027*luminosity[0];                          //pp cross sections

        Double_t err = 0;
        for (Int_t n = 0; n < 9; n++) {
            err = err + TMath::Power(errSystLuminPart[n],2);
            errSystLuminPartTotal[n] = errSystLuminPartTotal[n] + errSystLuminPart[n];
        }

        errSystLumin = TMath::Sqrt(err);

        hLuminosity_syst->SetBinContent(l,luminosity[0]);
        hLuminosity_syst->SetBinError(l,errSystLumin);

        errSystLuminTotal = errSystLuminTotal + errSystLumin;

    }

    Double_t luminosityTotal = hLuminosity->Integral();      //total luminosity
    Double_t errStatLuminTotal = TMath::Sqrt(errStatLuminSq);   //total statistical error

    Double_t errStatLuminAngleCut = TMath::Sqrt(TMath::Power(errSystLuminPartTotal[3],2) + TMath::Power(errSystLuminPartTotal[4],2));
    Double_t errStatLuminCuts = TMath::Sqrt(TMath::Power(errSystLuminPartTotal[2],2) + TMath::Power(errSystLuminPartTotal[3],2) + TMath::Power(errSystLuminPartTotal[4],2));
/*
    cout<<"Integrated luminosity is equal to "<<luminosityTotal<<" nb^{-1}"<<endl;
    cout<<"Statistical uncertainties: "<<errStatLuminTotal<<" nb^{-1} ("<<(errStatLuminTotal/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Systematic uncertainties: "<<errSystLuminTotal<<" nb^{-1} ("<<(errSystLuminTotal/luminosityTotal)*100<<"%)"<<endl;
    cout<<""<<endl;
    cout<<"Error from deltaPhi: "<<errSystLuminPartTotal[0]<<" ("<<(errSystLuminPartTotal[0]/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from interpol: "<<errSystLuminPartTotal[1]<<" ("<<(errSystLuminPartTotal[1]/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from energy cut: "<<errSystLuminPartTotal[2]<<" ("<<(errSystLuminPartTotal[2]/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from angle cut: "<<errStatLuminAngleCut<<" ("<<(errStatLuminAngleCut/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from cuts: "<<errStatLuminCuts<<" ("<<(errStatLuminCuts/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from Fermi distr: "<<errSystLuminPartTotal[5]<<" ("<<(errSystLuminPartTotal[5]/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from bkgnd in thetaCD: "<<errSystLuminPartTotal[6]<<" ("<<(errSystLuminPartTotal[6]/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from shadow eff: "<<errSystLuminPartTotal[7]<<" ("<<(errSystLuminPartTotal[7]/luminosityTotal)*100<<"%)"<<endl;
    cout<<"Error from pp cross sect: "<<errSystLuminPartTotal[8]<<" ("<<(errSystLuminPartTotal[8]/luminosityTotal)*100<<"%)"<<endl;
*/
    ofstream results;
    results.open("output/files/LuminosityQuasyFreeReaction.dat", ios::app);
    results<<Form("Integrated luminosity is equal to %.0f",luminosityTotal)<<endl;
    results<<Form("Statistical uncertainties: %.0f (%.1f%)",errStatLuminTotal,(errStatLuminTotal/luminosityTotal)*100)<<endl;
    results<<Form("Systematic from cuts: %.0f (%.1f%)",errStatLuminCuts,(errStatLuminCuts/luminosityTotal)*100)<<endl;
    results<<Form("Systematic from bkgnd rejection in thetaCD: %.0f (%.1f%)",errSystLuminPartTotal[6],(errSystLuminPartTotal[6]/luminosityTotal)*100)<<endl;
    results<<Form("Systematic from bkgnd rejection in deltaPhi: %.0f (%.1f%)",errSystLuminPartTotal[0],(errSystLuminPartTotal[0]/luminosityTotal)*100)<<endl;
    results<<Form("Systematic from Fermi mom distr: %.0f (%.1f%)",errSystLuminPartTotal[5],(errSystLuminPartTotal[5]/luminosityTotal)*100)<<endl;
    results<<Form("Systematic from interpolation: %.0f (%.1f%)",errSystLuminPartTotal[1],(errSystLuminPartTotal[1]/luminosityTotal)*100)<<endl;
    results<<Form("Systematic from shadow effect: %.0f (%.1f%)",errSystLuminPartTotal[7],(errSystLuminPartTotal[7]/luminosityTotal)*100)<<endl;
    results<<Form("Systematic from pp cross sections: %.0f (%.1f%)",errSystLuminPartTotal[8],(errSystLuminPartTotal[8]/luminosityTotal)*100)<<endl;
    results<<Form("Total systematic uncertainties: %.0f (%.1f%)",errSystLuminTotal,(errSystLuminTotal/luminosityTotal)*100)<<endl;
    results<<Form("Normalization error (from EDDA): %.0f (4%)",0.04*luminosityTotal)<<endl;
    results.close();

    //create new output file.root
    TFile* newFile = new TFile("output/files/LuminosityQuasyFreeReaction.root","RECREATE");
    newFile->cd();

    hLuminosity->Write("");
    hLuminosity_syst->Write("");
    hSignalThetaCD->Write("");
    hBkgndThetaCD->Write("");
    hGenerated_WMC[0]->Write("hGenerated_WMC_PARIS");
    hGenerated_WMC[11]->Write("hGenerated_WMC_CDBONN");
    for (Int_t j = 0; j < 12; j++) {
        hExcessEnergy_WMC_XS[j]->Write(Form("hExcessEnergy_WMC_XS_%d",j));
        hExcessEnergy_DATA[j]->Write(Form("hExcessEnergy_DATA_%d",j));
        hExcessEnergy_DATA_clear[j]->Write(Form("hExcessEnergy_DATA_clear_%d",j));
        hExcessEnergy_DATA_accept[j]->Write(Form("hExcessEnergy_DATA_accept_%d",j));
        hDeltaPhi[j]->Write(Form("hDeltaPhi_%d",j));
        hDeltaPhi_clear[j]->Write(Form("hDeltaPhi_clear_%d",j));

        for (Int_t i = 1; i < 41; i++) {
            hDeltaPhi_Q[j][i]->Write(Form("hDeltaPhi_Q%d_%d",i,j));
            hDeltaPhi_clear_Q[j][i]->Write(Form("hDeltaPhi_clear_Q%d_%d",i,j));
        }
    }
    newFile->Close();

    //fit function to intefrated luminosity
    Double_t beginFitLum = -70.; //ustawić
    Double_t endFitLum = 30.;
    Double_t beginDrawLum = -70;
    Double_t endDrawLum = 30.;

    TF1 *fitLum = new TF1("fitLum",polynomialFunc,beginFitLum,endFitLum,nParam);

    hLuminosity->Fit(fitLum,"","",beginFitLum,endFitLum);

    Double_t P[4];
    Double_t errP[4];

    P[0] = fitLum->GetParameter(0)*10e4;
    P[1] = fitLum->GetParameter(1)*10e3;
    P[2] = fitLum->GetParameter(2)*10e1;
    P[3] = fitLum->GetParameter(3);

    errP[0] = fitLum->GetParError(0)*10e4;
    errP[1] = fitLum->GetParError(1)*10e3;
    errP[2] = fitLum->GetParError(2)*10e1;
    errP[3] = fitLum->GetParError(3);

    //cout<<"p_0 = "<<P[0]<<" +- "<<errP[0]<<endl;
    //cout<<"p_1 = "<<P[1]<<" +- "<<errP[1]<<endl;
    //cout<<"p_2 = "<<P[2]<<" +- "<<errP[2]<<endl;
    //cout<<"p_3 = "<<P[3]<<" +- "<<errP[3]<<endl;


    ////
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    //gStyle->SetPadTopMargin(0.10);
    gStyle->SetPalette(55);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTextFont(62);

    //
    TCanvas* MyCanvas00 = new TCanvas;

    Double_t maxY00 = 1.2*(hExcessEnergy_DATA_accept[0]->GetMaximum());

    hExcessEnergy_DATA_accept[0]->SetTitle("Accepted experimental events");
    hExcessEnergy_DATA_accept[0]->GetXaxis()->SetTitle("excess energy [MeV]");
    hExcessEnergy_DATA_accept[0]->GetXaxis()->SetTitleOffset(1.);
    hExcessEnergy_DATA_accept[0]->GetXaxis()->SetTitleSize(0.06);
    hExcessEnergy_DATA_accept[0]->GetXaxis()->SetLabelSize(0.05);
    hExcessEnergy_DATA_accept[0]->GetYaxis()->SetTitle("counts");
    hExcessEnergy_DATA_accept[0]->GetYaxis()->SetTitleOffset(0.9);
    hExcessEnergy_DATA_accept[0]->GetYaxis()->SetTitleSize(0.06);
    hExcessEnergy_DATA_accept[0]->GetYaxis()->SetLabelSize(0.05);
    //hExcessEnergy_DATA_accept[0]->GetXaxis()->SetRangeUser(-70.,30.);
    hExcessEnergy_DATA_accept[0]->GetYaxis()->SetRangeUser(0.,maxY00);

    hExcessEnergy_DATA_accept[0]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[0]->SetLineColor(1);
    hExcessEnergy_DATA_accept[0]->SetMarkerStyle(23);
    hExcessEnergy_DATA_accept[0]->SetMarkerColor(1);
    hExcessEnergy_DATA_accept[0]->SetMarkerSize(0.8);
    hExcessEnergy_DATA_accept[0]->Draw("E1");

    hExcessEnergy_DATA_accept[1]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[1]->SetLineColor(51);
    hExcessEnergy_DATA_accept[1]->Draw("same C");

    hExcessEnergy_DATA_accept[2]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[2]->SetLineColor(56);
    hExcessEnergy_DATA_accept[2]->Draw("same C");

    hExcessEnergy_DATA_accept[5]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[5]->SetLineColor(63);
    hExcessEnergy_DATA_accept[5]->Draw("same C");

    hExcessEnergy_DATA_accept[6]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[6]->SetLineColor(68);
    hExcessEnergy_DATA_accept[6]->Draw("same C");

    hExcessEnergy_DATA_accept[7]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[7]->SetLineColor(77);
    hExcessEnergy_DATA_accept[7]->Draw("same C");

    hExcessEnergy_DATA_accept[8]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[8]->SetLineColor(90);
    hExcessEnergy_DATA_accept[8]->Draw("same C");

    hExcessEnergy_DATA_accept[9]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[9]->SetLineColor(94);
    hExcessEnergy_DATA_accept[9]->Draw("same C");

    hExcessEnergy_DATA_accept[10]->SetLineWidth(1);
    hExcessEnergy_DATA_accept[10]->SetLineColor(98);
    hExcessEnergy_DATA_accept[10]->Draw("same C");

    TLegend *legend00 = new TLegend(0.585, 0.215, 0.885, 0.600);
    legend00->SetFillStyle(1); legend00->SetFillColor(0); legend00->SetLineColor(0); legend00->SetTextSize(0.04);
    legend00->AddEntry( hExcessEnergy_DATA_accept[0], "main", "pe");
    legend00->AddEntry( hExcessEnergy_DATA_accept[1], "#Delta#phi #in (80#circ,280#circ)", "l");
    legend00->AddEntry( hExcessEnergy_DATA_accept[2], "#Delta#phi #in (100#circ,260#circ)", "l");
    legend00->AddEntry( hExcessEnergy_DATA_accept[5], "energy cut +10%", "l");
    legend00->AddEntry( hExcessEnergy_DATA_accept[6], "energy cut -10%", "l");
    legend00->AddEntry( hExcessEnergy_DATA_accept[7], "#theta_{CD} #in (34#circ,100#circ)", "l");
    legend00->AddEntry( hExcessEnergy_DATA_accept[8], "#theta_{CD} #in (46#circ,100#circ)", "l");
    legend00->AddEntry( hExcessEnergy_DATA_accept[9], "#theta_{CD} #in (40#circ,106#circ)", "l");
    legend00->AddEntry( hExcessEnergy_DATA_accept[10], "#theta_{CD} #in (40#circ,94#circ)", "l");
    legend00->Draw("same");

    MyCanvas00->Print("output/plots/hExcessEnergy_DATA_accept.png","png");
    MyCanvas00->Print("output/plots/hExcessEnergy_DATA_accept.eps","eps");

    //
    TCanvas* MyCanvas01 = new TCanvas;

    Double_t maxY01 = 1.2*(hGenerated_WMC[0]->GetMaximum());

    hGenerated_WMC[0]->SetTitle("Generated events");
    hGenerated_WMC[0]->GetXaxis()->SetTitle("excess energy [MeV]");
    hGenerated_WMC[0]->GetXaxis()->SetTitleOffset(1.);
    hGenerated_WMC[0]->GetXaxis()->SetTitleSize(0.06);
    hGenerated_WMC[0]->GetXaxis()->SetLabelSize(0.05);
    hGenerated_WMC[0]->GetYaxis()->SetTitle("counts");
    hGenerated_WMC[0]->GetYaxis()->SetTitleOffset(0.9);
    hGenerated_WMC[0]->GetYaxis()->SetTitleSize(0.06);
    hGenerated_WMC[0]->GetYaxis()->SetLabelSize(0.05);
    //hGenerated_WMC[0]->GetXaxis()->SetRangeUser(-70.,30.);
    hGenerated_WMC[0]->GetYaxis()->SetRangeUser(0.,maxY01);

    hGenerated_WMC[0]->SetLineWidth(1);
    hGenerated_WMC[0]->SetLineColor(kRed+1);
    hGenerated_WMC[0]->SetMarkerStyle(23);
    hGenerated_WMC[0]->SetMarkerColor(kRed+1);
    hGenerated_WMC[0]->SetMarkerSize(0.8);
    hGenerated_WMC[0]->Draw("E1");

    hGenerated_WMC[11]->SetLineWidth(1);
    hGenerated_WMC[11]->SetLineColor(kAzure-3);
    hGenerated_WMC[11]->SetMarkerStyle(24);
    hGenerated_WMC[11]->SetMarkerColor(kAzure-3);
    hGenerated_WMC[11]->SetMarkerSize(0.8);
    hGenerated_WMC[11]->Draw("same E1");

    TLegend *legend01 = new TLegend(0.710, 0.530, 0.885, 0.675);
    legend01->SetFillStyle(1); legend01->SetFillColor(0); legend01->SetLineColor(0); legend01->SetTextSize(0.04);
    legend01->AddEntry( hGenerated_WMC[0], "PARIS", "pe");
    legend01->AddEntry( hGenerated_WMC[11], "CDBONN", "pe");
    legend01->Draw("same");

    MyCanvas01->Print("output/plots/hGenerated_WMC.png","png");
    MyCanvas01->Print("output/plots/hGenerated_WMC.eps","eps");

    //
    TCanvas* MyCanvas02 = new TCanvas;

    Double_t maxY02 = 1.2*(hExcessEnergy_WMC_XS[0]->GetMaximum());

    hExcessEnergy_WMC_XS[0]->SetTitle("Accepted and weighted simulated events");
    hExcessEnergy_WMC_XS[0]->GetXaxis()->SetTitle("excess energy [MeV]");
    hExcessEnergy_WMC_XS[0]->GetXaxis()->SetTitleOffset(1.);
    hExcessEnergy_WMC_XS[0]->GetXaxis()->SetTitleSize(0.06);
    hExcessEnergy_WMC_XS[0]->GetXaxis()->SetLabelSize(0.05);
    hExcessEnergy_WMC_XS[0]->GetYaxis()->SetTitle("counts");
    hExcessEnergy_WMC_XS[0]->GetYaxis()->SetTitleOffset(0.9);
    hExcessEnergy_WMC_XS[0]->GetYaxis()->SetTitleSize(0.06);
    hExcessEnergy_WMC_XS[0]->GetYaxis()->SetLabelSize(0.05);
    //hExcessEnergy_WMC_XS[0]->GetXaxis()->SetRangeUser(-70.,30.);
    hExcessEnergy_WMC_XS[0]->GetYaxis()->SetRangeUser(0.,maxY02);

    hExcessEnergy_WMC_XS[0]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[0]->SetLineColor(1);
    hExcessEnergy_WMC_XS[0]->SetMarkerStyle(23);
    hExcessEnergy_WMC_XS[0]->SetMarkerColor(1);
    hExcessEnergy_WMC_XS[0]->SetMarkerSize(0.8);
    hExcessEnergy_WMC_XS[0]->Draw("E1");

    hExcessEnergy_WMC_XS[3]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[3]->SetLineColor(51);
    hExcessEnergy_WMC_XS[3]->Draw("same C");

    hExcessEnergy_WMC_XS[4]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[4]->SetLineColor(56);
    hExcessEnergy_WMC_XS[4]->Draw("same C");

    hExcessEnergy_WMC_XS[5]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[5]->SetLineColor(63);
    hExcessEnergy_WMC_XS[5]->Draw("same C");

    hExcessEnergy_WMC_XS[6]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[6]->SetLineColor(68);
    hExcessEnergy_WMC_XS[6]->Draw("same C");

    hExcessEnergy_WMC_XS[7]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[7]->SetLineColor(77);
    hExcessEnergy_WMC_XS[7]->Draw("same C");

    hExcessEnergy_WMC_XS[8]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[8]->SetLineColor(94);
    hExcessEnergy_WMC_XS[8]->Draw("same C");

    hExcessEnergy_WMC_XS[9]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[9]->SetLineColor(90);
    hExcessEnergy_WMC_XS[9]->Draw("same C");

    hExcessEnergy_WMC_XS[10]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[10]->SetLineColor(98);
    hExcessEnergy_WMC_XS[10]->Draw("same C");

    hExcessEnergy_WMC_XS[11]->SetLineWidth(1);
    hExcessEnergy_WMC_XS[11]->SetLineColor(30);
    hExcessEnergy_WMC_XS[11]->Draw("same C");

    TLegend *legend02 = new TLegend(0.555, 0.215, 0.885, 0.625);
    legend02->SetFillStyle(1); legend02->SetFillColor(0); legend02->SetLineColor(0); legend02->SetTextSize(0.04);
    legend02->AddEntry( hExcessEnergy_WMC_XS[0], "main", "pe");
    legend02->AddEntry( hExcessEnergy_WMC_XS[3], "linear interpolation", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[4], "closest points", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[5], "energy cut +10%", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[6], "energy cut -10%", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[7], "#theta_{CD} #in (34#circ,100#circ)", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[8], "#theta_{CD} #in (46#circ,100#circ)", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[8], "#theta_{CD} #in (40#circ,106#circ)", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[8], "#theta_{CD} #in (40#circ,94#circ)", "l");
    legend02->AddEntry( hExcessEnergy_WMC_XS[11], "CDBonn", "l");
    legend02->Draw("same");

    MyCanvas02->Print("output/plots/hExcessEnergy_WMC_XS.png","png");
    MyCanvas02->Print("output/plots/hExcessEnergy_WMC_XS.eps","eps");

    //
    TCanvas* MyCanvas03 = new TCanvas;

    Double_t maxY03 = 1.2*(hLuminosity->GetMaximum());

    hLuminosity_syst->SetTitle("Luminosity");
    hLuminosity_syst->GetXaxis()->SetTitle("excess energy [MeV]");
    hLuminosity_syst->GetXaxis()->SetTitleOffset(1.);
    hLuminosity_syst->GetXaxis()->SetTitleSize(0.06);
    hLuminosity_syst->GetXaxis()->SetLabelSize(0.05);
    hLuminosity_syst->GetYaxis()->SetTitle("integrated luminosity [nb^{-1}]");
    hLuminosity_syst->GetYaxis()->SetTitleOffset(0.8);
    hLuminosity_syst->GetYaxis()->SetTitleSize(0.06);
    hLuminosity_syst->GetYaxis()->SetLabelSize(0.05);
    //hLuminosity_syst->GetXaxis()->SetRangeUser(-70.,30.);
    hLuminosity_syst->GetYaxis()->SetRangeUser(0.,maxY03);

    hLuminosity_syst->SetLineWidth(1);
    hLuminosity_syst->SetLineColor(2);
    hLuminosity_syst->SetMarkerStyle(2);
    hLuminosity_syst->SetMarkerColor(2);
    hLuminosity_syst->SetMarkerSize(1);
    //hLuminosity_syst->SetFillStyle(3345);
    //hLuminosity_syst->SetFillColor(2);
    hLuminosity_syst->Draw("E1");

    hLuminosity->SetLineWidth(1);
    hLuminosity->SetLineColor(1);
    hLuminosity->SetMarkerStyle(23);
    hLuminosity->SetMarkerColor(1);
    hLuminosity->SetMarkerSize(0.8);
    hLuminosity->Draw("same E1");

    fitLum->SetLineWidth(2);
    fitLum->SetLineColor(kCyan-3);
    fitLum->SetLineStyle(1);
    fitLum->Draw("same C");

    TString Param[4];

    Param[0]=Form("a = (%.2f #pm %.2f)#upoint10^{-5} [nb^{-1}#upointMeV^{-3}]",P[0],errP[0]);
    Param[1]=Form("b = (%.2f #pm %.2f)#upoint10^{-4} [nb^{-1}#upointMeV^{-2}]",P[1],errP[1]);
    Param[2]=Form("c = (%.2f #pm %.2f)#upoint10^{-2} [nb^{-1}#upointMeV^{-1}]",P[2],errP[2]);
    Param[3]=Form("d = (%.2f #pm %.2f) [nb^{-1}]",P[3],errP[3]);

    TLegend *legend03 = new TLegend(0.150, 0.345, 0.480, 0.575);
    legend03->SetFillStyle(1); legend03->SetFillColor(0); legend03->SetLineColor(0); legend03->SetTextSize(0.04);
    legend03->AddEntry( hLuminosity, "luminosity", "pe");
    legend03->AddEntry( hLuminosity_syst, "systematic errors", "pe");
    legend03->AddEntry( fitLum, "fitting function", "l");
    legend03->AddEntry((TObject*)0, "aQ^{3} + bQ^{2} + cQ + d", "");
    //legend03->AddEntry( fitFcn, "fitting function (aQ^{3}+bQ^{2}+cQ+d)", "l");
    legend03->Draw("same");

    TPaveText *descr = new TPaveText(-19.,11.6,27.3,45.4);
    TText *tl01 = descr->AddText(Form("L = (%.0f #pm %.0f #pm %.0f #pm %.0f) [nb^{-1}]",luminosityTotal,errStatLuminTotal,errSystLuminTotal,0.04*luminosityTotal));
    TText *tl02 = descr->AddText(Param[0]);
    TText *tl03 = descr->AddText(Param[1]);
    TText *tl04 = descr->AddText(Param[2]);
    TText *tl05 = descr->AddText(Param[3]);
    descr->SetFillColor(0);
    tl01->SetTextFont(42);
    tl02->SetTextFont(42);
    tl03->SetTextFont(42);
    tl04->SetTextFont(42);
    tl05->SetTextFont(42);
    tl01->SetTextSize(0.035);
    tl02->SetTextSize(0.035);
    tl03->SetTextSize(0.035);
    tl04->SetTextSize(0.035);
    tl05->SetTextSize(0.035);
    tl01->SetTextAlign(22);
    tl02->SetTextAlign(12);
    tl03->SetTextAlign(12);
    tl04->SetTextAlign(12);
    tl05->SetTextAlign(12);
    descr->Draw("same");

    TPaveText *text_lev01 = new TPaveText(-62.,67.,-62.,67.,"capt");
    text_lev01->SetTextFont(42); text_lev01->SetTextSize(0.06);
    text_lev01->SetTextAlign(22);
    text_lev01->SetFillStyle(0);
    text_lev01->SetShadowColor(0); text_lev01->SetFillColor(0);
    text_lev01->SetBorderSize(0);
    text_lev01->AddText("(b)");
    text_lev01->Draw();

    MyCanvas03->Print("output/plots/hLuminosity.png","png");
    MyCanvas03->Print("output/plots/hLuminosity.eps","eps");

    //
    TCanvas* MyCanvas03a = new TCanvas;

    hLuminosity_syst->SetTitle("");
    hLuminosity_syst->GetXaxis()->SetTitle("\\hbox{energia dostępna [MeV]}");
    hLuminosity_syst->GetYaxis()->SetTitle("\\hbox{świetlność całkowita [} \\hbox{nb}^{\\hbox{-1}}\\hbox{]}");
    hLuminosity_syst->Draw("E1");
    hLuminosity->Draw("same E1");
    fitLum->Draw("same C");

    TLegend *legend03a = new TLegend(0.425, 0.385, 0.885, 0.615);
    legend03a->SetFillStyle(1001); legend03a->SetFillColor(19); legend03a->SetLineColor(1); legend03a->SetTextSize(0.04); legend03a->SetBorderSize(5);
    legend03a->AddEntry( hLuminosity, "\\hbox{świetlność}", "pe");
    legend03a->AddEntry( hLuminosity_syst, "\\hbox{niepewności systematyczne}", "pe");
    legend03a->AddEntry( fitLum, "dopasowana funkcja", "l");
    legend03a->AddEntry((TObject*)0, "aQ^{3} + bQ^{2} + cQ + d", "");
    legend03a->Draw();

    //descr->Draw("same");

    MyCanvas03a->Print("output/plots/hLuminosity_pl.png","png");
    MyCanvas03a->Print("output/plots/hLuminosity_pl.eps","eps");

    //
    TCanvas* MyCanvas04=new TCanvas;

    for(Int_t i = 1; i < 41; i++) {

        hDeltaPhi_Q[0][i]->SetTitle(Form("Q #in (%G,%G) MeV",-70.+(i-1)*2.5,-67.5+(i-1)*2.5));
        hDeltaPhi_Q[0][i]->GetXaxis()->SetTitle("(2#pi + #Delta#phi)mod2#pi [deg]");
        hDeltaPhi_Q[0][i]->GetYaxis()->SetTitle("counts");
        hDeltaPhi_Q[0][i]->GetXaxis()->SetTitleSize(0.06);
        hDeltaPhi_Q[0][i]->GetXaxis()->SetTitleOffset(1.0);
        hDeltaPhi_Q[0][i]->GetXaxis()->SetLabelSize(0.05);
        hDeltaPhi_Q[0][i]->GetYaxis()->SetTitleSize(0.06);
        hDeltaPhi_Q[0][i]->GetYaxis()->SetTitleOffset(1.1);
        hDeltaPhi_Q[0][i]->GetYaxis()->SetLabelSize(0.05);
        //hDeltaPhi_Q[0][i]->GetXaxis()->SetRangeUser(90,270);
        hDeltaPhi_Q[0][i]->SetLineWidth(2);
        hDeltaPhi_Q[0][i]->SetLineColor(1);
        hDeltaPhi_Q[0][i]->SetMarkerColor(1);
        hDeltaPhi_Q[0][i]->SetMarkerSize(1);
        hDeltaPhi_Q[0][i]->SetMarkerStyle(2);
        hDeltaPhi_Q[0][i]->Draw("C");

        hDeltaPhi_clear_Q[0][i]->SetLineColor(kOrange+7);
        hDeltaPhi_clear_Q[0][i]->SetLineWidth(1);
        hDeltaPhi_clear_Q[0][i]->SetFillColor(kOrange+7);
        hDeltaPhi_clear_Q[0][i]->SetFillStyle(3354);
        hDeltaPhi_clear_Q[0][i]->Draw("same LF2");

        fitBkgnd[0][i]->SetLineColor(kCyan-3);
        fitBkgnd[0][i]->SetLineWidth(2);
        fitBkgnd[0][i]->Draw("same");

        TLegend *legend04 = new TLegend(0.585, 0.700, 0.890, 0.885);
        legend04->SetFillStyle(1); legend04->SetFillColor(0); legend04->SetLineColor(0); legend04->SetTextSize(0.04);
        legend04->AddEntry( hDeltaPhi_Q[0][i], "experimental data", "l");
        legend04->AddEntry( fitBkgnd[0][i], "fitting function", "l");
        legend04->AddEntry( hDeltaPhi_clear_Q[0][i], "signal", "f");
        legend04->Draw("same");

        MyCanvas04->Print(Form("output/plots/hDeltaPhi_Q%d.png",i),"png");
        MyCanvas04->Print(Form("output/plots/hDeltaPhi_Q%d.eps",i),"eps");

    }

    //
    TCanvas* MyCanvas04a=new TCanvas;

    for(Int_t i = 1; i < 41; i++) {

        hDeltaPhi_Q[0][i]->SetTitle(Form("Q #in (%G,%G) MeV",-70.+(i-1)*2.5,-67.5+(i-1)*2.5));
        hDeltaPhi_Q[0][i]->GetXaxis()->SetTitle("(2#pi + #Delta#phi)mod2#pi [#circ]");
        hDeltaPhi_Q[0][i]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
        hDeltaPhi_Q[0][i]->Draw("C");

        hDeltaPhi_clear_Q[0][i]->Draw("same LF2");
        fitBkgnd[0][i]->Draw("same");

        TLegend *legend04a = new TLegend(0.535, 0.700, 0.890, 0.885);
        legend04a->SetFillStyle(1001); legend04a->SetFillColor(19); legend04a->SetLineColor(1); legend04a->SetTextSize(0.04); legend04a->SetBorderSize(5);
        legend04a->AddEntry( hDeltaPhi_Q[0][i], "dane eksperymentalne", "l");
        legend04a->AddEntry( fitBkgnd[0][i], "dopasowana funkcja", "l");
        legend04a->AddEntry( hDeltaPhi_clear_Q[0][i], "\\hbox{sygnał}", "f");
        legend04a->Draw();

        MyCanvas04a->Print(Form("output/plots/hDeltaPhi_Q%d_pl.png",i),"png");
        MyCanvas04a->Print(Form("output/plots/hDeltaPhi_Q%d_pl.eps",i),"eps");

    }

    //
    TCanvas* MyCanvas05=new TCanvas;

    hDeltaPhi[0]->SetTitle("Complanarity");
    hDeltaPhi[0]->GetXaxis()->SetTitle("(2#pi + #Delta#phi)mod2#pi [deg]");
    hDeltaPhi[0]->GetXaxis()->SetTitleOffset(0.9);
    hDeltaPhi[0]->GetXaxis()->SetTitleSize(0.06);
    hDeltaPhi[0]->GetXaxis()->SetLabelSize(0.05);
    hDeltaPhi[0]->GetYaxis()->SetTitle("counts");
    hDeltaPhi[0]->GetYaxis()->SetTitleOffset(0.9);
    hDeltaPhi[0]->GetYaxis()->SetTitleSize(0.06);
    hDeltaPhi[0]->GetYaxis()->SetLabelSize(0.05);
    hDeltaPhi[0]->GetYaxis()->SetRangeUser(0,1.07*(hDeltaPhi[0]->GetMaximum()));
    //hDeltaPhi[0]->GetXaxis()->SetRangeUser(90,270);

    hDeltaPhi[0]->SetLineWidth(2);
    hDeltaPhi[0]->SetLineColor(1);
    hDeltaPhi[0]->SetMarkerColor(1);
    hDeltaPhi[0]->SetMarkerStyle(2);
    hDeltaPhi[0]->SetMarkerSize(1);
    hDeltaPhi[0]->Draw("C");

    hDeltaPhi_clear[0]->SetLineWidth(1);
    hDeltaPhi_clear[0]->SetLineColor(kOrange+7);
    hDeltaPhi_clear[0]->SetFillColor(kOrange+7);
    hDeltaPhi_clear[0]->SetFillStyle(3354);
    hDeltaPhi_clear[0]->Draw("same LF2");

    TLegend *legend05 = new TLegend(0.585, 0.760, 0.890, 0.885);
    legend05->SetFillStyle(1); legend05->SetFillColor(0); legend05->SetLineColor(0); legend05->SetTextSize(0.04);
    legend05->AddEntry( hDeltaPhi[0], "experimental data", "l");
    legend05->AddEntry( hDeltaPhi_clear[0], "signal", "f");
    legend05->Draw();

    MyCanvas05->Print("output/plots/hDeltaPhi.png","png");
    MyCanvas05->Print("output/plots/hDeltaPhi.eps","eps");

    //
    TCanvas* MyCanvas05a=new TCanvas;

    hDeltaPhi[0]->SetTitle("");
    hDeltaPhi[0]->GetXaxis()->SetTitle("(2#pi + #Delta#phi)mod2#pi [#circ]");
    hDeltaPhi[0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hDeltaPhi[0]->Draw("C");

    hDeltaPhi_clear[0]->Draw("same LF2");

    TLegend *legend05a = new TLegend(0.535, 0.768, 0.890, 0.885);
    legend05a->SetFillStyle(1001); legend05a->SetFillColor(19); legend05a->SetLineColor(1); legend05a->SetTextSize(0.04); legend05a->SetBorderSize(5);
    legend05a->AddEntry( hDeltaPhi[0], "dane eksperymentalne", "l");
    legend05a->AddEntry( hDeltaPhi_clear[0], "\\hbox{sygnał}", "f");
    legend05a->Draw();

    MyCanvas05a->Print("output/plots/hDeltaPhi_pl.png","png");
    MyCanvas05a->Print("output/plots/hDeltaPhi_pl.eps","eps");

    //
    TCanvas* MyCanvas06=new TCanvas;

    //hDeltaPhi_MC->SetTitle("Complanarity");
    hDeltaPhi_MC->GetXaxis()->SetTitle("(2#pi + #Delta#phi)mod2#pi [deg]");
    hDeltaPhi_MC->GetXaxis()->SetTitleOffset(0.9);
    hDeltaPhi_MC->GetXaxis()->SetTitleSize(0.06);
    hDeltaPhi_MC->GetXaxis()->SetLabelSize(0.05);
    hDeltaPhi_MC->GetYaxis()->SetTitle("counts");
    hDeltaPhi_MC->GetYaxis()->SetTitleOffset(0.9);
    hDeltaPhi_MC->GetYaxis()->SetTitleSize(0.06);
    hDeltaPhi_MC->GetYaxis()->SetLabelSize(0.05);
    hDeltaPhi_MC->GetYaxis()->SetRangeUser(0,1.07*(hDeltaPhi_MC->GetMaximum()));
    //hDeltaPhi_MC->GetXaxis()->SetRangeUser(90,270);

    hDeltaPhi_MC->SetLineWidth(2);
    hDeltaPhi_MC->SetLineColor(kAzure-3);
    hDeltaPhi_MC->SetFillColor(kAzure-3);
    hDeltaPhi_MC->SetFillStyle(3354);
    hDeltaPhi_MC->Draw("LF2");

    TLegend *legend06 = new TLegend(0.560, 0.830, 0.890, 0.885);
    legend06->SetFillStyle(1); legend06->SetFillColor(0); legend06->SetLineColor(0); legend06->SetTextSize(0.04);
    legend06->AddEntry( hDeltaPhi_MC, "numerical data WMC", "f");
    legend06->Draw();

    MyCanvas06->Print("output/plots/hDeltaPhi_MC.png","png");
    MyCanvas06->Print("output/plots/hDeltaPhi_MC.eps","eps");

    //
    TCanvas* MyCanvas06a=new TCanvas;

    hDeltaPhi_MC->SetTitle("");
    hDeltaPhi_MC->GetXaxis()->SetTitle("(2#pi + #Delta#phi)mod2#pi [#circ]");
    hDeltaPhi_MC->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hDeltaPhi_MC->Draw("LF2");

    hDeltaPhi_MC->SetLineWidth(1);

    TLegend *legend06a = new TLegend(0.565, 0.815, 0.890, 0.885);
    legend06a->SetFillStyle(1001); legend06a->SetFillColor(19); legend06a->SetLineColor(1); legend06a->SetTextSize(0.04); legend06a->SetBorderSize(5);
    legend06a->AddEntry( hDeltaPhi_MC, "WMC: pd #rightarrow ppn_{sp}", "f");
    legend06a->Draw();

    MyCanvas06a->Print("output/plots/hDeltaPhi_MC_pl.png","png");
    MyCanvas06a->Print("output/plots/hDeltaPhi_MC_pl.eps","eps");

}
