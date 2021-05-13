#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TText.h"

using namespace std;

//Setting color and width of Histogram
void format_h(TH1D* h, int linecolor){
h->SetLineWidth(3);
h->SetLineColor(linecolor);
}

//This code is the third try of Solar neutrino trigger system//
//Adding hep solar neutrino//

int fprompt_compare(){

//////////////////////
//Reading parameter///
//////////////////////

////////////////////////
//Define needed number//
////////////////////////

    int tot_trigger = 0;//Used to count events we consider signal
    int solar_events_mc = 0;//Used to count signal events
    int solar_events_mc_hep = 0;//Used to count hep signal events
    int num_dataS = 0;//Used to count number of photon of signal in each step
    int num_data_Rn222 = 0;//Used to count number of photon of Rn222 alaph decay in each step
    int num_data_hep =0;//Used to count number of photon of hep solar neutrino in each step
    int num_dataS2 = 0;//Used to count number of photon on adjacent PMTs(Sig)
    int num_data_Rn222_2 = 0;
    int num_data_hep_2 = 0;
    int ana_time = 0;//Time prograss of analysis(in the loop)//You could make it 500000
    int reso = 1;
    int fake_trig = 0;//To count number of fake triggers
    int which_event=0;

///////////////////
//Open root files//
///////////////////

    //Solar//
    TFile *fS = new TFile("output/pulse_shape_merged_a-gamma_test_5ns_slow_alpha.root");//500 ns
    //TFile *fS = new TFile("output/pulse_shape_merged_a-gamma_100000_5ns_slow_alpha.root");//5 ns
    TH1D* fprompt_B8 = (TH1D*)fS->Get("Fprompt_S"); //Haven't consider energy yet
    TH1D* fprompt_S_vuv = (TH1D*)fS->Get("Fprompt_S_vuv");
    TH1D* fprompt_Rn222 = (TH1D*)fS->Get("Fprompt_Rn222");
    TH1D* fprompt_Rn222_vuv = (TH1D*)fS->Get("Fprompt_Rn222_vuv");
    TH2D* fprompt_vs_photons_per_event = (TH2D*)fS->Get("Fprompt_vs_photons_per_event_Rn222_vuv");

    //Rn222 alpha//
    //TFile *fRn222 = new TFile("pulse_shape_Rn222_210ns_5ns_reso_vuv.root");
    //TH1D* fprompt_Rn222 = (TH1D*)fRn222->Get("Fprompt_Rn222"); //Haven't consider energy yet

//////////
//Visual//
//////////

//Analysis//
////////////

    double max_different = 0;
    double different = 0;

    double significance[100];
    double significance_vuv[100];
    double max_significance = 0;
    double max_significance_vuv = 0;
    double s_over_b[100];
    double s_over_b_vuv[100];
    double max_s_over_b = 0;
    double max_s_over_b_vuv = 0;
    int max_i_significance = 0;
    int max_i_significance_vuv = 0;
    int max_i_s_over_b = 0;
    int max_i_s_over_b_vuv = 0;

    for(int i=0;i<100;i++){
      significance[i] = fprompt_B8->Integral(0,i)*fprompt_B8->Integral(0,i)/fprompt_Rn222->Integral(0,i);
      s_over_b[i] = fprompt_B8->Integral(0,i)/fprompt_Rn222->Integral(0,i);

      significance_vuv[i] = fprompt_S_vuv->Integral(0,i)*fprompt_S_vuv->Integral(0,i)/fprompt_Rn222_vuv->Integral(0,i);
      s_over_b_vuv[i] = fprompt_S_vuv->Integral(0,i)/fprompt_Rn222_vuv->Integral(0,i);

      if(significance[i]>max_significance){max_significance = significance[i];max_i_significance = i;}
      if(s_over_b[i]>max_s_over_b){max_s_over_b = s_over_b[i];max_i_s_over_b = i;}

      if(significance_vuv[i]>max_significance_vuv){max_significance_vuv = significance_vuv[i];max_i_significance_vuv = i;}
      if(s_over_b_vuv[i]>max_s_over_b_vuv){max_s_over_b_vuv = s_over_b_vuv[i];max_i_s_over_b_vuv = i;}
    }
    cout<<"Max significance (with foils) is "<<max_significance<<" at "<<max_i_significance<<"%"<<endl;
    cout<<"Max signal over BKG (with foils) is "<<max_s_over_b<<" at "<<max_i_s_over_b<<"%"<<endl;

    cout<<"Max significance (without foils) is "<<max_significance_vuv<<" at "<<max_i_significance_vuv<<"%"<<endl;
    cout<<"Max signal over BKG (without foils) is "<<max_s_over_b_vuv<<" at "<<max_i_s_over_b_vuv<<"%"<<endl;

//////////
//Output//
//////////

    //TGraph *running_integral_560KeV_electron_g = new TGraph(400,xaxis,running_integral_B8);
    //running_integral_560KeV_electron_g->Draw();
    //TGraph *running_integral_Rn222_g = new TGraph(400,xaxis,running_integral_Rn222);

    //TFile *fO = new TFile("running_integral.root","RECREATE");

    //running_integral_560KeV_electron_g->Write();
    //running_integral_Rn222_g->Write();
    //fO->Close();

////////
//Draw//
////////

    double scale = 0.2;
    fprompt_Rn222->Scale(scale);
    fprompt_B8->Scale(1.0);
    //TLatex *t1 = new TLatex(0.52,74000,"8B CC electron");
    //TLatex *t2 = new TLatex(0.52,67000,"#alpha from ^{222}Rn");
    TLatex *t1 = new TLatex(0.60,9000,"8B CC electron");
    //TLatex *t2 = new TLatex(0.10,6700,"#alpha from ^{222}Rn");
    TLatex *t2 = new TLatex(0.60,8000,"#alpha -gamma event");
    TLatex *t3 = new TLatex(0.60,7000,"time window = 135ns");
    t1->SetTextColor(kBlue);
    t2->SetTextColor(kRed);
    t3->SetTextColor(kGreen);
    fprompt_B8->SetTitle(";F_{prompt};Events");
    fprompt_Rn222_vuv->SetTitle(";F_{prompt};Events");
    fprompt_B8->SetLineColor(kBlue);
    fprompt_S_vuv->SetLineColor(kBlue);
    fprompt_Rn222->SetLineColor(kRed);
    fprompt_Rn222_vuv->SetLineColor(kRed);
    fprompt_B8->SetLineWidth(3);
    fprompt_S_vuv->SetLineWidth(3);
    fprompt_Rn222_vuv->SetLineWidth(3);
    fprompt_Rn222->SetLineWidth(3);
    //fprompt_B8->Scale(1./fprompt_B8->Integral());
    //fprompt_Rn222->Scale(1./fprompt_Rn222->Integral());
    
    gStyle->SetOptStat(000000000);
    fprompt_B8->SetLabelSize(0.05,"xy");
    fprompt_B8->SetTitleSize(0.05,"xy");
    fprompt_B8->SetTitleOffset(1.1,"y");
    fprompt_S_vuv->SetLabelSize(0.05,"xy");
    fprompt_S_vuv->SetTitleSize(0.05,"xy");

    fprompt_Rn222->SetLabelSize(0.05,"xy");
    fprompt_Rn222->SetTitleSize(0.05,"xy");
    fprompt_Rn222_vuv->SetLabelSize(0.05,"xy");
    fprompt_Rn222_vuv->SetTitleSize(0.05,"xy");

    TCanvas *c = new TCanvas("c","c",1920,1080);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.11);
    c->SetRightMargin(0.11);

    fprompt_vs_photons_per_event->SetLabelSize(0.05,"xy");
    fprompt_vs_photons_per_event->SetTitleSize(0.05,"xy");
    fprompt_vs_photons_per_event->SetTitle(";F_{prompt};Photons collected from each event");
    fprompt_vs_photons_per_event->Draw("Colz");

    fprompt_B8->Draw();
    fprompt_Rn222->Draw("same");
    //fprompt_B8->Draw("same");
    //fprompt_Rn222_vuv->Draw();
    //fprompt_S_vuv->Draw("same");
    t1->Draw("same");t2->Draw("same");t3->Draw("same");
    //c->SaveAs("../plots/timing/PSD_8B_Rn222_100000_5ns.png");

    return 0;
}
