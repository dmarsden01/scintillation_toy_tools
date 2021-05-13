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

using namespace std;

//This code is to calculate the Fprompt of Rn222 from it's pulse shape//

int merge_pulse_shape_upgrade(){

//////////////////////
//Reading parameter///
//////////////////////

////////////////////////
//Define needed number//
////////////////////////

    const int NUM_BINS = 2000;//10000[ns]/NUM_BINS is the timing resolution
    int reso = 10000/NUM_BINS;//Timing resolution
    bool photon_cut = false;
    int photon_cut_num = 10;

///////////////////
//Open root files//
///////////////////


    //Solar//
    TFile *fS = new TFile("../new_data/8B_100000_uniform.root");
    TTree* data_treeS = (TTree*)fS->Get("data_tree");//If it works on total photon, it works on vuv
    TTree* data_treeS_vuv = (TTree*)fS->Get("data_tree_vuv");//Without foils
    TTree* event_treeS = (TTree*)fS->Get("event_tree"); //Haven't consider energy yet

    int data_eventS;int data_eventS_vuv;int data_pmtS; double data_timeS;double data_timeS_vuv;double event_decay_timeS;double event_ES;double event_x_posS;//[cm]

    data_treeS->SetBranchAddress("data_time", &data_timeS);
    data_treeS->SetBranchAddress("data_event", &data_eventS);
    data_treeS_vuv->SetBranchAddress("data_time_vuv", &data_timeS_vuv);//VUV only
    data_treeS_vuv->SetBranchAddress("data_event_vuv", &data_eventS_vuv);
    event_treeS->SetBranchAddress("event_decay_time", &event_decay_timeS);
    event_treeS->SetBranchAddress("event_E", &event_ES);
    event_treeS->SetBranchAddress("event_x_pos", &event_x_posS);

    //Rn222//
    //TFile *fRn222 = new TFile("../new_data/a_gamma_100000_uniform_full_scint.root");
    TFile *fRn222 = new TFile("../new_data/a-gamma_test_3.root");
    TTree* data_treeRn222 = (TTree*)fRn222->Get("data_tree");//If it works on total photon, it works on vuv
    TTree* data_treeRn222_vuv = (TTree*)fRn222->Get("data_tree_vuv");//Without foils
    TTree* event_treeRn222 = (TTree*)fRn222->Get("event_tree"); //Haven't consider energy yet

    int data_eventRn222;int data_eventRn222_vuv;int data_pmtRn222; double data_timeRn222;double data_timeRn222_vuv;double event_decay_timeRn222;double event_ERn222;double event_x_posRn222;

    data_treeRn222->SetBranchAddress("data_time", &data_timeRn222);
    data_treeRn222->SetBranchAddress("data_event", &data_eventRn222);
    data_treeRn222_vuv->SetBranchAddress("data_time_vuv", &data_timeRn222_vuv);//VUV only
    data_treeRn222_vuv->SetBranchAddress("data_event_vuv", &data_eventRn222_vuv);
    event_treeRn222->SetBranchAddress("event_decay_time", &event_decay_timeRn222);
    event_treeRn222->SetBranchAddress("event_E", &event_ERn222);
    event_treeRn222->SetBranchAddress("event_x_pos", &event_x_posRn222);

    int new_data_eventRn222; 
    int new_data_eventRn222_vuv;

//////////
//Visual//
//////////

    TH1D *pulse_shape_S = new TH1D("pulse_shape_S","pulse_shape_S;time_of_photons_from_event[us];entries",NUM_BINS,0,10);
    TH1D *pulse_shape_S_vuv = new TH1D("pulse_shape_S_vuv","pulse_shape_S_vuv;time_of_photons_from_event[us];entries",NUM_BINS,0,10);
    TH1D *fprompt_S_h = new TH1D("Fprompt_S",";Fprompt;entries",100,0,1);
    TH1D *fprompt_S_vuv_h = new TH1D("Fprompt_S_vuv",";Fprompt;entries",100,0,1);
    TH2D *fprompt_vs_energy_S = new TH2D("Fprompt_vs_energy_S","Fprompt_vs_energy_S;Fprompt;event energy[MeV]",100,0,1,100,0,20);
    TH2D *fprompt_vs_photons_per_event_S = new TH2D("Fprompt_vs_photons_per_event_S","Fprompt_vs_photons_per_event_S;Fprompt;Number of photons",100,0,1,300,0,300);
    TH2D *fprompt_vs_energy_S_vuv = new TH2D("Fprompt_vs_energy_S_vuv","Fprompt_vs_energy_S_vuv;Fprompt;event energy[MeV]",100,0,1,100,0,20);
    TH2D *fprompt_vs_photons_per_event_S_vuv = new TH2D("Fprompt_vs_photons_per_event_S_vuv","Fprompt_vs_photons_per_event_S_vuv;Fprompt;Number of photons",100,0,1,300,0,300);

    TH1D *pulse_shape_Rn222 = new TH1D("pulse_shape_Rn222","pulse_shape_Rn222;time_of_photons_from_event[us];entries",NUM_BINS,0,10);//Each bin is 5ns
    TH1D *pulse_shape_Rn222_vuv = new TH1D("pulse_shape_Rn222_vuv","pulse_shape_Rn222_vuv;time_of_photons_from_event[us];entries",NUM_BINS,0,10);//Each bin is 5ns
    TH1D *fprompt_Rn222_h = new TH1D("Fprompt_Rn222",";Fprompt;entries",100,0,1);
    TH1D *fprompt_Rn222_vuv_h = new TH1D("Fprompt_Rn222_vuv",";Fprompt;entries",100,0,1);
    TH2D *fprompt_vs_energy_Rn222 = new TH2D("Fprompt_vs_energy_Rn222","Fprompt_vs_energy_Rn222;Fprompt;event energy[MeV]",100,0,1,100,0,20);
    TH2D *fprompt_vs_photons_per_event_Rn222 = new TH2D("Fprompt_vs_photons_per_event_Rn222","Fprompt_vs_photons_per_event_Rn222;Fprompt;Number of photons",100,0,1,300,0,300);
    TH2D *fprompt_vs_energy_Rn222_vuv = new TH2D("Fprompt_vs_energy_Rn222_vuv","Fprompt_vs_energy_Rn222_vuv;Fprompt;event energy[MeV]",100,0,1,100,0,20);
    TH2D *fprompt_vs_photons_per_event_Rn222_vuv = new TH2D("Fprompt_vs_photons_per_event_Rn222_vuv","Fprompt_vs_photons_per_event_Rn222_vuv;Fprompt;Number of photons",100,0,1,300,0,300);

////////////
//Analysis//
////////////

    cout<<"New way of finding first photon thanks to Patrick..."<<endl;//debug time
    if(photon_cut == true){cout<<"photon cut is at "<<photon_cut_num<<" photons per events"<<endl;}
    int number_events_S = event_treeS->GetEntries();//Number of events for 5.6 MeV electrons
    int number_events_Rn222 = event_treeRn222->GetEntries();//Number of events for Rn222,for now it's all 1 million events

    int new_number_events_Rn222 = number_events_Rn222/2;

    //5.6 MeV electron with foils//
    vector<double> earliest_time_both_S(number_events_S,1000000);
    vector<double> fprompt_S(number_events_S,0);
    vector<int> quick_photon_S(number_events_S,0);
    vector<int> all_photon_S(number_events_S,0);

    //5.6 MeV electron without foils//
    vector<double> earliest_time_vuv_S(number_events_S,1000000);
    vector<double> fprompt_S_vuv(number_events_S,0);
    vector<int> quick_photon_S_vuv(number_events_S,0);
    vector<int> all_photon_S_vuv(number_events_S,0);

    //Rn222 with foils//
    vector<double> earliest_time_both_Rn222(new_number_events_Rn222,1000000);
    vector<double> fprompt_Rn222(new_number_events_Rn222,0);
    vector<int> quick_photon_Rn222(new_number_events_Rn222,0);
    vector<int> all_photon_Rn222(new_number_events_Rn222,0);

    //Rn222 without foils//
    vector<double> earliest_time_vuv_Rn222(new_number_events_Rn222,1000000);
    vector<double> fprompt_Rn222_vuv(new_number_events_Rn222,0);
    vector<int> quick_photon_Rn222_vuv(new_number_events_Rn222,0);
    vector<int> all_photon_Rn222_vuv(new_number_events_Rn222,0);
    
    //Finding first photons for 5.6 MeV electron with foils//
    int size_data_tree_S = data_treeS->GetEntries();
    for (int i = 0; i < size_data_tree_S; i++) {
        data_treeS->GetEntry(i);
        all_photon_S[data_eventS] ++;
        if (data_timeS < earliest_time_both_S[data_eventS]){
            earliest_time_both_S[data_eventS] = data_timeS;
        }
    }

    //Finding first photons for 5.6 MeV electron without foils//
    int size_data_tree_S_vuv = data_treeS_vuv->GetEntries();
    for (int i = 0; i < size_data_tree_S_vuv; i++) {
        data_treeS_vuv->GetEntry(i);
        all_photon_S_vuv[data_eventS_vuv] ++;
        if (data_timeS_vuv < earliest_time_vuv_S[data_eventS_vuv]){
            earliest_time_vuv_S[data_eventS_vuv] = data_timeS_vuv;
        }
    }

    //Finding first photon for Rn222 alphas with foils//
    int size_data_tree_Rn222 = data_treeRn222->GetEntries();
    for (int i = 0; i < size_data_tree_Rn222; i++) {
        data_treeRn222->GetEntry(i);
	if (data_eventRn222 % 2 == 0) {//even
	    new_data_eventRn222 = data_eventRn222/2;
	}
	if (data_eventRn222 % 2 == 1) {//odd
            new_data_eventRn222 = (data_eventRn222-1)/2;
        }

        all_photon_Rn222[new_data_eventRn222] ++;
        if (data_timeRn222 < earliest_time_both_Rn222[new_data_eventRn222]){
            earliest_time_both_Rn222[new_data_eventRn222] = data_timeRn222;
        }
    }

    //Finding first photon for Rn222 alphas without foils//
    int size_data_tree_Rn222_vuv = data_treeRn222_vuv->GetEntries();
    for (int i = 0; i < size_data_tree_Rn222_vuv; i++) {
        data_treeRn222_vuv->GetEntry(i);
        if (data_eventRn222_vuv % 2 == 0) {//even
            new_data_eventRn222_vuv = data_eventRn222_vuv/2;
        }
        if (data_eventRn222_vuv % 2 == 1) {//odd
            new_data_eventRn222_vuv = (data_eventRn222_vuv-1)/2;
        }

        all_photon_Rn222_vuv[new_data_eventRn222_vuv] ++;
        if (data_timeRn222_vuv < earliest_time_vuv_Rn222[new_data_eventRn222_vuv]){
            earliest_time_vuv_Rn222[new_data_eventRn222_vuv] = data_timeRn222_vuv;
        }
    }
    cout<<"Got them all"<<endl;//debug time

    //#pragma omp parallel for//Wang shi bei ding zhong yuan ri, Jia ji wu wang gao nai weng
    //Pulse shape for 5.6 MeV electron with foils//
    if(photon_cut == true){

      pulse_shape_S->SetTitle("pulse_shape_S_photon_cut");
      pulse_shape_S_vuv->SetTitle("pulse_shape_S_vuv_photon_cut");
      pulse_shape_Rn222->SetTitle("pulse_shape_Rn222_photon_cut");
      pulse_shape_Rn222_vuv->SetTitle("pulse_shape_Rn222_vuv_photon_cut");

      cout<<"Finding pulse shape with photon cut..."<<endl;
      for(int i = 0; i < size_data_tree_S; i++){
         data_treeS->GetEntry(i);
         if(all_photon_S[data_eventS] > photon_cut_num){
           pulse_shape_S->Fill(data_timeS - earliest_time_both_S[data_eventS]);
         }
      }
      //Pulse shape for 5.6 MeV electron without foils//
      for(int i = 0; i < size_data_tree_S_vuv; i++){//Let's try 6000 photons first

         data_treeS_vuv->GetEntry(i);
         if(all_photon_S_vuv[data_eventS_vuv] > photon_cut_num){
           pulse_shape_S_vuv->Fill(data_timeS_vuv - earliest_time_vuv_S[data_eventS_vuv]);
         }
      }
      //Pulse shape for Rn222 alphas with foils//
      for(int i = 0; i < size_data_tree_Rn222; i++){//Let's try 6000 photons first

         data_treeRn222->GetEntry(i);
	 if (data_eventRn222 % 2 == 0) {//even
            new_data_eventRn222 = data_eventRn222/2;
         }
         if (data_eventRn222 % 2 == 1) {//odd
             new_data_eventRn222 = (data_eventRn222-1)/2;
         }

         if(all_photon_Rn222[new_data_eventRn222] > photon_cut_num){
           pulse_shape_Rn222->Fill(data_timeRn222 - earliest_time_both_Rn222[new_data_eventRn222]);
         }
      }
      //Pulse shape for Rn222 alphas without foils//
      for(int i = 0; i < size_data_tree_Rn222_vuv; i++){//Let's try 6000 photons first

         data_treeRn222_vuv->GetEntry(i);
	 if (data_eventRn222_vuv % 2 == 0) {//even
            new_data_eventRn222_vuv = data_eventRn222_vuv/2;
         }
         if (data_eventRn222_vuv % 2 == 1) {//odd
             new_data_eventRn222_vuv = (data_eventRn222_vuv-1)/2;
         }

         if(all_photon_Rn222_vuv[new_data_eventRn222_vuv] > photon_cut_num){
           pulse_shape_Rn222_vuv->Fill(data_timeRn222_vuv - earliest_time_vuv_Rn222[new_data_eventRn222_vuv]);
         }
      }
    }
    if(photon_cut == false){
      cout<<"Finding pulse shape without photon cut..."<<endl;
      for(int i = 0; i < size_data_tree_S; i++){
         data_treeS->GetEntry(i);
         pulse_shape_S->Fill(data_timeS - earliest_time_both_S[data_eventS]);
      }
      for(int i = 0; i < size_data_tree_S_vuv; i++){//Let's try 6000 photons first
         data_treeS_vuv->GetEntry(i);
         pulse_shape_S_vuv->Fill(data_timeS_vuv - earliest_time_vuv_S[data_eventS_vuv]);
      }
      for(int i = 0; i < size_data_tree_Rn222; i++){//Let's try 6000 photons first
         data_treeRn222->GetEntry(i);
	 if (data_eventRn222 % 2 == 0) {//even
            new_data_eventRn222 = data_eventRn222/2;
         }
         if (data_eventRn222 % 2 == 1) {//odd
            new_data_eventRn222 = (data_eventRn222-1)/2;
         }

         pulse_shape_Rn222->Fill(data_timeRn222 - earliest_time_both_Rn222[new_data_eventRn222]);
      }
      for(int i = 0; i < size_data_tree_Rn222_vuv; i++){//Let's try 6000 photons first
         data_treeRn222_vuv->GetEntry(i);
	 if (data_eventRn222_vuv % 2 == 0) {//even
            new_data_eventRn222_vuv = data_eventRn222_vuv/2;
         }
         if (data_eventRn222_vuv % 2 == 1) {//odd
             new_data_eventRn222_vuv = (data_eventRn222_vuv-1)/2;
         }

         pulse_shape_Rn222_vuv->Fill(data_timeRn222_vuv - earliest_time_vuv_Rn222[new_data_eventRn222_vuv]);
      }
    }

    //Running integral//
    cout<<"Looking at running integral..."<<endl;
    double diff = 0;
    double diff_vuv = 0;
    double biggest_diff = 0;
    double biggest_diff_vuv = 0;
    int max_i = 0;
    int max_i_vuv = 0;
    double biggest_diff_time = 0;
    double biggest_diff_time_vuv = 0;

    double running_integral_S[NUM_BINS];
    double running_integral_Rn222[NUM_BINS];
    double running_integral_S_vuv[NUM_BINS];
    double running_integral_Rn222_vuv[NUM_BINS];
/*
    double running_inte_S = pulse_shape_S->Integral();
    double running_inte_S_vuv = pulse_shape_S_vuv->Integral();
    double running_inte_Rn222 = pulse_shape_Rn222->Integral();
    double running_inte_Rn222_vuv = pulse_shape_Rn222_vuv->Integral();
*/
    for(int i = 0;i<NUM_BINS;i++){

      running_integral_S[i] = (pulse_shape_S->Integral(0,i)*1.0)/(pulse_shape_S->Integral());
      running_integral_Rn222[i] = (pulse_shape_Rn222->Integral(0,i)*1.0)/(pulse_shape_Rn222->Integral());
      running_integral_S_vuv[i] = (pulse_shape_S_vuv->Integral(0,i)*1.0)/(pulse_shape_S_vuv->Integral());
      running_integral_Rn222_vuv[i] = (pulse_shape_Rn222_vuv->Integral(0,i)*1.0)/(pulse_shape_Rn222_vuv->Integral());

      diff = running_integral_Rn222[i] - running_integral_S[i];//This "abs" may need header file
      diff_vuv = running_integral_Rn222_vuv[i] - running_integral_S_vuv[i];

      if(diff>biggest_diff){
        biggest_diff = diff;
        max_i = i;
      }
      if(diff_vuv>biggest_diff_vuv){
        biggest_diff_vuv = diff_vuv;
        max_i_vuv = i;
      }
    }
    biggest_diff_time = max_i*10000.0/NUM_BINS;//[ns]
    biggest_diff_time_vuv = max_i_vuv*10000.0/NUM_BINS;//Now we have the biggest difference
    cout<<"Biggest difference with foils is at "<<biggest_diff_time<<" ns"<<endl;//debug,for 5 ns resolution, it should be 290ns
    cout<<"Biggest difference without foils is at "<<biggest_diff_time_vuv<<" ns"<<endl;//debug

    cout<<"Calculating Fprompt..."<<endl;//debug
    //Fprompt for 5.6 MeV electron with foils//
    for(int i = 0; i < size_data_tree_S; i++){
      data_treeS->GetEntry(i);
      if((data_timeS - earliest_time_both_S[data_eventS])<biggest_diff_time*0.001){
        quick_photon_S[data_eventS] ++;
      }
    }
    //Fprompt for 5.6 MeV electron without foils//
    for(int i = 0; i < size_data_tree_S_vuv; i++){
      data_treeS_vuv->GetEntry(i);
      if((data_timeS_vuv - earliest_time_vuv_S[data_eventS_vuv])<biggest_diff_time_vuv*0.001){
        quick_photon_S_vuv[data_eventS_vuv] ++;
      }
    }
    //Fprompt for Rn222 with foils//
    for(int i = 0; i < size_data_tree_Rn222; i++){
      data_treeRn222->GetEntry(i);
      if (data_eventRn222 % 2 == 0) {//even
            new_data_eventRn222 = data_eventRn222/2;
         }
         if (data_eventRn222 % 2 == 1) {//odd
            new_data_eventRn222 = (data_eventRn222-1)/2;
         }

      if((data_timeRn222 - earliest_time_both_Rn222[new_data_eventRn222])<biggest_diff_time*0.001){
        quick_photon_Rn222[new_data_eventRn222] ++;
      }
    }
    std::cout << "size of quick photon Rn222 is " << quick_photon_Rn222.size() << endl;
    //Fprompt for Rn222 without foils//
    for(int i = 0; i < size_data_tree_Rn222_vuv; i++){
      data_treeRn222_vuv->GetEntry(i);
      if (data_eventRn222_vuv % 2 == 0) {//even
          new_data_eventRn222_vuv = data_eventRn222_vuv/2;
          }
      if (data_eventRn222_vuv % 2 == 1) {//odd
          new_data_eventRn222_vuv = (data_eventRn222_vuv-1)/2;
          }

      if((data_timeRn222_vuv - earliest_time_vuv_Rn222[new_data_eventRn222_vuv])<biggest_diff_time_vuv*0.001){
        quick_photon_Rn222_vuv[new_data_eventRn222_vuv] ++;
      }
    }

    if(photon_cut == false){
      cout<<"Filling histograms for 5.6 MeV electrons..."<<endl;
      for(int i = 0;i<event_treeS->GetEntries();i++){
        event_treeS->GetEntry(i);
        if(all_photon_S[i] == 0){
          fprompt_S[i] = 0;
        }
        else{
          fprompt_S[i] = quick_photon_S[i]*1.0/all_photon_S[i];
        }
        fprompt_S_h->Fill(fprompt_S[i]);
        fprompt_vs_energy_S->Fill(fprompt_S[i],event_ES);
        fprompt_vs_photons_per_event_S->Fill(fprompt_S[i],all_photon_S[i]);

        if(all_photon_S_vuv[i] == 0){
          fprompt_S_vuv[i] = 0;
        }
        else{
          fprompt_S_vuv[i] = quick_photon_S_vuv[i]*1.0/all_photon_S_vuv[i];
        }
        fprompt_S_vuv_h->Fill(fprompt_S_vuv[i]);
        fprompt_vs_energy_S_vuv->Fill(fprompt_S_vuv[i],event_ES);
        fprompt_vs_photons_per_event_S_vuv->Fill(fprompt_S_vuv[i],all_photon_S_vuv[i]);
      }
      cout<<"Filling histograms for Rn222..."<<endl;
      for(int i = 0;i<event_treeRn222->GetEntries();i++){
        event_treeRn222->GetEntry(i);
	if (i % 2 == 0) {//even
            new_data_eventRn222 = i/2;
        }
        if (i % 2 == 1) {//odd
            new_data_eventRn222 = (i-1)/2;
        }

        if(all_photon_Rn222[new_data_eventRn222]==0){
          fprompt_Rn222[new_data_eventRn222]=0;
        }
        else{
          fprompt_Rn222[new_data_eventRn222] = quick_photon_Rn222[new_data_eventRn222]*1.0/all_photon_Rn222[new_data_eventRn222];
        }
	if (i % 2 == 1) {//odd
            fprompt_Rn222_h->Fill(fprompt_Rn222[new_data_eventRn222]);
            fprompt_vs_energy_Rn222->Fill(fprompt_Rn222[new_data_eventRn222],event_ERn222);
            fprompt_vs_photons_per_event_Rn222->Fill(fprompt_Rn222[new_data_eventRn222],all_photon_Rn222[new_data_eventRn222]);
        }
	/*
        fprompt_Rn222_h->Fill(fprompt_Rn222[new_data_eventRn222]);
        fprompt_vs_energy_Rn222->Fill(fprompt_Rn222[new_data_eventRn222],event_ERn222);
        fprompt_vs_photons_per_event_Rn222->Fill(fprompt_Rn222[new_data_eventRn222],all_photon_Rn222[new_data_eventRn222]);
*/
	if (i % 2 == 0) {//even
          new_data_eventRn222_vuv = i/2;
          }
        if (i % 2 == 1) {//odd
          new_data_eventRn222_vuv = (i-1)/2;
          }

        if(all_photon_Rn222_vuv[new_data_eventRn222_vuv]==0){
          fprompt_Rn222_vuv[new_data_eventRn222_vuv]=0;
        }
        else{
          fprompt_Rn222_vuv[new_data_eventRn222_vuv] = quick_photon_Rn222[new_data_eventRn222_vuv]*1.0/all_photon_Rn222[new_data_eventRn222_vuv];
        }
        if (i % 2 == 1) {//odd
            //fprompt_Rn222_vuv[new_data_eventRn222_vuv] = quick_photon_Rn222_vuv[new_data_eventRn222_vuv]*1.0/all_photon_Rn222_vuv[new_data_eventRn222_vuv];
            fprompt_Rn222_vuv_h->Fill(fprompt_Rn222_vuv[new_data_eventRn222_vuv]);
            fprompt_vs_energy_Rn222_vuv->Fill(fprompt_Rn222_vuv[new_data_eventRn222_vuv],event_ERn222);
            fprompt_vs_photons_per_event_Rn222_vuv->Fill(fprompt_Rn222_vuv[new_data_eventRn222_vuv],all_photon_Rn222_vuv[new_data_eventRn222_vuv]);
	}
      }
    }

    if(photon_cut == true){ //Photon number cut

      fprompt_S_h->SetTitle("Fprompt_S_photon_cut");
      fprompt_vs_energy_S->SetTitle("Fprompt_vs_energy_S_photon_cut");
      fprompt_vs_photons_per_event_S->SetTitle("Fprompt_vs_photons_per_event_S_photon_cut");
      fprompt_S_vuv_h->SetTitle("Fprompt_S_vuv_photon_cut");
      fprompt_vs_energy_S_vuv->SetTitle("Fprompt_vs_energy_S_vuv_photon_cut");
      fprompt_vs_photons_per_event_S_vuv->SetTitle("Fprompt_vs_photons_per_event_S_vuv_photon_cut");
      fprompt_Rn222_h->SetTitle("Fprompt_Rn222_photon_cut");
      fprompt_vs_energy_Rn222->SetTitle("Fprompt_vs_energy_Rn222_photon_cut");
      fprompt_vs_photons_per_event_Rn222->SetTitle("Fprompt_vs_photons_per_event_Rn222_photon_cut");
      fprompt_Rn222_vuv_h->SetTitle("Fprompt_Rn222_vuv_photon_cut");
      fprompt_vs_energy_Rn222_vuv->SetTitle("Fprompt_vs_energy_Rn222_vuv_photon_cut");
      fprompt_vs_photons_per_event_Rn222_vuv->SetTitle("Fprompt_vs_photons_per_event_Rn222_vuv_photon_cut");

      cout<<"Filling histograms for electrons with photon number cut..."<<endl;
      for(int i = 0;i<event_treeS->GetEntries();i++){
        if(all_photon_S[i] > photon_cut_num){
          event_treeS->GetEntry(i);
          fprompt_S[i] = quick_photon_S[i]*1.0/all_photon_S[i];
          fprompt_S_h->Fill(fprompt_S[i]);
          fprompt_vs_energy_S->Fill(fprompt_S[i],event_ES);
          fprompt_vs_photons_per_event_S->Fill(fprompt_S[i],all_photon_S[i]);
        }

        if(all_photon_S_vuv[i] > photon_cut_num){
          fprompt_S_vuv[i] = quick_photon_S_vuv[i]*1.0/all_photon_S_vuv[i];
          fprompt_S_vuv_h->Fill(fprompt_S_vuv[i]);
          fprompt_vs_energy_S_vuv->Fill(fprompt_S_vuv[i],event_ES);
          fprompt_vs_photons_per_event_S_vuv->Fill(fprompt_S_vuv[i],all_photon_S_vuv[i]);
        }
      }
      cout<<"Filling histograms for Rn222 with photon number cut..."<<endl;
      for(int i = 0;i<event_treeRn222->GetEntries();i++){
        if(all_photon_Rn222[i] > photon_cut_num){
          event_treeRn222->GetEntry(i);
          fprompt_Rn222[i] = quick_photon_Rn222[i]*1.0/all_photon_Rn222[i];
          fprompt_Rn222_h->Fill(fprompt_Rn222[i]);
          fprompt_vs_energy_Rn222->Fill(fprompt_Rn222[i],event_ERn222);
          fprompt_vs_photons_per_event_Rn222->Fill(fprompt_Rn222[i],all_photon_Rn222[i]);
        }

        if(all_photon_Rn222_vuv[i] > photon_cut_num){
          fprompt_Rn222_vuv[i] = quick_photon_Rn222_vuv[i]*1.0/all_photon_Rn222_vuv[i];
          fprompt_Rn222_vuv_h->Fill(fprompt_Rn222_vuv[i]);
          fprompt_vs_energy_Rn222_vuv->Fill(fprompt_Rn222_vuv[i],event_ERn222);
          fprompt_vs_photons_per_event_Rn222_vuv->Fill(fprompt_Rn222_vuv[i],all_photon_Rn222_vuv[i]);
        }
      }
    }

//////////
//Output//
//////////

    std::string res;//Defining may take time
    std::stringstream ss;
    ss<<reso;ss>>res;
    std::string res1;
    std::stringstream ss1;
    ss1<<photon_cut_num;ss1>>res1;
    std::string out_front("output/pulse_shape_merged_a-gamma_test_");
    std::string out_back("ns_slow_alpha.root");
    std::string out_mid("_photon_cut_");
    if(photon_cut == true){out_front+=res1;out_front+=out_mid;}
    out_front+=res;out_front+=out_back;

    //TFile *fO = new TFile("pulse_shape_560KeV_Rn222_foils_no-foils_500ns_slow_alpha.root","RECREATE");
    TFile *fO = new TFile(out_front.c_str(),"RECREATE");

    pulse_shape_S->Write();
    pulse_shape_S_vuv->Write();
    pulse_shape_Rn222->Write();
    pulse_shape_Rn222_vuv->Write();
    fprompt_S_h->Write();
    fprompt_S_vuv_h->Write();
    fprompt_Rn222_h->Write();
    fprompt_Rn222_vuv_h->Write();
    fprompt_vs_energy_S->Write();
    fprompt_vs_energy_S_vuv->Write();
    fprompt_vs_energy_Rn222->Write();
    fprompt_vs_energy_Rn222_vuv->Write();
    fprompt_vs_photons_per_event_S->Write();
    fprompt_vs_photons_per_event_S_vuv->Write();
    fprompt_vs_photons_per_event_Rn222->Write();
    fprompt_vs_photons_per_event_Rn222_vuv->Write();
    fO->Close();

////////
//Draw//
////////

    return 0;
}
