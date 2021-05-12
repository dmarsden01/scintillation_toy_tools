#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include "TGraph.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TStyle.h"

int Filter_plotting(){

    //Setting up the canvas
    TCanvas *c = new TCanvas("c","Number of photons against x position of event",1000,800);//Final numbers are the number of pixels of the canvas
    c->SetGrid();
 
    //Which output file to read in
    TFile *f_Solar = new TFile("../new_data/8B_10000_uniform.root"); 
    //TFile *f_Solar = new TFile("../new_data/a_gamma_10000_uniform_full_scint.root");
    
    //This makes pointers to the TTrees
    TTree* event_tree_Solar = (TTree*)f_Solar->Get("event_tree"); //TTree for events 
    TTree* data_tree_Solar = (TTree*)f_Solar->Get("data_tree"); //TTree for individual photons detected
    TTree* data_tree_vuv_Solar = (TTree*)f_Solar->Get("data_tree_vuv"); //TTree for vuv photons
    TTree* data_tree_vis_Solar = (TTree*)f_Solar->Get("data_tree_vis"); //TTree for visible photons

    //-----------How many events to run over------------------//
    int number_events = event_tree_Solar->GetEntries();
    std::cout << "Number of events to run over is " << number_events << endl;

    //Storing the event a photon is associated with and the total number of detected photons
    int data_event, data_event_vuv, data_event_vis, pmt_number, pmt_number_vis, data_time;
    data_tree_Solar->SetBranchAddress("data_event", &data_event);//This stores the number of the event from which the photon was produced
    data_tree_Solar->SetBranchAddress("data_time", &data_time);
    data_tree_vuv_Solar->SetBranchAddress("data_event_vuv", &data_event_vuv);
    data_tree_vis_Solar->SetBranchAddress("data_event_vis", &data_event_vis);
    data_tree_Solar->SetBranchAddress("data_pmt", &pmt_number); //Gets the pmt each photon is detected at
    data_tree_vis_Solar->SetBranchAddress("data_pmt_vis", &pmt_number_vis); //Gets the pmt visible photons are detected at
    
    //Storing the number of photon hits on the two populations of optical detectors
    std::vector<int> coated_hits(number_events,0);
    std::vector<int> uncoated_hits(number_events,0);
    std::vector<int> true_vuv(number_events,0);
    std::vector<int> true_vis(number_events,0);    

    //Sets how many pmts per region of pmts
    int pmts_per_region = 6;
    int pmt_runs = ceil(480/pmts_per_region);
    std::vector<int> cell_number(pmt_runs,0);
    std::vector<double> frac_of_cell(number_events,0);
    int event_counter = 0;
    int photon_counter = 0;
    int cell_counter;

    //Fills a vector with the true number of vuv photon hits on all opdets    
    int size_data_tree_vuv = data_tree_vuv_Solar->GetEntries();
    std::cout << "Number of vuv photons is " << size_data_tree_vuv << endl;
        for (int i = 0; i < size_data_tree_vuv; i++) {
            data_tree_vuv_Solar->GetEntry(i);
            true_vuv[data_event_vuv]++;
        }

    //Fills a vector with the number of hits on coated opdets for each event
    int size_data_tree_total = data_tree_Solar->GetEntries();
        for (int i = 0; i < size_data_tree_total; i++) {
            data_tree_Solar->GetEntry(i);
	    if (pmt_number % 3 != 0) { //Here we set 2/3 of the optical detectors to be coated (i.e detect visible and vuv photons)
                coated_hits[data_event]++;
	    if (data_event == event_counter) {
		    
		    for (int n = 0; n < pmt_runs; n++) {
			if (pmt_number >= n*pmts_per_region && pmt_number < (n+1)*pmts_per_region) {
			    cell_number[n]++;
			}
	    }
	    }
	    else {
		photon_counter = 0;
		cell_counter = pmt_runs-1;
		std::sort(cell_number.begin(), cell_number.end());    
		while (photon_counter <= 0.8*coated_hits[event_counter]) {
		    photon_counter += cell_number[cell_counter];
		    cell_counter--;
		}
		frac_of_cell[event_counter] = (double(pmt_runs) - cell_counter)/pmt_runs;
	        event_counter = data_event;	
		std::fill(cell_number.begin(), cell_number.end(), 0);
		for (int n = 0; n < pmt_runs; n++) {
                    if (pmt_number >= n*pmts_per_region && pmt_number < (n+1)*pmts_per_region) {
                            cell_number[n]++;
                        }
	        }
        }
	}
	}
    //Fills a vector with the number of visible hits on uncoated opdets for each event
    int size_data_tree_vis = data_tree_vis_Solar->GetEntries();
    std::cout << "Number of visible photons is " << size_data_tree_vis << endl;
        for (int i = 0; i < size_data_tree_vis; i++) {
            data_tree_vis_Solar->GetEntry(i);
	    true_vis[data_event_vis]++;
            if (pmt_number_vis % 3 == 0) {
	        uncoated_hits[data_event_vis]++;
	    }
        }

    //Storing the positions and Energy of the events
    double event_x_Solar, event_y_Solar, event_z_Solar, event_E;
    event_tree_Solar->SetBranchAddress("event_x_pos", &event_x_Solar);
    event_tree_Solar->SetBranchAddress("event_y_pos", &event_y_Solar);
    event_tree_Solar->SetBranchAddress("event_z_pos", &event_z_Solar);
    event_tree_Solar->SetBranchAddress("event_E", &event_E);

    std::vector<double> x_pos(number_events,0);
    std::vector<double> y_pos(number_events,0);
    std::vector<double> Energy(number_events,0);    
    std::vector<double> ratio_coated_uncoated(number_events,0);
    std::vector<double> frac_of_pmts;
    std::vector<double> closest_y;
    std::vector<double> E_cut_vec, x_cut_vec, uncoated_x_cut_vec, coated_x_cut_vec, several_x_cut_vec, rat_cut_vec; //Stores the energy of the remaining events after a cut

    std::vector<int> num_events_per_bin (50); //Stores the number of events in each histogram bin
    int bin_number = 0; //Stores the bin number the event is associated with (by using the position)
    for(int i = 0; i < number_events; i++) {
        event_tree_Solar->GetEntry(i);
	bin_number = floor( event_x_Solar/7.0);
	num_events_per_bin[bin_number] ++;
    }

    //A 1D histogram for the energy of events
    //TH1D *h_uncut = new TH1D("histo", "", 100, 0, 25);
    //TH1D *h_cut_2 = new TH1D("histo", "", 100, 0, 15);
    //TH1D *h_uncut = new TH1D("histo", "", 100, 0, 5);
    //A 1D histogram for the number of uncoated photons detected
    //TH1D *h_vis = new TH1D("histo_vis", "", 50, 0, 350);
    //A 2D histogram for the number of uncoated photons detected
    //TH2D *h_vis_2D = new TH2D("2D_histo_vis", "10000 solar events", 125, 0, 350,100, 0, 100);
    //A 2D histogram for the number of coated minus 2*uncoated photons
    //TH2D *h_vuv_2D = new TH2D("2D_histo_vuv", "10000 solar events", 125, 0, 350,100, 0, 250);
    //A 1D histogram for the number of coated photons detected
    //TH1D *h_uncoated = new TH1D("histo_coated", "", 50, 0, 350);
    //A 1D histogram for approximating the number of vuv photons detected
    //TH1D *h_vuv = new TH1D("histo_vuv", "", 50, 0, 350);
    //A 1D histogram for the number of photons detected at each position
    //TH1D *h_total_uncut = new TH1D("histo_uncut", "", 100, 0, 350);
    //A 1D histogram for the number of uncoated photons from each event
    //TH1D *h_uncoated_num = new TH1D("histo_uncoated_num", "", 30, 0, 30);

    TH2D *h2 = new TH2D("", "", 125, 0, 370,125, 0, 0.9);    

    //TH2D *h_region = new TH2D("", "", 125, 0, 370, 100, 0, 0.7);

    int num_over_thresh = 0; //This is a counter for the number of events with a ratio over a threshold
    int num_passed_distance = 0; //This is a counter for the number of events with ratio over the threshold which are below a certain distance from the anode
    double ratio_vis_vuv, ratio_vuv_vis;
    double new_ratio;
    double threshold = 0.3; //This is the threshold ratio above which events will be cut out
    int num_plotted = 0;//This is to store the number of events which pass the x threshold
    double x_thresh = 0.0;
    double x_thresh_2 = 65.0;
    ////----E cut-----/////
    double E_thresh = 0.0; //This is the energy below which events can be cut out
    double ratio_thresh = 0.25;
    int photon_thresh = 0; //This is the number of photons per event below which events will be cut
    int photon_thresh_2 = 0;
    int uncoated_photon_thresh = 0;
    int coated_photon_thresh = 0;
    int difference_thresh = 0;
    double region_frac_thresh = 0.15;
    double frac_above;
    double frac_passed;
    double mean_ratio;
    double total_ratio;
    double y_distance;

    double scaling;

    int E_cut_2_num = 0;

    //Loop over each event
    for(int i = 0; i < number_events; i++) {
        event_tree_Solar->GetEntry(i);
	x_pos[i] = event_x_Solar;//Array filled with x positions of events
        y_pos[i] = event_y_Solar;//Array filled with y positions of events
        Energy[i] = event_E; //This stores the energy of each event    
	
	//h_uncoated_num->Fill(coated_hits[i]);	
       	//This fills a 1D histogram of the energy of events
        //h_uncut->Fill(event_E);
	//This fills a 1D histogram with the position of events
	//h_total_uncut->Fill(event_x_Solar);
	//Fills a 2D histogram with the position and number of photons per event
	//h_vis_2D->Fill(event_x_Solar, uncoated_hits[i]);
	//Fills a 2D histogram with the position and number of photons per event
	//h_vuv_2D->Fill(event_x_Solar, coated_hits[i]-2*uncoated_hits[i]);
	//h_region->Fill(event_x_Solar, frac_of_cell[i]);

	ratio_vis_vuv = ((double) true_vis[i])/true_vuv[i]; //Total number of visible photons detected over total number of vuv photons detected
        ratio_vuv_vis = ((double) true_vuv[i])/true_vis[i];
	new_ratio = ((double) uncoated_hits[i])/coated_hits[i]; //Photons detected by uncoated pmts over photons detected by coated pmts
	ratio_coated_uncoated[i] = new_ratio;
	total_ratio += ratio_vis_vuv;

        //h2->Fill(event_x_Solar, new_ratio);

	//Making photon cut on events
	//Cut on uncoated
        //if (frac_of_cell[i] <= region_frac_thresh) {
	if (new_ratio < ratio_thresh || frac_of_cell[i] <= region_frac_thresh) {
            //h_region->Fill(event_x_Solar, frac_of_cell[i]);
	    h2->Fill(event_x_Solar, new_ratio);
	    //if (new_ratio < ratio_thresh) {
	    //if (frac_of_cell[i] <= region_frac_thresh) {
                uncoated_x_cut_vec.push_back (event_x_Solar);
	    //}
        }
	//Cut on coated
        if (coated_hits[i] >= coated_photon_thresh) {
            coated_x_cut_vec.push_back (event_x_Solar);
        }
	//Chain of cuts 
        if (coated_hits[i]+uncoated_hits[i] >= photon_thresh) {
	    //rat_cut_vec.push_back (new_ratio);
	    E_cut_vec.push_back (event_E);
	    //if (uncoated_hits[i] >= uncoated_photon_thresh) {
	        //if (coated_hits[i]-2*uncoated_hits[i] >= difference_thresh) {
                    //if (coated_hits[i]+uncoated_hits[i] >= photon_thresh_2) {
		    if (new_ratio < ratio_thresh) {
		        several_x_cut_vec.push_back (event_x_Solar);
		        rat_cut_vec.push_back (new_ratio);
			//h_cut_2->Fill(event_E);
			E_cut_2_num++;
		        //h2->Fill(event_x_Solar, new_ratio);
		    //}
		    //}
	        }
	    //}
	}
    

	//Making a cut on the ratio
	if (new_ratio > threshold) {
	    num_over_thresh++;
	    if (x_pos[i] < x_thresh_2) {
	        num_passed_distance++;
	    }
        }

        if (x_pos[i]>x_thresh && x_pos[i]<x_thresh_2) {
            //The following two statements take the modulus of y[i], as the pmt positions are symetric about 0, and then outputs the y-distance from an event to the closest pmt. 
            if (y_pos[i] < 0) {
                y_pos[i] = -1.0*y_pos[i];
            }
            y_distance = std::fmod(y_pos[i],62.0)-31;
                if (y_distance < 0) {
                    y_distance = -1.0*y_distance;
                }
            closest_y.push_back (y_distance);
	}

    }//End of loop over events

    //----PLOTTING----//

    //This is for a coloured 2D histogram of the events
    //TH2D *h2 = new TH2D("", "", 125, 0, 350,125, 0, 0.9);

    //This is for a 1D histogram of events
    //TH1D *h_cut = new TH1D("histo_cut", "", 100, 0, 15);//For solar events
    //TH1D *h_cut = new TH1D("histo_cut", "", 100, 0, 5);
    //TH1D *h_total = new TH1D("histo_uncoated_cut", "", 100, 0, 350);
    //TH1D *h_coated = new TH1D("histo_coated_cut", "", 100, 0, 350);
    //TH1D *h_difference = new TH1D("histo_difference_cut", "", 100, 0, 350);
    //TH1D *h_several = new TH1D("histo_several_cut", "", 100, 0, 350);

    //This is for making an energy/number of photons cut
    auto size_E = E_cut_vec.size();
    auto size_total = uncoated_x_cut_vec.size();
    auto size_coated = coated_x_cut_vec.size();
    auto size_several = several_x_cut_vec.size();
    double fraction_plotted = (double)size_several/number_events;
    double fraction_plotted_E_2 = (double)E_cut_2_num/number_events;
    double fraction_plotted1 = (double)size_total/number_events;
    double fraction_plotted2 = (double)size_coated/number_events;
    double fraction_plotted_several = (double)size_several/number_events;
    double perc_rounded = roundf(fraction_plotted*1000)/10;
    double perc_rounded_E_2 = roundf(fraction_plotted_E_2*1000)/10;
    double perc_rounded1 = roundf(fraction_plotted1*1000)/10;
    double perc_rounded2 = roundf(fraction_plotted2*1000)/10;
    double perc_rounded_several = roundf(fraction_plotted_several*1000)/10;

/*    for(int i = 0; i < size_total; i++) {
        h_total->Fill(uncoated_x_cut_vec[i]);
    }

    for(int i = 0; i < size_coated; i++) {
        h_coated->Fill(coated_x_cut_vec[i]);
    }

    for(int i = 0; i < size_difference; i++) {
        h_difference->Fill(difference_x_cut_vec[i]);
    }
    
    for(int i = 0; i < size_several; i++) {
        h_several->Fill(several_x_cut_vec[i]);
	//h2->Fill(several_x_cut_vec[i], rat_cut_vec[i]);
	//h_cut->Fill(E_cut_vec[i]);
    }

    for(int i = 0; i < size_E; i++) {
	//E_cut[i] = E_cut_vec.at(i);
	//x_cut[i] = x_cut_vec.at(i);
	//rat_cut[i] = rat_cut_vec.at(i);
	//h2->Fill(x_cut[i], rat_cut[i]);
	h_cut->Fill(E_cut_vec[i]);
    }

    frac_above = ((double)num_over_thresh)/num_events;
    frac_passed = ((double)num_passed_distance)/num_over_thresh;
    mean_ratio = total_ratio/num_events;
*/
    //This plots a graph of the number of photons per event against the x position of the event
    //TGraph *g = new TGraph(size_E, x_cut, rat_cut);

    gStyle->SetOptStat(0);

    //Scatter plot for ratio against distance
    /*
    TGraph *g = new TGraph(num_events, x, rat_Solar_vis);
    g->SetTitle("Solar events with 2/3 detectors TPB coated;" "Distance of event from Anode (cm);" "Ratio of uncoated opdet photons to coated opdet photons;");
    g->SetMarkerSize(1.0);
    g->SetMarkerStyle(21);
    //g->GetXaxis()->SetLabelSize(0.03);
    //g->GetXaxis()->SetTitleSize(0.04);
    //g->GetXaxis()->SetTitleFont(132);
    //g->SetTitleOffset(1.2,"Y");
    g->Draw("AP");
    */

    //Scatter plot for y distance to closest opdet
    /*
    TGraph *g = new TGraph(num_plotted, y_plot_new, vuv_new_y);
    g->SetTitle("Solar events with photon cut of 100;" "Vertical distance of event from closest pmt (cm);" "Ratio of uncoated PMT photons to coated PMT photons;");
    g->SetMarkerSize(1.0);
    g->SetMarkerStyle(21);
    //g->GetXaxis()->SetLabelSize(0.03);
    //g->GetXaxis()->SetTitleSize(0.04);
    //g->GetXaxis()->SetTitleFont(132);
    //g->SetTitleOffset(1.2,"Y");
    g->Draw("AP");
    */

    //1D Energy histogram drawing
    
    ostringstream d;
    d << perc_rounded_several << "% remaining";
    ostringstream d1;
    d1 << "Combined cut 0.25 and 0.3: " << perc_rounded1 << "% remain" << endl;
    ostringstream d2;
    d2 << perc_rounded2 << "% remaining";
    ostringstream d_several;
    d_several << "15 estimated vuv: " <<  perc_rounded << "% remain" << endl;
    ostringstream d_several_E_2;
    d_several_E_2 << "7 uncoated photons: " << perc_rounded1 << "% remain" << endl;
/*
    //Energy 1D
    //h_uncut->SetTitle("10000 solar events with and without photon cuts");
    h_uncut->SetTitle("10000 alpha-gamma events");
    h_uncut->SetLabelSize(0.04,"xyz");
    h_uncut->SetTitleSize(0.045,"xyz"); // size of axis title font
    h_uncut->SetTitleFont(132,"xyz"); // font option
    h_uncut->SetLabelFont(132,"xyz");
    h_uncut->SetTitleOffset(0.95,"Y");
    h_uncut->GetXaxis()->SetTitle("Energy (MeV)");
    h_uncut->GetYaxis()->SetTitle("Number of events");
    h_uncut->SetLineColor( kBlack);
    //h_cut->SetLineColor( kBlue);
    //h_cut_2->SetLineColor(kRed);
    h_uncut->Draw();
    //h_cut->Draw("same");
    //h_cut_2->Draw("same");
    TLatex t_uncut(9,190,"#splitline{No photon cut}{}");
    //TLatex t_uncut(1.5,100,"#splitline{No photon cut}{}");
    t_uncut.SetTextSize(0.025);
    t_uncut.SetTextColor(kBlack);
    t_uncut.DrawClone();

    TLatex t_cut(10,180,"#splitline{50 total photons}{}");
    //TLatex t_cut(1.5,90,"#splitline{Several photon cuts}{}");
    t_cut.SetTextSize(0.025);
    t_cut.SetTextColor(kBlue);
    t_cut.DrawClone();

    TLatex t_cut_2(10,170,"#splitline{100 total photons}{}");
    t_cut_2.SetTextSize(0.025);
    t_cut_2.SetTextColor(kRed);
    t_cut_2.DrawClone();    
    
    TLatex t_portion(9,180,d_several.str().c_str());
    //TLatex t_portion(1.5,80,d_several.str().c_str());
    t_portion.SetTextSize(0.025);
    t_portion.SetTextColor(kBlue);
    t_portion.DrawClone();    

    TLatex t_portion_E_2(9,170,d_several_E_2.str().c_str());
    //TLatex t_portion(1.5,80,d_several_E_2.str().c_str());
    t_portion_E_2.SetTextSize(0.025);
    t_portion_E_2.SetTextColor(kRed);
    t_portion_E_2.DrawClone();

    h_total_uncut->SetLabelSize(0.04,"xyz");
    h_total_uncut->SetTitleSize(0.045,"xyz"); // size of axis title font
    h_total_uncut->SetTitleFont(132,"xyz"); // font option
    h_total_uncut->SetLabelFont(132,"xyz");
    h_total_uncut->SetTitleOffset(0.9,"Y");   

    h_total_uncut->SetTitle("");
    h_total_uncut->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_total_uncut->GetYaxis()->SetTitle("Number of events");
    h_total_uncut->SetLineColor( kBlack);
    h_total_uncut->SetLineWidth(1);
    h_total_uncut->Draw();
    
    //Uncoated cut
    
    h_total->SetLineColor( kBlue);    
    h_total->Draw("same");
    
    //Coated cut
    
    h_coated->SetLineColor( kRed);
    h_coated->Draw("same");
    */
    //Difference cut
/*
    h_difference->SetLineColor( kRed);
    h_difference->Draw("same");

    //Several cuts
    
    h_several->SetLineColor( kRed);
    h_several->Draw("same");

    TLatex t_uncut(100,45,"#splitline{No photon cut}{}");
    t_uncut.SetTextSize(0.035);
    t_uncut.SetTextColor(kBlack);
    t_uncut.DrawClone();

    //TLatex t_cut(275,60,"#splitline{Ratio and total photon cuts}{}");
    //Uncoated cut
    TLatex t_cut(100,40,d1.str().c_str());
    t_cut.SetTextSize(0.035);
    t_cut.SetTextColor(kBlue);
    t_cut.DrawClone();
    //Coated cut
    //TLatex t_portion(250,40,d2.str().c_str());
    //Difference cut
    //TLatex t_portion(250,40,diff.str().c_str());
    //Several cuts
    TLatex t_portion(100,30,d_several.str().c_str());
    t_portion.SetTextSize(0.035);
    t_portion.SetTextColor(kRed);
    t_portion.DrawClone();

    //1D Position histogram drawing
    
    h_vis->SetTitle("Photons detected by uncoated opdets for 10000 solar events");
    h_vis->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_vis->GetYaxis()->SetTitle("Number of entries divided by events per bin");
    //h_vis->SetLineColor( kBlue);
    h_vis->Draw();
    */
    //Distribution of coated/uncoated photons detected 
/*
    h_uncoated_num->SetTitle("Photons detected by coated opdets for 10000 40K Beta events");
    h_uncoated_num->GetXaxis()->SetTitle("Number of photons detected");
    h_uncoated_num->GetYaxis()->SetTitle("Number of entries");
    //h_vis->SetLineColor( kBlue);
    h_uncoated_num->Draw();
*/
    //1D for coated
/*
    h_coated->SetTitle("Photons detected by uncoated opdets for 10000 solar events");
    h_coated->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_coated->GetYaxis()->SetTitle("Number of entries divided by events per bin");
    //h_coated->SetLineColor( kBlue);
    h_coated->Draw();
*/
    //1D for vuv
/*
    h_vuv->Add(h_vis, -2);
    h_vuv->SetTitle("Estimate for vuv photons detected (coated - 2*uncoated)");
    h_vuv->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_vuv->GetYaxis()->SetTitle("Number of entries divided by events per bin");
    //h_vuv->SetLineColor( kBlue);
    h_vuv->Draw();

    //2D histogram drawing ratio and distance
*/    
    h2->SetLabelSize(0.04,"xyz");
    h2->SetTitleSize(0.045,"xyz"); // size of axis title font
    h2->SetTitleFont(132,"xyz"); // font option
    h2->SetLabelFont(132,"xyz");
    h2->SetTitleOffset(0.9,"Y");
    h2->GetXaxis()->SetTitle("Distance from anode (cm)");
    h2->GetYaxis()->SetTitle("Ratio of uncoated to coated photons");
    gStyle->SetOptStat(0);
    h2->Draw("COLZ");
    TLatex t_portion(25,0.65,d1.str().c_str());
    t_portion.SetTextSize(0.035);
    t_portion.DrawClone();
/*

    h_region->SetTitle("10000 8B events, 6 optical channels per region");
    h_region->SetLabelSize(0.04,"xyz");
    h_region->SetTitleSize(0.045,"xyz"); // size of axis title font
    h_region->SetTitleFont(132,"xyz"); // font option
    h_region->SetLabelFont(132,"xyz");
    h_region->SetTitleOffset(0.9,"Y");
    h_region->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_region->GetYaxis()->SetTitle("Fraction of regions for 80% of total light");
    gStyle->SetOptStat(0);
    h_region->Draw("COLZ");
    TLatex t_portion(25,0.5,d1.str().c_str());
    t_portion.SetTextSize(0.035);
    t_portion.DrawClone();

    //2D histogram of number of uncoated per event
    
    h_vis_2D->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_vis_2D->GetYaxis()->SetTitle("Number of photons detected by uncoated opdets");
    gStyle->SetOptStat(0);
    h_vis_2D->Draw("COLZ");
    */
    //TLatex t_portion(25,0.65,d_several.str().c_str());
    //t_portion.SetTextSize(0.035);
    //t_portion.DrawClone();
    //2D histogram of coated-2*uncoated
  /*  
    h_vuv_2D->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_vuv_2D->GetYaxis()->SetTitle("Estimated number of vuv photons detected");
    gStyle->SetOptStat(0);
    h_vuv_2D->Draw("COLZ");
*/
    //c->SaveAs("../plots/geometry/combined_OR_8B_10000_region_6_0.2_0.25.png"); //Saving the plot in the canvas

    return 0;

    }
