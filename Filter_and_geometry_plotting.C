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

int Filter_and_geometry_plotting(){

    //Setting up the canvas
    TCanvas *c = new TCanvas("c","Number of photons against x position of event",1000,800);//Final numbers are the number of pixels of the canvas
    c->SetGrid();
 
    //Which toy MC output file to read in
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
    int data_event, data_event_vuv, data_event_vis, pmt_number, pmt_number_vis;
    double data_time;
    data_tree_Solar->SetBranchAddress("data_event", &data_event);//This stores the number of the event from which the photon was produced
    data_tree_Solar->SetBranchAddress("data_time", &data_time);
    data_tree_vuv_Solar->SetBranchAddress("data_event_vuv", &data_event_vuv);
    data_tree_vis_Solar->SetBranchAddress("data_event_vis", &data_event_vis);
    data_tree_Solar->SetBranchAddress("data_pmt", &pmt_number); //Gets the pmt each photon is detected at
    data_tree_vis_Solar->SetBranchAddress("data_pmt_vis", &pmt_number_vis); //Gets the pmt visible photons are detected at
    
    //----------Filters----------//
    //Storing the number of photon hits on the two populations of opdets, i.e coated and uncoated
    std::vector<int> coated_hits(number_events,0);
    std::vector<int> uncoated_hits(number_events,0);
    std::vector<int> true_vuv(number_events,0);
    std::vector<int> true_vis(number_events,0);    

    //----------Geometry---------//
    //Splitting the opdets into "regions" this sets how many per region
    int pmts_per_region = 6;
    int pmt_runs = ceil(480/pmts_per_region);
    std::vector<int> cell_number(pmt_runs,0);
    std::vector<double> frac_of_cell(number_events,0); 
    int event_counter = 0;
    int photon_counter = 0;
    int cell_counter;
    double fraction_total = 0.8;

    //Fills a vector with the true number of vuv photon hits on all opdets    
    int size_data_tree_vuv = data_tree_vuv_Solar->GetEntries();
    std::cout << "Number of vuv photons is " << size_data_tree_vuv << endl;
        for (int i = 0; i < size_data_tree_vuv; i++) { //Loop through all VUV photons detected
            data_tree_vuv_Solar->GetEntry(i);
            true_vuv[data_event_vuv]++;
        }

    //Fills a vector with the number of hits on coated opdets for each event
    int size_data_tree_total = data_tree_Solar->GetEntries();
        for (int i = 0; i < size_data_tree_total; i++) { //Loop through all photons detected
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
		while (photon_counter <= fraction_total*coated_hits[event_counter]) { //Here we get the fraction of opdet regions needed to reach fraction of total light 
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
	} //End of loop over all photons
    //Fills a vector with the number of visible hits on uncoated opdets for each event
    int size_data_tree_vis = data_tree_vis_Solar->GetEntries();
    std::cout << "Number of visible photons is " << size_data_tree_vis << endl;
        for (int i = 0; i < size_data_tree_vis; i++) { //Loop through all visible photons detected
            data_tree_vis_Solar->GetEntry(i);
	    true_vis[data_event_vis]++;
            if (pmt_number_vis % 3 == 0) { //Here 1/3 are uncoated
	        uncoated_hits[data_event_vis]++;
	    }
        }

    //Variables for the positions and energies of the events
    double event_x_Solar, event_y_Solar, event_z_Solar, event_E;
    event_tree_Solar->SetBranchAddress("event_x_pos", &event_x_Solar);
    event_tree_Solar->SetBranchAddress("event_y_pos", &event_y_Solar);
    event_tree_Solar->SetBranchAddress("event_z_pos", &event_z_Solar);
    event_tree_Solar->SetBranchAddress("event_E", &event_E);

    std::vector<double> x_pos(number_events,0);
    std::vector<double> y_pos(number_events,0);
    std::vector<double> Energy(number_events,0);    
    std::vector<double> ratio_coated_uncoated(number_events,0);  //This stores the ratio of hits on coated to uncoated opdets
    std::vector<double> frac_of_pmts;
    std::vector<double> closest_y; //This is for a test to look at the y-dependence of light detected (found there was only a correlation for x < 50cm)
    std::vector<double> E_cut_vec, x_cut_vec, uncoated_x_cut_vec, coated_x_cut_vec, several_x_cut_vec, rat_cut_vec;

    std::vector<int> num_events_per_bin (50); //To make a histogram with 50 bins which displays the x-distribution of events
    int bin_number = 0; //Stores the bin number the event is associated with (by using the position)
    for(int i = 0; i < number_events; i++) {
        event_tree_Solar->GetEntry(i);
	bin_number = floor( event_x_Solar/7.0); //Here divided by 7 as old x range was 350, 350/7 = 50, i.e number of histogram bins
	num_events_per_bin[bin_number] ++;
    }

    //-----Initialising histograms to be filled-------//
    //1D histograms for the energy of events
    //TH1D *h_1D_E = new TH1D("histo_E", "", 100, 0, 25);
    //TH1D *h_1D_E_cut = new TH1D("histo_E_cut", "", 100, 0, 25);

    //A 1D histogram for number over x range
    //TH1D *h_1D_x = new TH1D("histo_num", "", 50, 0, 350);
    //TH1D *h_1D_x_cut = new TH1D("histo_num_cut", "", 50, 0, 350);

    //A 2D histogram for number over x range
    //TH2D *h_2D_x_num = new TH2D("2D_histo_vis", "10000 solar events", 125, 0, 350,100, 0, 100);

    //A 2D histogram for ratio over x range -- FILTER PLOTS
    TH2D *h_2D_x_ratio = new TH2D("", "", 125, 0, 370,125, 0, 0.9);    

    //A 2D histogram for the fraction over x range -- GEOMETRY PLOTS
    //TH2D *h_2D_x_frac = new TH2D("", "", 125, 0, 370, 100, 0, 0.7);

    //--------Setting thresholds------//
    int num_over_thresh = 0; //This is a counter for the number of events which pass certain threshold conditions
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

    //Loop over events
    for(int i = 0; i < number_events; i++) {
        event_tree_Solar->GetEntry(i);
	x_pos[i] = event_x_Solar;//Array filled with x positions of events
        y_pos[i] = event_y_Solar;//Array filled with y positions of events
        Energy[i] = event_E; //Array filled with the energy of each event    

        ratio_vis_vuv = ((double) true_vis[i])/true_vuv[i]; //Total number of visible photons detected over total number of vuv photons detected
        ratio_vuv_vis = ((double) true_vuv[i])/true_vis[i];
        new_ratio = ((double) uncoated_hits[i])/coated_hits[i]; //Photons detected by uncoated opdets over photons detected by coated opdets
        ratio_coated_uncoated[i] = new_ratio;
        total_ratio += ratio_vis_vuv;

	//----Filling histograms----//
       	//Filling a 1D histogram with event E
        //h_1D_E->Fill(event_E);
	//Filling a 1D histogram with the x position of events
	//h_1D_x->Fill(event_x_Solar);
	//Fills a 2D histogram with the x-position and number of some type of photons per event
	//h_2D_x_num->Fill(event_x_Solar, uncoated_hits[i]); //Hits on uncoated opdets, i.e only visible photons
	//h_2D_x_num->Fill(event_x_Solar, coated_hits[i]); //Hits on coated opdets, i.e visible and vuv
	//h_vuv_2D->Fill(event_x_Solar, coated_hits[i]-2*uncoated_hits[i]); //Estimated number of vuv hits, factor of 2 because twice as many coated as uncoated
	//Fills a 2D histogram with x-position and fraction of regions needed to obtain some percentage of total light --REGION PLOTTING
	//h_2D_x_frac->Fill(event_x_Solar, frac_of_cell[i]);
	//Fill a 2D histogram with x-position and ratio of uncoated/coated hits -- FILTER PLOTTING
        h_2D_x_ratio->Fill(event_x_Solar, new_ratio);

	//--------CUTS--------//
	//Chain of cuts -- Can comment these out individually, or set thresholds as required, if commenting out make sure to comment out end of if statement brackets too 
        if (coated_hits[i]+uncoated_hits[i] >= photon_thresh) { //Total photon cut
	    E_cut_vec.push_back (event_E);
	    if (uncoated_hits[i] >= uncoated_photon_thresh) { //Uncoated hits cut
		uncoated_x_cut_vec.push_back (event_x_Solar);
		if (coated_hits[i] >= coated_photon_thresh) { //Coated hits cut
	            coated_x_cut_vec.push_back (event_x_Solar);
	            if (coated_hits[i]-2*uncoated_hits[i] >= difference_thresh) { //Estimated vuv cut
		        if (new_ratio < ratio_thresh) { //FILTER ratio cut
			    rat_cut_vec.push_back (new_ratio);
			    if(frac_of_cell[i] <= region_frac_thresh) { //REGION fraction cut
		                several_x_cut_vec.push_back (event_x_Solar);
		                //h2->Fill(event_x_Solar, new_ratio);
				num_over_thresh++;
			        } //End of REGION fraction cut
		            } //End of FILTER ratio cut
		        } //End of Estimated vuv cut
		    } //End of Coated hits cut
	        } //End of Uncoated hits cut
	    }//End of Total photon cut
    

	//------y-dependence test------//
        if (x_pos[i]>x_thresh && x_pos[i]<x_thresh_2) {
            //The following two statements take the modulus of y[i], as the opdet positions are symetric about 0, and then outputs the y-distance from an event to the closest pmt. 
            if (y_pos[i] < 0) {
                y_pos[i] = -1.0*y_pos[i];
            }
            y_distance = std::fmod(y_pos[i],62.0)-31;
                if (y_distance < 0) {
                    y_distance = -1.0*y_distance;
                }
            closest_y.push_back (y_distance);
	}

    }//End of event loop

    //----------PLOTTING---------//

    //This calculates the percentage of events which pass various thresholds so that can be printed on the plots
    auto size_E = E_cut_vec.size();
    auto size_uncoated = uncoated_x_cut_vec.size();
    auto size_coated = coated_x_cut_vec.size();
    auto size_several = several_x_cut_vec.size();
    double fraction_plotted_E = (double)size_E/number_events;
    double fraction_plotted_uncoated = (double)size_uncoated/number_events;
    double fraction_plotted_coated = (double)size_coated/number_events;
    double fraction_plotted_several = (double)size_several/number_events;
    double perc_rounded_E = roundf(fraction_plotted_E*1000)/10;
    double perc_rounded_uncoated = roundf(fraction_plotted_uncoated*1000)/10;
    double perc_rounded_coated = roundf(fraction_plotted_coated*1000)/10;
    double perc_rounded_several = roundf(fraction_plotted_several*1000)/10;

    mean_ratio = total_ratio/num_events;

    gStyle->SetOptStat(0);

    //Making strings
    ostringstream d_E;
    d_E << perc_rounded_E << "% remaining";
    ostringstream d_uncoated;
    d_uncoated << perc_rounded_uncoated << "% remain" << endl;
    ostringstream d_coated;
    d_coated << perc_rounded_coated << "% remaining";
    ostringstream d_several;
    d_several << perc_rounded_several << "% remain" << endl;
/*
    //Energy 1D
    h_1D_E->SetTitle("10000 solar events with and without photon cuts");
    h_1D_E->SetLabelSize(0.04,"xyz");
    h_1D_E->SetTitleSize(0.045,"xyz"); // size of axis title font
    h_1D_E->SetTitleFont(132,"xyz"); // font option
    h_1D_E->SetLabelFont(132,"xyz");
    h_1D_E->SetTitleOffset(0.95,"Y");
    h_1D_E->GetXaxis()->SetTitle("Energy (MeV)");
    h_1D_E->GetYaxis()->SetTitle("Number of events");
    h_1D_E->SetLineColor( kBlack);
    //h_1D_E_cut->SetLineColor( kBlue);
    //h_1D_E_cut->SetLineColor(kRed);
    h_1D_E->Draw();
    //h_cut->Draw("same");
    TLatex t_uncut(9,190,"#splitline{No photon cut}{}"); //If using this, will probably have to play with the first two numbers so the text appears on the plot
    //TLatex t_uncut(1.5,100,"#splitline{No photon cut}{}");
    t_uncut.SetTextSize(0.025);
    t_uncut.SetTextColor(kBlack);
    t_uncut.DrawClone();

    TLatex t_cut(10,180,"#splitline{50 total photons}{}");
    //TLatex t_cut(1.5,90,"#splitline{Several photon cuts}{}");
    t_cut.SetTextSize(0.025);
    t_cut.SetTextColor(kBlue);
    t_cut.DrawClone();

    TLatex t_portion(9,180,d_E.str().c_str());
    //TLatex t_portion(1.5,80,d_E.str().c_str());
    t_portion.SetTextSize(0.025);
    t_portion.SetTextColor(kBlue);
    t_portion.DrawClone();    

    //x-position 1D
    h_1D_x->SetLabelSize(0.04,"xyz");
    h_1D_x->SetTitleSize(0.045,"xyz"); // size of axis title font
    h_1D_x->SetTitleFont(132,"xyz"); // font option
    h_1D_x->SetLabelFont(132,"xyz");
    h_1D_x->SetTitleOffset(0.9,"Y");   

    h_1D_x->SetTitle("");
    h_1D_x->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_1D_x->GetYaxis()->SetTitle("Number of events");
    h_1D_x->SetLineColor( kBlack);
    h_1D_x->SetLineWidth(1);
    h_1D_x->Draw();
    
    h_1D_x_cut->SetLineColor( kBlue);    
    h_1D_x_cut->Draw("same");
    
    */
/*

    //2D x-position and ratio - FILTER PLOT
    h_2D_x_ratio->SetLabelSize(0.04,"xyz");
    h_2D_x_ratio->SetTitleSize(0.045,"xyz"); // size of axis title font
    h_2D_x_ratio->SetTitleFont(132,"xyz"); // font option
    h_2D_x_ratio->SetLabelFont(132,"xyz");
    h_2D_x_ratio->SetTitleOffset(0.9,"Y");
    h_2D_x_ratio->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_2D_x_ratio->GetYaxis()->SetTitle("Ratio of uncoated to coated photons");
    h_2D_x_ratio->Draw("COLZ");
    TLatex t_portion(25,0.65,d_several.str().c_str());
    t_portion.SetTextSize(0.035);
    t_portion.DrawClone();

    //2D x-position and fraction - GEOMETRY PLOT
    h_2D_x_frac->SetTitle("10000 8B events, 6 optical channels per region");
    h_2D_x_frac->SetLabelSize(0.04,"xyz");
    h_2D_x_frac->SetTitleSize(0.045,"xyz"); // size of axis title font
    h_2D_x_frac->SetTitleFont(132,"xyz"); // font option
    h_2D_x_frac->SetLabelFont(132,"xyz");
    h_2D_x_frac->SetTitleOffset(0.9,"Y");
    h_2D_x_frac->GetXaxis()->SetTitle("Distance from anode (cm)");
    h_2D_x_frac->GetYaxis()->SetTitle("Fraction of regions for 80% of total light");
    gStyle->SetOptStat(0);
    h_region->Draw("COLZ");
    TLatex t_portion(25,0.5,d_several.str().c_str());
    t_portion.SetTextSize(0.035);
    t_portion.DrawClone();
*/
    
    //c->SaveAs("../plots/geometry/combined_OR_8B_10000_region_6_0.2_0.25.png"); //Saving the plot in the canvas

    return 0;

    }
