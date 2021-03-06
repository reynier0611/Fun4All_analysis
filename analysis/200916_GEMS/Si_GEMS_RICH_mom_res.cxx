#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <filesystem>

// Root includes
#include "TRint.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVectorT.h"
#include "TGraphErrors.h"

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max );
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetErrorX(0.0001);
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		"../../output/output_skimmed_pi-_det2_20x20_Beast_FastTrackingEval.root",			//  0
		"../../output/output_skimmed_pi-_det2_10x10_Beast_FastTrackingEval.root",			//  1
		"../../output/output_skimmed_pi-_det2_both_GEMs_20x20_Beast_FastTrackingEval.root",		//  2
		"../../output/output_skimmed_pi-_det2_both_GEMs_10x10_Beast_FastTrackingEval.root",		//  3
		"../../output/output_skimmed_pi-_det2_20x20_sPHENIX_FastTrackingEval.root",			//  4
		"../../output/output_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEval.root",			//  5
		"../../output/output_skimmed_pi-_det2_both_GEMs_20x20_sPHENIX_FastTrackingEval.root",		//  6
		"../../output/output_skimmed_pi-_det2_both_GEMs_10x10_sPHENIX_FastTrackingEval.root",		//  7
		"../../output/output_skimmed_pi-_det2_both_GEMs_RICH_20x20_Beast_FastTrackingEval.root",	//  8
		"../../output/output_skimmed_pi-_det2_both_GEMs_RICH_10x10_Beast_FastTrackingEval.root",	//  9
		"../../output/output_skimmed_pi-_det2_both_GEMs_RICH_20x20_sPHENIX_FastTrackingEval.root",	// 10
                "../../output/output_skimmed_pi-_det2_both_GEMs_RICH_10x10_sPHENIX_FastTrackingEval.root",	// 11
		"../../output/output_skimmed_combined_pi-_det2_10umGEM_RICH_20x20_Beast_FastTrackingEval.root",	// 12
		"../../output/output_skimmed_combined_pi-_det2_10umGEM_RICH_10x10_Beast_FastTrackingEval.root",	// 13
		"../../output/output_skimmed_combined_pi-_det2_LoResGEM_RICH_20x20_Beast_FastTrackingEval.root",// 14
		"../../output/output_skimmed_combined_pi-_det2_LoResGEM_RICH_10x10_Beast_FastTrackingEval.root"	// 15
	};
	// #######################################################################################################################################
	// YOU SHOULDN'T NEED TO MODIFY ANYTHING IN THE BLOCK OF CODE BELOW AND UNTIL AFTER THE NEXT LINE WITH ###...
	const int size_loaded = sizeof(fnames)/sizeof(*fnames);
	// ------------------------------------------------------------------------------
	// Preparing variables that will later on be filled from root files
	TVectorT<double> ** TVT_eta_bin = new TVectorT<double>*[size_loaded];
	TVectorT<double> ** TVT_mom_bin = new TVectorT<double>*[size_loaded];
	int num_eta_bin[size_loaded] = {0};
	int num_mom_bin[size_loaded] = {0};
	// ------------------------------------------------------------------------------
	// Loading root files and info therein
	TFile ** Fin = new TFile*[size_loaded];

	TH1F *** h1_dpp_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dth_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dph_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dpp_v_et_p_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dth_v_et_p_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dph_v_et_p_bins = new TH1F ** [size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		ifstream fin;
		fin.open(fnames[f]);
		if(!fin){ cout << "\033[1;31mCouldn't find input file '" << fnames[f] << "'. Bailing out!\033[0m" << endl; exit(0);}
		fin.close();

		Fin[f] = new TFile(fnames[f]);

		TVT_eta_bin[f] = (TVectorT<double> *) Fin[f] -> Get("TVT_eta_bin");
		TVT_mom_bin[f] = (TVectorT<double> *) Fin[f] -> Get("TVT_mom_bin");

		num_eta_bin[f] = (*TVT_eta_bin[f]).GetNoElements()-1;
		num_mom_bin[f] = (*TVT_mom_bin[f]).GetNoElements()-1;

		h1_dpp_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
		h1_dth_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
		h1_dph_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
		h1_dpp_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];
		h1_dth_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];
		h1_dph_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];

		for(int et = 0 ; et < num_eta_bin[f] ; et++){
			h1_dpp_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dpp_v_p_et_bins_%i",et));
			h1_dth_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dth_v_p_et_bins_%i",et));
			h1_dph_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dph_v_p_et_bins_%i",et));
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			h1_dpp_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dpp_v_et_p_bins_%i",p ));
			h1_dth_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dth_v_et_p_bins_%i",p ));
			h1_dph_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dph_v_et_p_bins_%i",p ));
		}
	}

	cout << "\neta bin boundaries:\n"; for(int et = 0 ; et < num_eta_bin[0]+1 ; et++) cout << (*TVT_eta_bin[0])[et] << ", "; cout << "\n";
	cout << "\np bin boundaries:\n"  ; for(int p  = 0 ; p  < num_mom_bin[0]+1 ; p ++) cout << (*TVT_mom_bin[0])[ p] << ", "; cout << "\n\n";

	// #######################################################################################################################################
        // EDIT THE CODE BELOW DEPENDING ON WHAT YOU WANT TO PLOT
	// ------------------------------------------------------------------------------
	// Copying and editing histograms
	double mom_bin[] = {4.,10.,25.};
	const int p_4GeV  = idx_from_vector( 4.,TVT_mom_bin[0]);
	const int p_10GeV = idx_from_vector(10.,TVT_mom_bin[0]);
	const int p_25GeV = idx_from_vector(25.,TVT_mom_bin[0]);
	int selected_bins[] = {p_4GeV,p_10GeV,p_25GeV};
	int size_selected_bins = sizeof(selected_bins)/sizeof(*selected_bins);

	TGraphErrors ** g_dpp_v_et_selected_20um_Beast_si              = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_10um_Beast_si              = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_20um_Beast_si_GEM          = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_10um_Beast_si_GEM          = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_20um_Beast_si_GEM_RICH     = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_10um_Beast_si_GEM_RICH     = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_20um_Beast_si_10umGEM_RICH = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_10um_Beast_si_10umGEM_RICH = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_20um_Beast_si_loReGEM_RICH = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_10um_Beast_si_loReGEM_RICH = new TGraphErrors * [size_selected_bins];

	TGraphErrors ** g_dpp_v_et_selected_20um_sPHENIX_si            = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_10um_sPHENIX_si            = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_20um_sPHENIX_si_GEM        = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_10um_sPHENIX_si_GEM        = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_20um_sPHENIX_si_GEM_RICH   = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_10um_sPHENIX_si_GEM_RICH   = new TGraphErrors * [size_selected_bins];
        
    	for(int i = 0 ; i < size_selected_bins ; i++){
		double max_val = 0.46*mom_bin[i]+5.33;

		g_dpp_v_et_selected_20um_Beast_si             [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 0][selected_bins[i]] ,66,21,0,max_val);
		g_dpp_v_et_selected_10um_Beast_si             [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 1][selected_bins[i]] ,94,21,0,max_val);
		g_dpp_v_et_selected_20um_Beast_si_GEM         [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 2][selected_bins[i]] ,62,20,0,max_val);
		g_dpp_v_et_selected_10um_Beast_si_GEM         [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 3][selected_bins[i]] , 2,20,0,max_val);
		g_dpp_v_et_selected_20um_Beast_si_GEM_RICH    [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 8][selected_bins[i]] , 4,22,0,max_val);
		g_dpp_v_et_selected_10um_Beast_si_GEM_RICH    [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 9][selected_bins[i]] ,50,22,0,max_val);
		g_dpp_v_et_selected_20um_Beast_si_10umGEM_RICH[i] = graph_from_histo( h1_dpp_v_et_p_bins[12][selected_bins[i]] , 1,20,0,max_val);
		g_dpp_v_et_selected_10um_Beast_si_10umGEM_RICH[i] = graph_from_histo( h1_dpp_v_et_p_bins[13][selected_bins[i]] , 1,20,0,max_val);
		g_dpp_v_et_selected_20um_Beast_si_loReGEM_RICH[i] = graph_from_histo( h1_dpp_v_et_p_bins[14][selected_bins[i]] , 8,23,0,max_val);
		g_dpp_v_et_selected_10um_Beast_si_loReGEM_RICH[i] = graph_from_histo( h1_dpp_v_et_p_bins[15][selected_bins[i]] , 8,23,0,max_val);

		g_dpp_v_et_selected_20um_sPHENIX_si           [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 4][selected_bins[i]] ,66,21,0,max_val);
		g_dpp_v_et_selected_10um_sPHENIX_si           [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 5][selected_bins[i]] ,94,21,0,max_val);
		g_dpp_v_et_selected_20um_sPHENIX_si_GEM       [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 6][selected_bins[i]] ,62,20,0,max_val);
		g_dpp_v_et_selected_10um_sPHENIX_si_GEM       [i] = graph_from_histo( h1_dpp_v_et_p_bins[ 7][selected_bins[i]] , 2,20,0,max_val);
		g_dpp_v_et_selected_20um_sPHENIX_si_GEM_RICH  [i] = graph_from_histo( h1_dpp_v_et_p_bins[10][selected_bins[i]] , 4,22,0,max_val);
		g_dpp_v_et_selected_10um_sPHENIX_si_GEM_RICH  [i] = graph_from_histo( h1_dpp_v_et_p_bins[11][selected_bins[i]] ,50,22,0,max_val);

		g_dpp_v_et_selected_20um_Beast_si      [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                g_dpp_v_et_selected_10um_Beast_si      [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                g_dpp_v_et_selected_20um_Beast_si_GEM  [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
	        g_dpp_v_et_selected_10um_Beast_si_GEM  [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));

                g_dpp_v_et_selected_20um_sPHENIX_si    [i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                g_dpp_v_et_selected_10um_sPHENIX_si    [i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                g_dpp_v_et_selected_20um_sPHENIX_si_GEM[i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
		g_dpp_v_et_selected_10um_sPHENIX_si_GEM[i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
	}

	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1200,900);
	gPad -> SetRightMargin(0.02); gPad -> SetGridx(); gPad -> SetGridy(); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
	g_dpp_v_et_selected_20um_Beast_si         [2] -> Draw(   "APL");
	g_dpp_v_et_selected_10um_Beast_si         [2] -> Draw("samePL");
	g_dpp_v_et_selected_20um_Beast_si_GEM     [2] -> Draw("samePL");
	g_dpp_v_et_selected_10um_Beast_si_GEM     [2] -> Draw("samePL");
	g_dpp_v_et_selected_20um_Beast_si_GEM_RICH[2] -> Draw("samePL");
	g_dpp_v_et_selected_10um_Beast_si_GEM_RICH[2] -> Draw("samePL");
	// ------------
	TLegend * leg1 = new TLegend(0.35,0.5,0.75,0.89);
	leg1 -> SetLineColor(0);
	leg1 -> AddEntry( g_dpp_v_et_selected_20um_Beast_si         [2] , "All-Si (20 #mum)"             );
	leg1 -> AddEntry( g_dpp_v_et_selected_10um_Beast_si         [2] , "All-Si (10 #mum)"             );
	leg1 -> AddEntry( g_dpp_v_et_selected_20um_Beast_si_GEM     [2] , "All-Si (20 #mum) + GEM"       );
	leg1 -> AddEntry( g_dpp_v_et_selected_10um_Beast_si_GEM     [2] , "All-Si (10 #mum) + GEM"       );
	leg1 -> AddEntry( g_dpp_v_et_selected_20um_Beast_si_GEM_RICH[2] , "All-Si (20 #mum) + GEM + RICH");
	leg1 -> AddEntry( g_dpp_v_et_selected_10um_Beast_si_GEM_RICH[2] , "All-Si (10 #mum) + GEM + RICH");
	// ------------
	leg1 -> Draw("same");
	c1 -> Modified();
	c1 -> Update();
	// ----------------------------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1200,900);
	gPad -> SetRightMargin(0.02); gPad -> SetGridx(); gPad -> SetGridy(); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
	g_dpp_v_et_selected_20um_sPHENIX_si         [2] -> Draw(   "APL");
	g_dpp_v_et_selected_10um_sPHENIX_si         [2] -> Draw("samePL");
	g_dpp_v_et_selected_20um_sPHENIX_si_GEM     [2] -> Draw("samePL");
	g_dpp_v_et_selected_10um_sPHENIX_si_GEM     [2] -> Draw("samePL");
	g_dpp_v_et_selected_20um_sPHENIX_si_GEM_RICH[2] -> Draw("samePL");
	g_dpp_v_et_selected_10um_sPHENIX_si_GEM_RICH[2] -> Draw("samePL");
	leg1 -> Draw("same");
	c2 -> Modified();
	c2 -> Update();
	// ----------------------------------------------
        TCanvas * c3 = new TCanvas("c3","c3",1200,900);
	c3 -> Divide(2,2);
	for(int i = 0 ; i < 3 ; i++){
		c3 -> cd(i+1); gPad -> SetGridx(); gPad -> SetGridy(); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
        	g_dpp_v_et_selected_20um_Beast_si         [i] -> Draw(   "APL");
                g_dpp_v_et_selected_10um_Beast_si         [i] -> Draw("samePL");
                g_dpp_v_et_selected_20um_Beast_si_GEM     [i] -> Draw("samePL");
                g_dpp_v_et_selected_10um_Beast_si_GEM     [i] -> Draw("samePL");
                g_dpp_v_et_selected_20um_Beast_si_GEM_RICH[i] -> Draw("samePL");
                g_dpp_v_et_selected_10um_Beast_si_GEM_RICH[i] -> Draw("samePL");
		g_dpp_v_et_selected_20um_Beast_si_10umGEM_RICH[i] -> Draw("samePL");
	}
	c3 -> cd(4);
	leg1 -> Draw("same");
	c3 -> Modified();
	c3 -> Update();
	// ----------------------------------------------
        TCanvas * c4 = new TCanvas("c4","c4",1200,900);
        c4 -> Divide(2,2);
        for(int i = 0 ; i < 3 ; i++){
                c4 -> cd(i+1); gPad -> SetGridx(); gPad -> SetGridy(); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
                g_dpp_v_et_selected_20um_sPHENIX_si         [i] -> Draw(   "APL");
                g_dpp_v_et_selected_10um_sPHENIX_si         [i] -> Draw("samePL");
                g_dpp_v_et_selected_20um_sPHENIX_si_GEM     [i] -> Draw("samePL");
                g_dpp_v_et_selected_10um_sPHENIX_si_GEM     [i] -> Draw("samePL");
                g_dpp_v_et_selected_20um_sPHENIX_si_GEM_RICH[i] -> Draw("samePL");
                g_dpp_v_et_selected_10um_sPHENIX_si_GEM_RICH[i] -> Draw("samePL");
        }
        c4 -> cd(4);
        leg1 -> Draw("same");
        c4 -> Modified();
        c4 -> Update();
	// ----------------------------------------------
	double xDum1[] = { 1.5, 3.6};
	double xDum2[] = {-3.6,-1.5};
	double yDum [] = { 0.0, 17.};
	TGraph * gDum1 = new TGraph(2.,xDum1,yDum);
	TGraph * gDum2 = new TGraph(2.,xDum2,yDum);

	int i = 2;
	gDum1 -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
	gDum1 -> GetXaxis() -> SetTitle("#eta");
	gDum1 -> GetXaxis() -> CenterTitle();
	gDum1 -> GetXaxis() -> SetNdivisions(108);
	gDum1 -> GetXaxis() -> SetLabelSize(0.06);
	gDum1 -> GetXaxis() -> SetTitleSize(0.06);
	gDum1 -> GetYaxis() -> SetTitle("dp / p [%]");
	gDum1 -> GetYaxis() -> CenterTitle();
	gDum1 -> GetYaxis() -> SetNdivisions(108);
	gDum1 -> GetYaxis() -> SetLabelSize(0.06);
        gDum1 -> GetYaxis() -> SetTitleSize(0.06);
	gDum2 -> SetTitle("");
        gDum2 -> GetXaxis() -> SetTitle("#eta");
        gDum2 -> GetXaxis() -> CenterTitle();
        gDum2 -> GetXaxis() -> SetNdivisions(108);
        gDum2 -> GetXaxis() -> SetLabelSize(0.06);
        gDum2 -> GetXaxis() -> SetTitleSize(0.06);
	gDum2 -> GetYaxis() -> SetTitle("dp / p [%]");
        gDum2 -> GetYaxis() -> CenterTitle();
        gDum2 -> GetYaxis() -> SetNdivisions(108);
	gDum2 -> GetYaxis() -> SetLabelSize(0.06);
        gDum2 -> GetYaxis() -> SetTitleSize(0.06);

        TCanvas * c5 = new TCanvas("c5","c5",1000,900);
        c5 -> Divide(1,2);
	c5 -> cd(1);
	gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetTopMargin(0.01);
	gDum1 -> Draw("AP");
        g_dpp_v_et_selected_20um_Beast_si         [2] -> Draw("samePL");
        //g_dpp_v_et_selected_20um_Beast_si_GEM     [2] -> Draw("samePL");
        g_dpp_v_et_selected_20um_Beast_si_GEM_RICH[2] -> Draw("samePL");
        g_dpp_v_et_selected_20um_Beast_si_10umGEM_RICH[2] -> Draw("samePL");
	//g_dpp_v_et_selected_20um_Beast_si_loReGEM_RICH[2] -> Draw("samePL");
	// ------------
        TLegend * leg2 = new TLegend(0.20,0.45,0.65,0.89);
        leg2 -> SetLineColor(0);
        leg2 -> AddEntry( g_dpp_v_et_selected_20um_Beast_si            [2] , "All-Si (20 #mum)"             );
	//leg2 -> AddEntry(g_dpp_v_et_selected_20um_Beast_si_loReGEM_RICH[2] , "All-Si (20 #mum) + GEM (#sigma_{r} = 1cm/#sqrt{12}, #sigma_{#phi} = 70 #mum)");
        leg2 -> AddEntry( g_dpp_v_et_selected_20um_Beast_si_GEM_RICH   [2] , "All-Si (20 #mum) + GEM (#sigma = 50 #mum)");
	leg2 -> AddEntry(g_dpp_v_et_selected_20um_Beast_si_10umGEM_RICH[2] , "All-Si (20 #mum) + Si disk (10 #mum)");
        // ------------
        leg2 -> Draw("same");
	// ------------
	c5 -> cd(2);
        gPad -> SetRightMargin(0.02); /*gPad -> SetGridx(); gPad -> SetGridy();*/ gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetTopMargin(0.02);
        gDum2 -> Draw("AP");
	g_dpp_v_et_selected_20um_Beast_si         [2] -> Draw("samePL");
        g_dpp_v_et_selected_20um_Beast_si_GEM_RICH[2] -> Draw("samePL");
        g_dpp_v_et_selected_20um_Beast_si_10umGEM_RICH[2] -> Draw("samePL"); 
	//g_dpp_v_et_selected_20um_Beast_si_loReGEM_RICH[2] -> Draw("samePL");
	c5 -> Modified();
        c5 -> Update();
	// ------------------
	TCanvas * c6 = new TCanvas("c6","c6",1000,900);
        c6 -> Divide(1,2);
        c6 -> cd(1);
        gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetTopMargin(0.01);
        gDum1 -> Draw("AP");
        g_dpp_v_et_selected_10um_Beast_si         [2] -> Draw("samePL");
        g_dpp_v_et_selected_10um_Beast_si_GEM_RICH[2] -> Draw("samePL");
        g_dpp_v_et_selected_10um_Beast_si_10umGEM_RICH[2] -> Draw("samePL");
        //g_dpp_v_et_selected_10um_Beast_si_loReGEM_RICH[2] -> Draw("samePL");
        // ------------
	TLegend * leg3 = new TLegend(0.20,0.45,0.65,0.89);
        leg3 -> SetLineColor(0);
        leg3 -> AddEntry( g_dpp_v_et_selected_10um_Beast_si            [2] , "All-Si (10 #mum)"             );
        leg3 -> AddEntry( g_dpp_v_et_selected_10um_Beast_si_GEM_RICH   [2] , "All-Si (10 #mum) + GEM (#sigma = 50 #mum)");
        leg3 -> AddEntry(g_dpp_v_et_selected_10um_Beast_si_10umGEM_RICH[2] , "All-Si (10 #mum) + Si disk (10 #mum)");
	// ------------
        leg3 -> Draw("same");
        // ------------
        c6 -> cd(2);
        gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetTopMargin(0.02);
        gDum2 -> Draw("AP");
        g_dpp_v_et_selected_10um_Beast_si         [2] -> Draw("samePL"); 
        g_dpp_v_et_selected_10um_Beast_si_GEM_RICH[2] -> Draw("samePL");
        g_dpp_v_et_selected_10um_Beast_si_10umGEM_RICH[2] -> Draw("samePL");
        //g_dpp_v_et_selected_10um_Beast_si_loReGEM_RICH[2] -> Draw("samePL");
        c6 -> Modified();
        c6 -> Update();
	// ------------------------------------------------------------------------------
	double etas_to_check[] = {-3.5,-3.0,-2.5,-2.0,-1.5,1.5,2.0,2.5,3.0,3.5};
	for(int i = 0 ; i < sizeof(etas_to_check)/sizeof(*etas_to_check) ; i++){
		cout << etas_to_check[i] << "\t";
		cout << 100*(1-g_dpp_v_et_selected_20um_Beast_si_GEM_RICH    [2]->Eval(etas_to_check[i])/g_dpp_v_et_selected_20um_Beast_si[2]->Eval(etas_to_check[i])) << "\t";
		cout << 100*(1-g_dpp_v_et_selected_20um_Beast_si_10umGEM_RICH[2]->Eval(etas_to_check[i])/g_dpp_v_et_selected_20um_Beast_si[2]->Eval(etas_to_check[i])) << "\t";
		cout << 100*(1-g_dpp_v_et_selected_20um_Beast_si_loReGEM_RICH[2]->Eval(etas_to_check[i])/g_dpp_v_et_selected_20um_Beast_si[2]->Eval(etas_to_check[i])) << "\n";
	}
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_Si_GEMS_RICH_mom_res.pdf(");
	c2 -> Print("results_Si_GEMS_RICH_mom_res.pdf" );
	c3 -> Print("results_Si_GEMS_RICH_mom_res.pdf" );
	c4 -> Print("results_Si_GEMS_RICH_mom_res.pdf" );
	c5 -> Print("results_Si_GEMS_RICH_mom_res.pdf)");

	myapp -> Run();
	return 0;
}
// ============================================================================================================================================
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max ){
	h1 -> SetLineWidth(2);
	h1 -> SetLineColor(color);
	h1 -> SetMarkerStyle(marker);
	h1 -> SetMarkerColor(color);

	if(min!=999) h1 -> SetMinimum(min);
	if(max!=999) h1 -> SetMaximum(max);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetNdivisions(108); // to draw less tick marks
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetNdivisions(108); // to draw less tick marks

	h1 -> SetMinimum(0.001);
}
// ============================================================================================================================================
int idx_from_vector( double value , TVectorT<double> * vec ){
	int size_vec = (*vec).GetNoElements();
	double diff = 999.;
	int idx = -1;
	for(int i = 0 ; i < size_vec ; i++){
		if(abs((*vec)[i] - value)<diff){
			diff = abs((*vec)[i] - value);
			idx = i;
		}
	}
	return idx;
}
// ============================================================================================================================================
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max ){
	// Gathering information from histogram
	float yval[500] = {0};
	float xval[500] = {0};
	float err [500] = {0};
	int ctr = 0;

	TString xax = h1 -> GetXaxis() -> GetTitle();
	TString yax = h1 -> GetYaxis() -> GetTitle();
	TString tit = h1 -> GetTitle();

	const int nbin = h1 -> GetSize()-2;

	for(int i = 0 ; i < nbin ; i++){
		if(h1 -> GetBinContent(i+1)>0){
			yval[ctr] = h1 -> GetBinContent(i+1);
			xval[ctr] = h1 -> GetXaxis() -> GetBinCenter(i+1);
			err[ctr] = h1 -> GetBinError(i+1);
			ctr++;
		}
	}

	// Creating a TGraphErrors object that mirrors the histogram
	TGraphErrors * g_l1 = new TGraphErrors(ctr,xval,yval,0,err);
	g_l1->SetLineColor(color);
	g_l1->SetMarkerColor(color);
	g_l1->SetMarkerStyle(marker);
	g_l1->SetMarkerSize(1.4);
	g_l1->SetLineWidth(2);

	g_l1->GetXaxis()->SetTitle(xax);
	g_l1->GetXaxis()->SetNdivisions(108);
	g_l1->GetXaxis()->SetTitleSize(0.05);
	g_l1->GetXaxis()->SetLabelSize(0.05);
	g_l1->GetXaxis()->CenterTitle();

	g_l1->GetYaxis()->SetTitle(yax);
	g_l1->GetYaxis()->SetNdivisions(108);
	g_l1->GetYaxis()->SetTitleSize(0.05);
	g_l1->GetYaxis()->SetLabelSize(0.05);
	g_l1->GetYaxis()->CenterTitle();

	g_l1->SetTitle(tit);

	if(min!=999) g_l1 -> SetMinimum(min);
        if(max!=999) g_l1 -> SetMaximum(max);

	float minxval = h1 -> GetBinLowEdge(1);
	float maxxval = h1 -> GetBinLowEdge(nbin) + h1 -> GetBinWidth(nbin); 

	g_l1->GetXaxis()->SetRangeUser( minxval , maxxval );

	return g_l1;
}
