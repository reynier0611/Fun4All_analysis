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
#include "TGraph.h"

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
TGraph * graph_from_histo( TH1F * h1 );
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
		"output/output_skimmed_pi-_det2_20x20_Beast_FastTrackingEval.root",
		"output/output_skimmed_pi-_det2_10x10_Beast_FastTrackingEval.root",
		"output/output_skimmed_pi-_det2_both_GEMs_20x20_Beast_FastTrackingEval.root",
		"output/output_skimmed_pi-_det2_both_GEMs_10x10_Beast_FastTrackingEval.root",
		"output/output_skimmed_pi-_det2_20x20_sPHENIX_FastTrackingEval.root",
		"output/output_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEval.root",
		"output/output_skimmed_pi-_det2_both_GEMs_20x20_sPHENIX_FastTrackingEval.root",
		"output/output_skimmed_pi-_det2_both_GEMs_10x10_sPHENIX_FastTrackingEval.root"
	};
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

	// ------------------------------------------------------------------------------
	// Copying and editing histograms
	double mom_bin[] = {4.,10.,25.};
	const int p_4GeV  = idx_from_vector( 4.,TVT_mom_bin[0]);
	const int p_10GeV = idx_from_vector(10.,TVT_mom_bin[0]);
	const int p_25GeV = idx_from_vector(25.,TVT_mom_bin[0]);
	int selected_bins[] = {p_4GeV,p_10GeV,p_25GeV};
	int size_selected_bins = sizeof(selected_bins)/sizeof(*selected_bins);

	TH1F ** h1_dpp_v_et_selected_20um_Beast_si       = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_20um_Beast_si       = new TGraph * [size_selected_bins];
        TH1F ** h1_dpp_v_et_selected_10um_Beast_si       = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_10um_Beast_si       = new TGraph * [size_selected_bins];
        TH1F ** h1_dpp_v_et_selected_20um_Beast_si_GEM   = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_20um_Beast_si_GEM   = new TGraph * [size_selected_bins];
        TH1F ** h1_dpp_v_et_selected_10um_Beast_si_GEM   = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_10um_Beast_si_GEM   = new TGraph * [size_selected_bins];

	TH1F ** h1_dpp_v_et_selected_20um_sPHENIX_si     = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_20um_sPHENIX_si     = new TGraph * [size_selected_bins];
        TH1F ** h1_dpp_v_et_selected_10um_sPHENIX_si     = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_10um_sPHENIX_si     = new TGraph * [size_selected_bins];
        TH1F ** h1_dpp_v_et_selected_20um_sPHENIX_si_GEM = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_20um_sPHENIX_si_GEM = new TGraph * [size_selected_bins];
        TH1F ** h1_dpp_v_et_selected_10um_sPHENIX_si_GEM = new TH1F *[size_selected_bins];	TGraph ** g_dpp_v_et_selected_10um_sPHENIX_si_GEM = new TGraph * [size_selected_bins];
        
    	for(int i = 0 ; i < size_selected_bins ; i++){
		double max_val = 0.46*mom_bin[i]+5.33;

		h1_dpp_v_et_selected_20um_Beast_si      [i] = (TH1F*) h1_dpp_v_et_p_bins[0][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_20um_Beast_si      [i],62,21,999,max_val);
		h1_dpp_v_et_selected_10um_Beast_si      [i] = (TH1F*) h1_dpp_v_et_p_bins[1][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_10um_Beast_si      [i],94,21,999,max_val);
		h1_dpp_v_et_selected_20um_Beast_si_GEM  [i] = (TH1F*) h1_dpp_v_et_p_bins[2][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_20um_Beast_si_GEM  [i], 4,20,999,max_val);
		h1_dpp_v_et_selected_10um_Beast_si_GEM  [i] = (TH1F*) h1_dpp_v_et_p_bins[3][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_10um_Beast_si_GEM  [i], 2,20,999,max_val);

		h1_dpp_v_et_selected_20um_sPHENIX_si    [i] = (TH1F*) h1_dpp_v_et_p_bins[4][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_20um_sPHENIX_si    [i],62,21,999,max_val);
                h1_dpp_v_et_selected_10um_sPHENIX_si    [i] = (TH1F*) h1_dpp_v_et_p_bins[5][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_10um_sPHENIX_si    [i],94,21,999,max_val);
                h1_dpp_v_et_selected_20um_sPHENIX_si_GEM[i] = (TH1F*) h1_dpp_v_et_p_bins[6][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_20um_sPHENIX_si_GEM[i], 4,20,999,max_val);
                h1_dpp_v_et_selected_10um_sPHENIX_si_GEM[i] = (TH1F*) h1_dpp_v_et_p_bins[7][selected_bins[i]] -> Clone();	prettyTH1F(h1_dpp_v_et_selected_10um_sPHENIX_si_GEM[i], 2,20,999,max_val);

		g_dpp_v_et_selected_20um_Beast_si      [i] = graph_from_histo( h1_dpp_v_et_selected_20um_Beast_si      [i] );
		g_dpp_v_et_selected_10um_Beast_si      [i] = graph_from_histo( h1_dpp_v_et_selected_10um_Beast_si      [i] );
		g_dpp_v_et_selected_20um_Beast_si_GEM  [i] = graph_from_histo( h1_dpp_v_et_selected_20um_Beast_si_GEM  [i] );
		g_dpp_v_et_selected_10um_Beast_si_GEM  [i] = graph_from_histo( h1_dpp_v_et_selected_10um_Beast_si_GEM  [i] );
	                                                                                                                        
		g_dpp_v_et_selected_20um_sPHENIX_si    [i] = graph_from_histo( h1_dpp_v_et_selected_20um_sPHENIX_si    [i] );
		g_dpp_v_et_selected_10um_sPHENIX_si    [i] = graph_from_histo( h1_dpp_v_et_selected_10um_sPHENIX_si    [i] );
		g_dpp_v_et_selected_20um_sPHENIX_si_GEM[i] = graph_from_histo( h1_dpp_v_et_selected_20um_sPHENIX_si_GEM[i] );
		g_dpp_v_et_selected_10um_sPHENIX_si_GEM[i] = graph_from_histo( h1_dpp_v_et_selected_10um_sPHENIX_si_GEM[i] );

		h1_dpp_v_et_selected_20um_Beast_si      [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                h1_dpp_v_et_selected_10um_Beast_si      [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                h1_dpp_v_et_selected_20um_Beast_si_GEM  [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
	        h1_dpp_v_et_selected_10um_Beast_si_GEM  [i] -> SetTitle(Form("Beast (3.0 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                                                                                                                                                                               
                h1_dpp_v_et_selected_20um_sPHENIX_si    [i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                h1_dpp_v_et_selected_10um_sPHENIX_si    [i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
                h1_dpp_v_et_selected_20um_sPHENIX_si_GEM[i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
		h1_dpp_v_et_selected_10um_sPHENIX_si_GEM[i] -> SetTitle(Form("BaBar (1.4 T), %.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[selected_bins[i]],(*TVT_mom_bin[0])[selected_bins[i]+1]));
	}

	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1200,900);
	gPad -> SetRightMargin(0.02); gPad -> SetGridx(); gPad -> SetGridy();
	h1_dpp_v_et_selected_20um_Beast_si    [2] -> Draw(      );	g_dpp_v_et_selected_20um_Beast_si    [2] -> Draw("sameL");
	h1_dpp_v_et_selected_10um_Beast_si    [2] -> Draw("same");	g_dpp_v_et_selected_10um_Beast_si    [2] -> Draw("sameL");
	h1_dpp_v_et_selected_20um_Beast_si_GEM[2] -> Draw("same");	g_dpp_v_et_selected_20um_Beast_si_GEM[2] -> Draw("sameL");
	h1_dpp_v_et_selected_10um_Beast_si_GEM[2] -> Draw("same");	g_dpp_v_et_selected_10um_Beast_si_GEM[2] -> Draw("sameL");
	// ------------
	TLegend * leg1 = new TLegend(0.35,0.5,0.75,0.89);
	leg1 -> SetLineColor(0);
	leg1 -> AddEntry( h1_dpp_v_et_selected_20um_Beast_si    [2] , "All-Si (20 #mum)"       );
	leg1 -> AddEntry( h1_dpp_v_et_selected_10um_Beast_si    [2] , "All-Si (10 #mum)"       );
	leg1 -> AddEntry( h1_dpp_v_et_selected_20um_Beast_si_GEM[2] , "All-Si (20 #mum) + GEM" );
	leg1 -> AddEntry( h1_dpp_v_et_selected_10um_Beast_si_GEM[2] , "All-Si (10 #mum) + GEM" );
	// ------------
	leg1 -> Draw("same");
	c1 -> Modified();
	c1 -> Update();
	// ----------------------------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1200,900);
	gPad -> SetRightMargin(0.02); gPad -> SetGridx(); gPad -> SetGridy();
	h1_dpp_v_et_selected_20um_sPHENIX_si    [2] -> Draw(      );	g_dpp_v_et_selected_20um_sPHENIX_si    [2] -> Draw("sameL");
	h1_dpp_v_et_selected_10um_sPHENIX_si    [2] -> Draw("same");	g_dpp_v_et_selected_10um_sPHENIX_si    [2] -> Draw("sameL");
	h1_dpp_v_et_selected_20um_sPHENIX_si_GEM[2] -> Draw("same");	g_dpp_v_et_selected_20um_sPHENIX_si_GEM[2] -> Draw("sameL");
	h1_dpp_v_et_selected_10um_sPHENIX_si_GEM[2] -> Draw("same");	g_dpp_v_et_selected_10um_sPHENIX_si_GEM[2] -> Draw("sameL");
	leg1 -> Draw("same");
	c2 -> Modified();
	c2 -> Update();
	// ----------------------------------------------
        TCanvas * c3 = new TCanvas("c3","c3",1200,900);
	c3 -> Divide(2,2);
	for(int i = 0 ; i < 3 ; i++){
		c3 -> cd(i+1); gPad -> SetGridx(); gPad -> SetGridy();
        	h1_dpp_v_et_selected_20um_Beast_si    [i] -> Draw(      );      g_dpp_v_et_selected_20um_Beast_si    [i] -> Draw("sameL");
        	h1_dpp_v_et_selected_10um_Beast_si    [i] -> Draw("same");      g_dpp_v_et_selected_10um_Beast_si    [i] -> Draw("sameL");
        	h1_dpp_v_et_selected_20um_Beast_si_GEM[i] -> Draw("same");      g_dpp_v_et_selected_20um_Beast_si_GEM[i] -> Draw("sameL");
        	h1_dpp_v_et_selected_10um_Beast_si_GEM[i] -> Draw("same");      g_dpp_v_et_selected_10um_Beast_si_GEM[i] -> Draw("sameL");
	}
	c3 -> cd(4);
	leg1 -> Draw("same");
	c3 -> Modified();
	c3 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_plotter1.pdf(");
	c2 -> Print("results_plotter1.pdf" );
	c3 -> Print("results_plotter1.pdf)");

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
TGraph * graph_from_histo( TH1F * h1 ){
	float res[100] = {0};
	float eta[100] = {0};
	int ctr = 0;

	int color = h1 -> GetMarkerColor();

	for(int i = 0 ; i < h1 -> GetSize()-2 ; i++){
		res[ctr] = h1 -> GetBinContent(i+1);
		eta[ctr] = h1 -> GetXaxis() -> GetBinCenter(i+1);
		ctr++;
	}

	TGraph * g_l1 = new TGraph(ctr,eta,res);
	g_l1->SetLineColor(color);
	g_l1->SetLineWidth(2);

	return g_l1;
}
