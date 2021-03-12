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
	TString mag_field = "3.0";
	TString mat_bud = "0.05";
	TString pix_size = "10";
	gStyle->SetErrorX(0.0001);
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
                "../../output/output_vtx_res_skimmed_combined_vtx_111_B_"+mag_field+"T_vtx_"+mat_bud+"XX0_"+pix_size+"um_FastSimEvalsigma_eta_3_p_27_.root",
                "../../output/output_vtx_res_skimmed_combined_vtx_101_B_"+mag_field+"T_vtx_"+mat_bud+"XX0_"+pix_size+"um_FastSimEvalsigma_eta_3_p_27_.root",
                "../../output/output_vtx_res_skimmed_combined_vtx_110_B_"+mag_field+"T_vtx_"+mat_bud+"XX0_"+pix_size+"um_FastSimEvalsigma_eta_3_p_27_.root",
                "../../output/output_vtx_res_skimmed_combined_vtx_011_B_"+mag_field+"T_vtx_"+mat_bud+"XX0_"+pix_size+"um_FastSimEvalsigma_eta_3_p_27_.root",
                "../../output/output_vtx_res_skimmed_combined_vtx_100_B_"+mag_field+"T_vtx_"+mat_bud+"XX0_"+pix_size+"um_FastSimEvalsigma_eta_3_p_27_.root",
                "../../output/output_vtx_res_skimmed_combined_vtx_010_B_"+mag_field+"T_vtx_"+mat_bud+"XX0_"+pix_size+"um_FastSimEvalsigma_eta_3_p_27_.root",
                "../../output/output_vtx_res_skimmed_combined_vtx_001_B_"+mag_field+"T_vtx_"+mat_bud+"XX0_"+pix_size+"um_FastSimEvalsigma_eta_3_p_27_.root"
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

	TH1F *** h1_dvl_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dvt_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dvl_v_et_p_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dvt_v_et_p_bins = new TH1F ** [size_loaded];

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

		h1_dvl_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
		h1_dvt_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
		h1_dvl_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];
		h1_dvt_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];

		for(int et = 0 ; et < num_eta_bin[f] ; et++){
			h1_dvl_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dvl_v_p_et_bins_%i",et));
			h1_dvt_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dvt_v_p_et_bins_%i",et));
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			h1_dvl_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dvl_v_et_p_bins_%i",p ));
			h1_dvt_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dvt_v_et_p_bins_%i",p ));
		}
	}

	cout << "\neta bin boundaries:\n"; for(int et = 0 ; et < num_eta_bin[0]+1 ; et++) cout << (*TVT_eta_bin[0])[et] << ", "; cout << "\n";
	cout << "\np bin boundaries:\n"  ; for(int p  = 0 ; p  < num_mom_bin[0]+1 ; p ++) cout << (*TVT_mom_bin[0])[ p] << ", "; cout << "\n\n";

	// #######################################################################################################################################
        // EDIT THE CODE BELOW DEPENDING ON WHAT YOU WANT TO PLOT
	// ------------------------------------------------------------------------------
	// Copying and editing histograms
	const int eta_0 = idx_from_vector(0.,TVT_eta_bin[0]);
	int selected_bins[] = {eta_0};
	int size_selected_bins = sizeof(selected_bins)/sizeof(*selected_bins);
	
	TGraphErrors ** g_dvl_v_p_111_005 = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dvl_v_p_101_005 = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dvl_v_p_110_005 = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dvl_v_p_011_005 = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dvl_v_p_100_005 = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dvl_v_p_010_005 = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dvl_v_p_001_005 = new TGraphErrors * [size_selected_bins];

    	for(int i = 0 ; i < size_selected_bins ; i++){
		g_dvl_v_p_111_005[i] = graph_from_histo(h1_dvl_v_p_et_bins[0][selected_bins[i]],  1,20,4,100);	g_dvl_v_p_111_005[i] -> GetYaxis() -> SetMoreLogLabels();
		g_dvl_v_p_101_005[i] = graph_from_histo(h1_dvl_v_p_et_bins[1][selected_bins[i]],  2,21,4,100);	g_dvl_v_p_101_005[i] -> GetYaxis() -> SetMoreLogLabels();
		g_dvl_v_p_110_005[i] = graph_from_histo(h1_dvl_v_p_et_bins[2][selected_bins[i]], 50,21,4,100);	g_dvl_v_p_110_005[i] -> GetYaxis() -> SetMoreLogLabels();
		g_dvl_v_p_011_005[i] = graph_from_histo(h1_dvl_v_p_et_bins[3][selected_bins[i]], 96,21,4,100);	g_dvl_v_p_011_005[i] -> GetYaxis() -> SetMoreLogLabels();
		g_dvl_v_p_100_005[i] = graph_from_histo(h1_dvl_v_p_et_bins[4][selected_bins[i]],  4,22,4,100);	g_dvl_v_p_100_005[i] -> GetYaxis() -> SetMoreLogLabels();
		g_dvl_v_p_010_005[i] = graph_from_histo(h1_dvl_v_p_et_bins[5][selected_bins[i]], 62,22,4,100);	g_dvl_v_p_010_005[i] -> GetYaxis() -> SetMoreLogLabels();
		g_dvl_v_p_001_005[i] = graph_from_histo(h1_dvl_v_p_et_bins[6][selected_bins[i]], 65,22,4,100);	g_dvl_v_p_001_005[i] -> GetYaxis() -> SetMoreLogLabels();

		g_dvl_v_p_111_005[i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f, "+mat_bud+"%% X/X_{0}, "+pix_size+"#mum pixel",(*TVT_eta_bin[0])[selected_bins[i]],(*TVT_eta_bin[0])[selected_bins[i]+1]));
                g_dvl_v_p_101_005[i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f, "+mat_bud+"%% X/X_{0}, "+pix_size+"#mum pixel",(*TVT_eta_bin[0])[selected_bins[i]],(*TVT_eta_bin[0])[selected_bins[i]+1]));
                g_dvl_v_p_110_005[i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f, "+mat_bud+"%% X/X_{0}, "+pix_size+"#mum pixel",(*TVT_eta_bin[0])[selected_bins[i]],(*TVT_eta_bin[0])[selected_bins[i]+1]));
                g_dvl_v_p_011_005[i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f, "+mat_bud+"%% X/X_{0}, "+pix_size+"#mum pixel",(*TVT_eta_bin[0])[selected_bins[i]],(*TVT_eta_bin[0])[selected_bins[i]+1]));
                g_dvl_v_p_100_005[i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f, "+mat_bud+"%% X/X_{0}, "+pix_size+"#mum pixel",(*TVT_eta_bin[0])[selected_bins[i]],(*TVT_eta_bin[0])[selected_bins[i]+1]));
                g_dvl_v_p_001_005[i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f, "+mat_bud+"%% X/X_{0}, "+pix_size+"#mum pixel",(*TVT_eta_bin[0])[selected_bins[i]],(*TVT_eta_bin[0])[selected_bins[i]+1]));

		g_dvl_v_p_111_005[i] -> GetXaxis() -> SetRangeUser(0,7);
	}

	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1200,900);
	gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13); gPad -> SetLogy();
	g_dvl_v_p_111_005[0] -> Draw(   "APL");
	g_dvl_v_p_101_005[0] -> Draw("samePL");
	g_dvl_v_p_110_005[0] -> Draw("samePL");
	g_dvl_v_p_011_005[0] -> Draw("samePL");
	g_dvl_v_p_100_005[0] -> Draw("samePL");
	g_dvl_v_p_010_005[0] -> Draw("samePL");
	g_dvl_v_p_001_005[0] -> Draw("samePL");
	// ------------
	TLegend * leg1 = new TLegend(0.85,0.4,0.93,0.89);
	leg1 -> SetLineColor(0);
	leg1 -> AddEntry( g_dvl_v_p_111_005[0] , "1,1,1" );
	leg1 -> AddEntry((TObject*)0, "", "");
	leg1 -> AddEntry( g_dvl_v_p_101_005[0] , "1,0,1" );
	leg1 -> AddEntry( g_dvl_v_p_110_005[0] , "1,1,0" );
	leg1 -> AddEntry( g_dvl_v_p_011_005[0] , "0,1,1" );
	leg1 -> AddEntry( (TObject*)0, "", "");
	leg1 -> AddEntry( g_dvl_v_p_100_005[0] , "1,0,0" );
	leg1 -> AddEntry( g_dvl_v_p_010_005[0] , "0,1,0" );
	leg1 -> AddEntry( g_dvl_v_p_001_005[0] , "0,0,1" );
	// ------------
	leg1 -> Draw("same");
	c1 -> Modified();
	c1 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_simple_vtx_res_comparisons.pdf");

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
