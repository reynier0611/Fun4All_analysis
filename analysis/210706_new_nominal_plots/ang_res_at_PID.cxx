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
#include "TLatex.h"
#include "TLine.h"

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
		"output_mom_res_projection_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEvalsigma_eta_16_p_10_.root",
		"output_mom_res_projection_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_RICH_GEM_res_50.0um_B_ATHENA_210507_FastSimEvalsigma_eta_16_p_10_.root",
		"output_mom_res_projection_skimmed_combined_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_BaBar_FastSimEvalsigma_eta_16_p_10_.root"
	};

	TString label[] = {
                "B = 3.0 T (ATHENA 05/07/21) (all-si)",
                "B = 3.0 T (ATHENA 05/07/21) (all-si+GEMs)",
                "B = 1.4 T (BaBar) (all-si)"
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
		fin.open("../../output/"+fnames[f]);
		if(!fin){ cout << "\033[1;31mCouldn't find input file '" << fnames[f] << "'. Bailing out!\033[0m" << endl; exit(0);}
		fin.close();

		Fin[f] = new TFile("../../output/"+fnames[f]);

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
	TGraphErrors *** g_dpp_v_p_et_bins   = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dth_v_p_et_bins   = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dph_v_p_et_bins   = new TGraphErrors ** [size_loaded]; 

	TGraphErrors *** g_dpp_v_et_p_bins   = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dth_v_et_p_bins   = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dph_v_et_p_bins   = new TGraphErrors ** [size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		std::string temp_lab = (std::string)(label[f]);
		double max_dpp = temp_lab.find("1.4") != std::string::npos ? 13.:6.;

		g_dpp_v_p_et_bins  [f] = new TGraphErrors * [num_eta_bin[f]];
		g_dth_v_p_et_bins  [f] = new TGraphErrors * [num_eta_bin[f]];
		g_dph_v_p_et_bins  [f] = new TGraphErrors * [num_eta_bin[f]];

		g_dpp_v_et_p_bins  [f] = new TGraphErrors * [num_mom_bin[f]];
		g_dth_v_et_p_bins  [f] = new TGraphErrors * [num_mom_bin[f]];
		g_dph_v_et_p_bins  [f] = new TGraphErrors * [num_mom_bin[f]];

		for(int et = 0 ; et < num_eta_bin[f] ; et++){
			g_dpp_v_p_et_bins  [f][et] = graph_from_histo( h1_dpp_v_p_et_bins  [f][et] , 55 + (et-num_eta_bin[f]/2)*7 , 20 , 0.0 , max_dpp );
			g_dth_v_p_et_bins  [f][et] = graph_from_histo( h1_dth_v_p_et_bins  [f][et] , 55 + (et-num_eta_bin[f]/2)*7 , 20 , 2e-2, 100 );
			g_dph_v_p_et_bins  [f][et] = graph_from_histo( h1_dph_v_p_et_bins  [f][et] , 55 + (et-num_eta_bin[f]/2)*7 , 20 , 2e-2, 1e4 );
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			g_dpp_v_et_p_bins  [f][p] = graph_from_histo( h1_dpp_v_et_p_bins  [f][p] , 51 + p*5 , 20 , 0.0 , max_dpp );
			g_dth_v_et_p_bins  [f][p] = graph_from_histo( h1_dth_v_et_p_bins  [f][p] , 51 + p*5 , 20 , 2e-2, 100 );
			g_dph_v_et_p_bins  [f][p] = graph_from_histo( h1_dph_v_et_p_bins  [f][p] , 51 + p*5 , 20 , 2e-2, 1e4 );
		}
	}

	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas ** c1 = new TCanvas*[size_loaded];
	TLegend ** leg21 = new TLegend*[size_loaded];
	TLegend ** leg22 = new TLegend*[size_loaded];
	TLegend ** leg23 = new TLegend*[size_loaded];
	TLegend ** leg24 = new TLegend*[size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		c1[f] = new TCanvas(Form("c1_%i",f),Form("c1_%i",f),1200,900);
		c1[f] -> Divide(2,2);

		c1[f] -> cd(1); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetLogy();
		g_dth_v_p_et_bins[f][num_eta_bin[f]/2] -> GetXaxis() -> SetRangeUser(0,30);
		g_dth_v_p_et_bins[f][num_eta_bin[f]/2] -> Draw("APL");
		for(int et = num_eta_bin[f]/2 ; et < num_eta_bin[f] ; et++){
			g_dth_v_p_et_bins[f][et] -> Draw("samePL");
		}

		c1[f] -> cd(2); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetLogy();
		g_dph_v_p_et_bins[f][num_eta_bin[f]/2] -> GetXaxis() -> SetRangeUser(0,30);
		g_dph_v_p_et_bins[f][num_eta_bin[f]/2] -> Draw("APL");
		for(int et = num_eta_bin[f]/2 ; et < num_eta_bin[f] ; et++){
			g_dph_v_p_et_bins[f][et] -> Draw("samePL");
		}

		c1[f] -> cd(3); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetLogy();
		g_dth_v_et_p_bins[f][num_mom_bin[f]-1] -> GetXaxis() -> SetRangeUser(0,4);
		g_dth_v_et_p_bins[f][num_mom_bin[f]-1] -> Draw("APL");
		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			g_dth_v_et_p_bins[f][p] -> Draw("samePL");
		}

		c1[f] -> cd(4); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.15); gPad -> SetLogy();
		g_dph_v_et_p_bins[f][0] -> GetXaxis() -> SetRangeUser(0,4);
		g_dph_v_et_p_bins[f][0] -> Draw("APL");
		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			g_dph_v_et_p_bins[f][p] -> Draw("samePL");
		}

		// ------------
		leg21[f] = new TLegend(0.20,0.6,0.95,0.89);
		leg21[f] -> SetLineColor(0);
		leg21[f] -> SetNColumns(3);
		leg21[f] -> SetTextSize(0.05);
		leg21[f] -> SetHeader(label[f]+", 10 #mum pixel","C");
		leg21[f] -> AddEntry( (TObject*)0, "#eta", "");
		leg21[f] -> AddEntry( (TObject*)0, "" , "");
		leg21[f] -> AddEntry( (TObject*)0, "" , "");
		for(int et = num_eta_bin[f]/2 ; et < num_eta_bin[f]-1 ; et++){
			leg21[f] -> AddEntry( g_dth_v_p_et_bins[f][et] , Form("( %.1f , %.1f )", (*TVT_eta_bin[f])[et] , (*TVT_eta_bin[f])[et+1] ) );
		}
		c1[f] -> cd(1);
		leg21[f] -> Draw("same");
		// ------------
		leg22[f] = new TLegend(0.20,0.6,0.95,0.89);
		leg22[f] -> SetLineColor(0);
		leg22[f] -> SetNColumns(3);
		leg22[f] -> SetTextSize(0.05);
		leg22[f] -> SetHeader(label[f]+", 10 #mum pixel","C");
		leg22[f] -> AddEntry( (TObject*)0, "#eta", "");
		leg22[f] -> AddEntry( (TObject*)0, "" , "");
		leg22[f] -> AddEntry( (TObject*)0, "" , "");
		for(int et = num_eta_bin[f]/2 ; et < num_eta_bin[f]-1 ; et++){
			leg22[f] -> AddEntry( g_dth_v_p_et_bins[f][et] , Form("( %.1f , %.1f )", (*TVT_eta_bin[f])[et] , (*TVT_eta_bin[f])[et+1] ) );
		}
		c1[f] -> cd(2);
		leg22[f] -> Draw("same");
		// ------------
		leg23[f] = new TLegend(0.20,0.6,0.95,0.89);
		leg23[f] -> SetLineColor(0);
		leg23[f] -> SetNColumns(3);
		leg23[f] -> SetTextSize(0.05);
		leg23[f] -> SetHeader(label[f]+", 10 #mum pixel","C");
		leg23[f] -> AddEntry( (TObject*)0, "p [ GeV/#it{c} ]", "");
		leg23[f] -> AddEntry( (TObject*)0, "" , "");
		leg23[f] -> AddEntry( (TObject*)0, "" , "");
		for(int p = 0 ; p < num_mom_bin[0] ; p++){
			leg23[f] -> AddEntry( g_dth_v_et_p_bins[f][p] , Form("( %.0f , %.0f )", (*TVT_mom_bin[f])[p] , (*TVT_mom_bin[f])[p+1] ) );
		}
		c1[f] -> cd(3);
		leg23[f] -> Draw("same");
		// ------------
		leg24[f] = new TLegend(0.20,0.6,0.95,0.89);
		leg24[f] -> SetLineColor(0);
		leg24[f] -> SetNColumns(3);
		leg24[f] -> SetTextSize(0.05);
		leg24[f] -> SetHeader(label[f]+", 10 #mum pixel","C");
		leg24[f] -> AddEntry( (TObject*)0, "p [ GeV/#it{c} ]", "");
		leg24[f] -> AddEntry( (TObject*)0, "" , "");
		leg24[f] -> AddEntry( (TObject*)0, "" , "");
		for(int p = 0 ; p < num_mom_bin[0] ; p++){
			leg24[f] -> AddEntry( g_dph_v_et_p_bins[f][p] , Form("( %.0f , %.0f )", (*TVT_mom_bin[f])[p] , (*TVT_mom_bin[f])[p+1] ) );
		}
		c1[f] -> cd(4);
		leg24[f] -> Draw("same");

		// ------------
		c1[f] -> Modified();
		c1[f] -> Update();
	}
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	TString out_name = "results_ang_res_at_PID_";
	for(int f = 0 ; f < size_loaded ; f++){
		c1[f] -> Print( out_name + Form("%i.pdf",f+1)); 
	}

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
	g_l1->SetMarkerSize(2.2);
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
