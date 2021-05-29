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
#include "TLatex.h"
#include "TLegend.h"
#include "TVectorT.h"
#include "TGraphErrors.h"

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max );
void pretty_axis_graph( TGraph * g, TString xtit, TString ytit, double pmin, double pmax, TString Bfield , double miny, double maxy );
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetErrorX(0.0001);
	gStyle->SetTitleSize(0.08,"t");
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		"../../output/output_skimmed_pi-_det2_10x10_Beast_FastTrackingEval.root",			
		"../../output/output_skimmed_pi-_det2_both_GEMs_RICH_10x10_Beast_FastTrackingEval.root",        
		"../../output/output_skimmed_combined_pi-_det2_10umGEM_RICH_10x10_Beast_FastTrackingEval.root", 

		"../../output/output_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEval.root",			
                "../../output/output_skimmed_pi-_det2_both_GEMs_RICH_10x10_sPHENIX_FastTrackingEval.root",	
		"../../output/output_skimmed_combined_pi-_det2_10umGEM_RICH_10x10_sPHENIX_FastTrackingEval.root"
	
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
	int color [] = { 1, 2,62, 1, 2,62};
        int marker[] = {20,21,22,20,21,22};

        TGraphErrors *** g_dpp_v_p_et_bins = new TGraphErrors ** [size_loaded];
        TGraphErrors *** g_dpp_v_et_p_bins = new TGraphErrors ** [size_loaded];

        for(int f = 0 ; f < size_loaded ; f++){
                g_dpp_v_p_et_bins[f] = new TGraphErrors * [num_eta_bin[f]];
                g_dpp_v_et_p_bins[f] = new TGraphErrors * [num_mom_bin[f]];

                for(int et = 0 ; et < num_eta_bin[f] ; et++){
                        g_dpp_v_p_et_bins[f][et] = graph_from_histo( h1_dpp_v_p_et_bins[f][et],color[f],marker[f],0,999);
                        g_dpp_v_p_et_bins[f][et] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[et],(*TVT_eta_bin[0])[et+1]));
                }

                for(int p = 0 ; p < num_mom_bin[f] ; p++){
                        g_dpp_v_et_p_bins[f][p] = graph_from_histo( h1_dpp_v_et_p_bins[f][p],color[f],marker[f],0,999);
                        g_dpp_v_et_p_bins[f][p] -> SetTitle(Form("%.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[p],(*TVT_mom_bin[0])[p+1]));
                }
        }
	// ------------------------------------------------------------------------------
	TLegend * leg1 = new TLegend(0.35,0.5,0.75,0.89);
        leg1 -> SetLineColor(0);
        leg1 -> AddEntry( g_dpp_v_et_p_bins[0][0] , "All-Si (10 #mum)"                          );
	leg1 -> AddEntry( g_dpp_v_et_p_bins[1][0] , "All-Si (10 #mum) + GEM (#sigma = 50 #mum)" );
	leg1 -> AddEntry( g_dpp_v_et_p_bins[2][0] , "All-Si (10 #mum) + Si disk (10 #mum)"      );
	// ------------------------------------------------------------------------------
	// Plotting graphs
	// ----------------------------------------------
	const int bin_1 = 9;
	const int bin_2 = 6;
	// ----------------------------------------------
	TCanvas * c1 = new TCanvas("c1","c1",1100,900);
	gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
	g_dpp_v_et_p_bins[0][bin_1] -> Draw("ALP");
	for(int f = 0 ; f < size_loaded/2 ; f++){
		g_dpp_v_et_p_bins[f][bin_1] -> Draw("sameLP");
	}
	leg1 -> Draw("same");
	c1 -> Modified();
	c1 -> Update();
	// ----------------------------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1100,900);
        gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
        g_dpp_v_et_p_bins[size_loaded/2][bin_1] -> Draw("ALP");
        for(int f = size_loaded/2 ; f < size_loaded ; f++){
        	g_dpp_v_et_p_bins[f][bin_1] -> Draw("sameLP");
        } 
        c2 -> Modified();
        c2 -> Update();
	// ----------------------------------------------
	double xDum1[] = { 1.5, 3.6};
	double xDum2[] = {-3.6,-1.5};
	double yDum [] = { 0.0, 17.};
	TGraph * gDum_Forw_Beast = new TGraph(2.,xDum1,yDum);
	TGraph * gDum_Back_Beast = new TGraph(2.,xDum2,yDum);
	TGraph * gDum_Forw_BaBar = new TGraph(2.,xDum1,yDum);
        TGraph * gDum_Back_BaBar = new TGraph(2.,xDum2,yDum);

	pretty_axis_graph( gDum_Forw_Beast, "#eta" , "dp/p [%] (Forward)"  , (*TVT_mom_bin[0])[bin_1] , (*TVT_mom_bin[0])[bin_1+1] , "Beast (3.0 T)" , 0 , 10 );
	pretty_axis_graph( gDum_Back_Beast, "#eta" , "dp/p [%] (Backward)" , (*TVT_mom_bin[0])[bin_2] , (*TVT_mom_bin[0])[bin_2+1] , "Beast (3.0 T)" , 0 , 10 );
	pretty_axis_graph( gDum_Forw_BaBar, "#eta" , "dp/p [%] (Forward)"  , (*TVT_mom_bin[0])[bin_1] , (*TVT_mom_bin[0])[bin_1+1] , "BaBar (1.4 T)" , 0 , 22 );
	pretty_axis_graph( gDum_Back_BaBar, "#eta" , "dp/p [%] (Backward)" , (*TVT_mom_bin[0])[bin_2] , (*TVT_mom_bin[0])[bin_2+1] , "BaBar (1.4 T)" , 0 , 22 );
	// ----------------------------------------------
	TLatex * tex_Forw_Beast = new TLatex(  3.5, 2.5 , "z = 300 cm" );
	TLatex * tex_Back_Beast = new TLatex( -1.6, 2.5 , "z =-180 cm" );
	TLatex * tex_Forw_BaBar = new TLatex(  3.5, 5.0 , "z = 300 cm" );
	TLatex * tex_Back_BaBar = new TLatex( -1.6, 5.0 , "z =-180 cm" );
	// ----------------------------------------------
	TCanvas * c3 = new TCanvas("c3","c3",1100,900);
	c3 -> Divide(1,2);
	c3 -> cd(1); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
       	gDum_Forw_Beast -> Draw("AP");
        for(int f = 0 ; f < size_loaded/2 ; f++){
                g_dpp_v_et_p_bins[f][bin_1] -> Draw("sameLP");
        }
	leg1 -> Draw("same");
	tex_Forw_Beast -> Draw("same");
	c3 -> cd(2); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.15); gPad -> SetLeftMargin(0.13);
	gDum_Back_Beast -> Draw("AP");
       	for(int f = 0 ; f < size_loaded/2 ; f++){
                g_dpp_v_et_p_bins[f][bin_2] -> Draw("sameLP");
        }
	tex_Back_Beast -> Draw("same");
        c3 -> Modified();
        c3 -> Update();
	// ----------------------------------------------
        TCanvas * c4 = new TCanvas("c4","c4",1100,900);
        c4 -> Divide(1,2);
        c4 -> cd(1); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
        gDum_Forw_BaBar -> Draw("AP");
        for(int f = size_loaded/2 ; f < size_loaded ; f++){
                g_dpp_v_et_p_bins[f][bin_1] -> Draw("sameLP");
        }
	leg1 -> Draw("same");
	tex_Forw_BaBar -> Draw("same");
        c4 -> cd(2); gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.15); gPad -> SetLeftMargin(0.13);
       	gDum_Back_BaBar -> Draw("AP");
        for(int f = size_loaded/2 ; f < size_loaded ; f++){
                g_dpp_v_et_p_bins[f][bin_2] -> Draw("sameLP");
        }
	tex_Back_BaBar -> Draw("same");
        c4 -> Modified();
        c4 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c3 -> Print("results_newest_code.pdf(");
	c4 -> Print("results_newest_code.pdf)");

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
	//g_l1->SetFillColorAlpha(color,0.2);
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
// ============================================================================================================================================
void pretty_axis_graph( TGraph * g, TString xtit, TString ytit, double pmin, double pmax, TString Bfield , double miny, double maxy ){
	g -> SetTitle(Form(Bfield+", %.1f < p < %.1f GeV/#it{c}",pmin,pmax));

        g -> GetXaxis() -> SetTitle(xtit);
        g -> GetXaxis() -> CenterTitle();
        g -> GetXaxis() -> SetNdivisions(108);
        g -> GetXaxis() -> SetLabelSize(0.07);
        g -> GetXaxis() -> SetTitleSize(0.07);
	
        g -> GetYaxis() -> SetTitle(ytit);
        g -> GetYaxis() -> CenterTitle();
        g -> GetYaxis() -> SetNdivisions(108);
        g -> GetYaxis() -> SetLabelSize(0.07);
        g -> GetYaxis() -> SetTitleSize(0.07);
	g -> SetMinimum( miny );
	g -> SetMaximum( maxy );

	g -> SetMarkerColor(0);
}
