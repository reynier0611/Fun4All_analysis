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

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
void prettyTH1( TH1F * h1 , int color , int marker , float min , float max );
double sq(double x){ return x*x;}
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	//gStyle->SetErrorX(0.0001);
	gStyle->SetTitleSize(0.08,"t");
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		//"../../output/output_mom_res_skimmed_pi-_det2_10x10_Beast_FastTrackingEvalsigma_eta_7_p_10_.root"                 ,	//  0
		//"../../output/output_mom_res_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEvalsigma_eta_7_p_10_.root"               	//  1
		"../../output/output_mom_res_skimmed_pi-_det2_10x10_Beast_FastTrackingEvalsigma_eta_8_p_16_.root"                 ,   //  0
                "../../output/output_mom_res_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEvalsigma_eta_8_p_16_.root"                   //  1
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

	TH1F *** h1_dpp_v_p_et_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dth_v_p_et_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dph_v_p_et_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dppT_v_pT_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dpp_v_et_p_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dth_v_et_p_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dph_v_et_p_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dppT_v_et_pT_bins = new TH1F ** [size_loaded];

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

		h1_dpp_v_p_et_bins  [f] = new TH1F * [num_eta_bin[f]];
		h1_dth_v_p_et_bins  [f] = new TH1F * [num_eta_bin[f]];
		h1_dph_v_p_et_bins  [f] = new TH1F * [num_eta_bin[f]];
		h1_dppT_v_pT_et_bins[f] = new TH1F * [num_eta_bin[f]];
		h1_dpp_v_et_p_bins  [f] = new TH1F * [num_mom_bin[f]];
		h1_dth_v_et_p_bins  [f] = new TH1F * [num_mom_bin[f]];
		h1_dph_v_et_p_bins  [f] = new TH1F * [num_mom_bin[f]];
		h1_dppT_v_et_pT_bins[f] = new TH1F * [num_mom_bin[f]];

		for(int et = 0 ; et < num_eta_bin[f] ; et++){
			h1_dpp_v_p_et_bins  [f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dpp_v_p_et_bins_%i"  ,et));
			h1_dth_v_p_et_bins  [f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dth_v_p_et_bins_%i"  ,et));
			h1_dph_v_p_et_bins  [f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dph_v_p_et_bins_%i"  ,et));
			h1_dppT_v_pT_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dppT_v_pT_et_bins_%i",et));
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			h1_dpp_v_et_p_bins  [f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dpp_v_et_p_bins_%i"  ,p ));
			h1_dth_v_et_p_bins  [f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dth_v_et_p_bins_%i"  ,p ));
			h1_dph_v_et_p_bins  [f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dph_v_et_p_bins_%i"  ,p ));
			h1_dppT_v_et_pT_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dppT_v_et_pT_bins_%i",p ));
		}
	}

	cout << "\neta bin boundaries:\n"; for(int et = 0 ; et < num_eta_bin[0]+1 ; et++) cout << (*TVT_eta_bin[0])[et] << ", "; cout << "\n";
	cout << "\np bin boundaries:\n"  ; for(int p  = 0 ; p  < num_mom_bin[0]+1 ; p ++) cout << (*TVT_mom_bin[0])[ p] << ", "; cout << "\n\n";

	// #######################################################################################################################################
        // EDIT THE CODE BELOW DEPENDING ON WHAT YOU WANT TO PLOT
	// ------------------------------------------------------------------------------
	double max[] = {1,1,1,1,2,3,5,10};
	// Editing histograms
	for(int i = 0 ; i < num_eta_bin[0] ; i++){ 
		h1_dpp_v_p_et_bins[0][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dpp_v_p_et_bins[1][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));

		prettyTH1( h1_dpp_v_p_et_bins[0][i] , 96 , 20 , 999 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[1][i] , 62 , 24 , 999 , max[i] );
        
		h1_dppT_v_pT_et_bins[0][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dppT_v_pT_et_bins[1][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));

                prettyTH1( h1_dppT_v_pT_et_bins[0][i] , 96 , 20 , 999 , max[i] );
                prettyTH1( h1_dppT_v_pT_et_bins[1][i] , 62 , 24 , 999 , max[i] );
	}
	// ------------------------------------------------------------------------------
	TF1 ** f_3_0_T_fits = new TF1 * [8];
	TF1 ** f_1_4_T_fits = new TF1 * [8];

	f_3_0_T_fits[0] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.018,2)*x*x+pow(0.382,2))",0,5);
	f_1_4_T_fits[0] = new TF1("f_1_4_T_fits_0","sqrt(pow(0.041,2)*x*x+pow(0.773,2))",0,5);
                                                    
	f_3_0_T_fits[1] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.016,2)*x*x+pow(0.431,2))",0,5);
        f_1_4_T_fits[1] = new TF1("f_1_4_T_fits_7","sqrt(pow(0.034,2)*x*x+pow(0.906,2))",0,5);
                                                  
	f_3_0_T_fits[2] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.016,2)*x*x+pow(0.424,2))",0,5);
        f_1_4_T_fits[2] = new TF1("f_1_4_T_fits_7","sqrt(pow(0.034,2)*x*x+pow(0.922,2))",0,5);
                                                    
	f_3_0_T_fits[3] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.012,2)*x*x+pow(0.462,2))",0,5);
        f_1_4_T_fits[3] = new TF1("f_1_4_T_fits_7","sqrt(pow(0.026,2)*x*x+pow(1.000,2))",0,5);
                                                    
	f_3_0_T_fits[4] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.018,2)*x*x+pow(0.721,2))",0,5);
        f_1_4_T_fits[4] = new TF1("f_1_4_T_fits_7","sqrt(pow(0.041,2)*x*x+pow(1.551,2))",0,5);
                                                    
	f_3_0_T_fits[5] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.039,2)*x*x+pow(1.331,2))",0,5);
        f_1_4_T_fits[5] = new TF1("f_1_4_T_fits_7","sqrt(pow(0.085,2)*x*x+pow(2.853,2))",0,5);
                                                    
	f_3_0_T_fits[6] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.103,2)*x*x+pow(2.441,2))",0,5);
        f_1_4_T_fits[6] = new TF1("f_1_4_T_fits_7","sqrt(pow(0.215,2)*x*x+pow(5.254,2))",0,5);
                                                    
	f_3_0_T_fits[7] = new TF1("f_3_0_T_fits_0","sqrt(pow(0.281,2)*x*x+pow(4.716,2))",0,5);
        f_1_4_T_fits[7] = new TF1("f_1_4_T_fits_7","sqrt(pow(0.642,2)*x*x+pow(9.657,2))",0,5);

	// ------------------------------------------------------------------------------
	// Doing fits

	TF1 ** f_dpp_10um_30T  = new TF1*[num_eta_bin[0]];
	TF1 ** f_dpp_10um_14T  = new TF1*[num_eta_bin[0]];
        TF1 ** f_dppT_10um_30T = new TF1*[num_eta_bin[0]];
        TF1 ** f_dppT_10um_14T = new TF1*[num_eta_bin[0]];

	for(int i = 0 ; i < num_eta_bin[0] ; i++){
		f_dpp_10um_30T[i] = new TF1(Form("f_dpp_10um_30T_%i",i),"pol0",0,5);
                f_dpp_10um_14T[i] = new TF1(Form("f_dpp_10um_14T_%i",i),"pol0",0,5);
		//f_dpp_10um_30T[i] = new TF1(Form("f_dpp_10um_30T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,5);
		//f_dpp_10um_14T[i] = new TF1(Form("f_dpp_10um_14T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,5);
                f_dpp_10um_30T[i] -> SetLineColor(2); f_dpp_10um_30T[i] -> SetLineStyle(3);
                f_dpp_10um_14T[i] -> SetLineColor(4); f_dpp_10um_14T[i] -> SetLineStyle(3);
		h1_dpp_v_p_et_bins[0][i] -> Fit(Form("f_dpp_10um_30T_%i",i),"R");
		h1_dpp_v_p_et_bins[1][i] -> Fit(Form("f_dpp_10um_14T_%i",i),"R");
		
                f_dppT_10um_30T[i] = new TF1(Form("f_dppT_10um_30T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
                f_dppT_10um_14T[i] = new TF1(Form("f_dppT_10um_14T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
                f_dppT_10um_30T[i] -> SetLineColor(2); f_dppT_10um_30T[i] -> SetLineStyle(3);
                f_dppT_10um_14T[i] -> SetLineColor(4); f_dppT_10um_14T[i] -> SetLineStyle(3);
                h1_dppT_v_pT_et_bins[0][i] -> Fit(Form("f_dppT_10um_30T_%i",i),"R");
                h1_dppT_v_pT_et_bins[1][i] -> Fit(Form("f_dppT_10um_14T_%i",i),"R");
	}
	// ------------------------------------------------------------------------------
	// Legend	
	TLegend ** leg = new TLegend*[num_eta_bin[0]];
	for(int i = 0 ; i < num_eta_bin[0] ; i++){
		leg[i] = new TLegend(0.20,0.7,0.89,0.89);
		leg[i] -> SetLineColor(0);	
		leg[i] -> AddEntry(h1_dpp_v_p_et_bins[0][i],Form("3.0T, B = %.3f",f_dpp_10um_30T[i]->GetParameter(0)));
		leg[i] -> AddEntry(h1_dpp_v_p_et_bins[1][i],Form("1.4T, B = %.3f",f_dpp_10um_14T[i]->GetParameter(0)));
	}
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1400,900);
	c1 -> Divide(4,2);

	for(int i = 0 ; i < num_eta_bin[0] ; i++){
		c1 -> cd(i+1);
		gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);
		h1_dpp_v_p_et_bins[0][i] -> Draw(      ); // 10um pixel, 3.0T
		h1_dpp_v_p_et_bins[1][i] -> Draw("same"); // 10um pixel, 1.4T
		leg[i] -> Draw("same");
		f_3_0_T_fits[i] -> Draw("same");
		f_1_4_T_fits[i] -> Draw("same");
	}
	c1 -> Modified();
	c1 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_mom_res_fits.pdf");
	
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
void prettyTH1( TH1F * h1 , int color , int marker , float min , float max ){
	h1->SetLineColor(color);
	h1->SetMarkerColor(color);
	h1->SetMarkerStyle(marker);
	h1->SetMarkerSize(1.4);
	h1->SetLineWidth(1);

	h1->GetXaxis()->SetNdivisions(108);
	h1->GetXaxis()->SetTitleSize(0.07);
	h1->GetXaxis()->SetLabelSize(0.07);
	h1->GetXaxis()->CenterTitle();

	h1->GetYaxis()->SetNdivisions(108);
	h1->GetYaxis()->SetTitleSize(0.07);
	h1->GetYaxis()->SetLabelSize(0.07);
	h1->GetYaxis()->CenterTitle();

	if(min!=999) h1 -> SetMinimum(min);
        if(max!=999) h1 -> SetMaximum(max);
}
