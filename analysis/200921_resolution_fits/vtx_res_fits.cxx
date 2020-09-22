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
void prettyTH1( TH1F * h1 , int color , int marker , float min , float max );
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	//gStyle->SetErrorX(0.0001);
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		"../../output/output_vtx_res_skimmed_combined_vtx_new_pi-_det2_20x20_Beast_FastTrackingEvalsigma_eta_16_p_16_.root",
		"../../output/output_vtx_res_skimmed_combined_vtx_new_pi-_det2_10x10_Beast_FastTrackingEvalsigma_eta_16_p_16_.root"
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
	// Editing histograms
	double max_dvl[] = {10000,8000,4000,1000,200,100,100,60,60,100,100,200,1000,4000,8000,10000};
	double max_dvt[] = { 3000,3000,1500, 800,200,100,100,60,60,100,100,200, 800,1500,3000, 3000};

	for(int i = 0 ; i < num_eta_bin[0] ; i++){
		prettyTH1( h1_dvl_v_p_et_bins[1][i] , 96 , 20 , 0 , max_dvl[i] );
		prettyTH1( h1_dvt_v_p_et_bins[1][i] , 96 , 20 , 0 , max_dvt[i] );

		h1_dvl_v_p_et_bins[1][i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
		h1_dvt_v_p_et_bins[1][i] -> SetTitle(Form("Beast (3.0 T), %.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
	}
	// ------------------------------------------------------------------------------
        // Doing fits
        TF1 ** f_dvl_20um_30T = new TF1*[num_eta_bin[0]];
	TF1 ** f_dvt_20um_30T = new TF1*[num_eta_bin[0]];
	TF1 ** f_dvl_10um_30T = new TF1*[num_eta_bin[0]];
  	TF1 ** f_dvt_10um_30T = new TF1*[num_eta_bin[0]];

        for(int i = 0 ; i < num_eta_bin[0] ; i++){
                f_dvl_20um_30T[i] = new TF1(Form("f_dvl_20um_30T_%i",i),"sqrt(sq([0]/x)+sq([1]))",0,30);
		f_dvt_20um_30T[i] = new TF1(Form("f_dvt_20um_30T_%i",i),"sqrt(sq([0]/x)+sq([1]))",0,30);
                f_dvl_10um_30T[i] = new TF1(Form("f_dvl_10um_30T_%i",i),"sqrt(sq([0]/x)+sq([1]))",0,30);
  		f_dvt_10um_30T[i] = new TF1(Form("f_dvt_10um_30T_%i",i),"sqrt(sq([0]/x)+sq([1]))",0,30);

		//f_dvl_20um_30T[i] -> SetParLimits(1,0,1E+08);
                //f_dvt_20um_30T[i] -> SetParLimits(1,0,1E+08);
                //f_dvl_10um_30T[i] -> SetParLimits(1,0,1E+08);
                //f_dvt_10um_30T[i] -> SetParLimits(1,0,1E+08);

                f_dvl_20um_30T[i] -> SetLineColor(2); f_dvl_20um_30T[i] -> SetLineStyle(3);
		f_dvt_20um_30T[i] -> SetLineColor(2); f_dvt_20um_30T[i] -> SetLineStyle(3);
                f_dvl_10um_30T[i] -> SetLineColor(2); f_dvl_10um_30T[i] -> SetLineStyle(3);
 	 	f_dvt_10um_30T[i] -> SetLineColor(2); f_dvt_10um_30T[i] -> SetLineStyle(3);

                h1_dvl_v_p_et_bins[0][i] -> Fit(Form("f_dvl_20um_30T_%i",i),"R");
                h1_dvt_v_p_et_bins[0][i] -> Fit(Form("f_dvt_20um_30T_%i",i),"R");
		h1_dvl_v_p_et_bins[1][i] -> Fit(Form("f_dvl_10um_30T_%i",i),"R");
        	h1_dvt_v_p_et_bins[1][i] -> Fit(Form("f_dvt_10um_30T_%i",i),"R");
	}
	// ------------------------------------------------------------------------------
	// Preparing legends
	TLegend ** leg_dvl = new TLegend*[num_eta_bin[0]];
        for(int i = 8 ; i < num_eta_bin[0] ; i++){
                leg_dvl[i] = new TLegend(0.35,0.7,0.96,0.89);
                leg_dvl[i] -> SetLineColor(0);
                leg_dvl[i] -> SetHeader("#sigma(DCA_{z}) = #frac{A}{p} #oplus B","C");
                leg_dvl[i] -> AddEntry(h1_dvl_v_p_et_bins[1][i],Form("A = %.1f, B = %.2f",f_dvl_10um_30T[i]->GetParameter(0),abs(f_dvl_10um_30T[i]->GetParameter(1))));
        }
	// -----------------------------------
	TLegend ** leg_dvt = new TLegend*[num_eta_bin[0]];
        for(int i = 8 ; i < num_eta_bin[0] ; i++){
                leg_dvt[i] = new TLegend(0.35,0.7,0.96,0.89);
                leg_dvt[i] -> SetLineColor(0);
                leg_dvt[i] -> SetHeader("#sigma(DCA_{T}) = #frac{A}{p} #oplus B","C");
                leg_dvt[i] -> AddEntry(h1_dvt_v_p_et_bins[1][i],Form("A = %.1f, B = %.2f",f_dvt_10um_30T[i]->GetParameter(0),abs(f_dvt_10um_30T[i]->GetParameter(1))));
        }
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1400,900);
	c1 -> Divide(4,2);
	for(int i = 8 ; i < num_eta_bin[0] ; i++){
		c1 -> cd(i-7);
                gPad -> SetRightMargin(0.01); gPad -> SetLeftMargin(0.24); gPad -> SetBottomMargin(0.17);
		h1_dvl_v_p_et_bins[1][i] -> Draw();
		leg_dvl[i] -> Draw("same");
	}
	c1 -> Modified();
	c1 -> Update();
	// -----------------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1400,900);
        c2 -> Divide(4,2);
        for(int i = 8 ; i < num_eta_bin[0] ; i++){
                c2 -> cd(i-7);
                gPad -> SetRightMargin(0.01); gPad -> SetLeftMargin(0.24); gPad -> SetBottomMargin(0.17);
                h1_dvt_v_p_et_bins[1][i] -> Draw();
        	leg_dvt[i] -> Draw("same");
	}
        c2 -> Modified();
        c2 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_vtx_res_fits.pdf(");
	c2 -> Print("results_vtx_res_fits.pdf)");
	
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
        h1->GetXaxis()->SetTitleSize(0.06);
        h1->GetXaxis()->SetLabelSize(0.06);
        h1->GetXaxis()->CenterTitle();

        h1->GetYaxis()->SetNdivisions(108);
        h1->GetYaxis()->SetTitleSize(0.06);
        h1->GetYaxis()->SetLabelSize(0.06);
        h1->GetYaxis()->CenterTitle();

        if(min!=999) h1 -> SetMinimum(min);
        if(max!=999) h1 -> SetMaximum(max);
}

