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
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max );
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	//gStyle->SetErrorX(0.0001);
	gStyle->SetTitleSize(0.07,"t");
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		// I was previously using files eta_16_p_10
		"../../output/output_mom_res_skimmed_pi-_det2_10x10_Beast_FastTrackingEvalsigma_eta_8_p_15_.root",
		"../../output/output_mom_res_skimmed_mom_pi+_det2_Beastsigma_eta_8_p_15_.root",
		"../../output/output_mom_res_skimmed_mom_e-_det2_Beastsigma_eta_8_p_15_.root",
		"../../output/output_mom_res_skimmed_mom_mu-_det2_Beastsigma_eta_8_p_15_.root",
		"../../output/output_mom_res_skimmed_mom_proton_det2_Beastsigma_eta_8_p_15_.root"
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
	double max[] = {.7,.7,.7,.7,1,1.9,6,13,29};
	int color[] = {2,94,8,62,51};
	int marker[] = {21,20,22,29,23};

	TGraphErrors *** g_dpp_v_p_et_bins = new TGraphErrors ** [size_loaded];
	for(int f = 0 ; f < size_loaded ; f++){
		g_dpp_v_p_et_bins[f] = new TGraphErrors * [num_eta_bin[0]];
		for(int et = 0 ; et < num_eta_bin[0] ; et++){
			g_dpp_v_p_et_bins[f][et] = graph_from_histo( h1_dpp_v_p_et_bins[f][et] , color[f] , marker[f] , 0 ,max[et]); 
			g_dpp_v_p_et_bins[f][et] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[et],(*TVT_eta_bin[0])[et+1]));
		}
	}

	TGraphErrors *** g_dpp_v_et_p_bins = new TGraphErrors ** [size_loaded];
	for(int f = 0 ; f < size_loaded ; f++){
		g_dpp_v_et_p_bins[f] = new TGraphErrors * [num_mom_bin[0]];
                for(int p = 0 ; p < num_mom_bin[0] ; p++){
			g_dpp_v_et_p_bins[f][p] = graph_from_histo( h1_dpp_v_et_p_bins[f][p] , color[f] , marker[f] , 0 , 999 );
			g_dpp_v_et_p_bins[f][p] -> SetTitle(Form("%.1f < #it{p} < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[p],(*TVT_mom_bin[0])[p+1]));
		}
	}

	// Editing histograms
	for(int i = 0 ; i < num_eta_bin[0] ; i++){ 
		h1_dpp_v_p_et_bins[0][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
		h1_dpp_v_p_et_bins[1][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
		h1_dpp_v_p_et_bins[2][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
		h1_dpp_v_p_et_bins[3][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1])); 

		prettyTH1( h1_dpp_v_p_et_bins[0][i] ,  2 , 21 , 0 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[1][i] , 94 , 20 , 0 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[2][i] ,  8 , 22 , 0 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[3][i] , 62 , 29 , 0 , max[i] );

		h1_dppT_v_pT_et_bins[0][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
		h1_dppT_v_pT_et_bins[1][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
		h1_dppT_v_pT_et_bins[2][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
		h1_dppT_v_pT_et_bins[3][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));

		prettyTH1( h1_dppT_v_pT_et_bins[0][i] ,  2 , 21 , 0 , max[i] );
		prettyTH1( h1_dppT_v_pT_et_bins[1][i] , 94 , 20 , 0 , max[i] );
		prettyTH1( h1_dppT_v_pT_et_bins[2][i] ,  8 , 22 , 0 , max[i] );
		prettyTH1( h1_dppT_v_pT_et_bins[3][i] , 62 , 29 , 0 , max[i] );
	}
	// ------------------------------------------------------------------------------
	// Legend
	TLegend * leg = new TLegend(0.3,0.2,0.93,0.5);
	leg -> SetLineColor(0);
	leg -> SetNColumns(2);
	leg -> AddEntry(g_dpp_v_p_et_bins[0][0],"#pi-");
	leg -> AddEntry(g_dpp_v_p_et_bins[1][0],"#pi+");
	leg -> AddEntry(g_dpp_v_p_et_bins[2][0],"e-"  );
	leg -> AddEntry(g_dpp_v_p_et_bins[3][0],"#mu-");
	leg -> AddEntry(g_dpp_v_p_et_bins[4][0],"p"   );
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1400,900);
	c1 -> Divide(4,2);
	for(int i = 0 ; i < num_eta_bin[0] ; i++){
		c1 -> cd(i+1); 
		gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);
		g_dpp_v_p_et_bins[0][i] -> Draw("APL"   );
		g_dpp_v_p_et_bins[1][i] -> Draw("PLsame");
		g_dpp_v_p_et_bins[2][i] -> Draw("PLsame");
		g_dpp_v_p_et_bins[3][i] -> Draw("PLsame");
		g_dpp_v_p_et_bins[4][i] -> Draw("PLsame");
	}
	c1 -> cd(1);
	leg -> Draw("same");
	c1 -> Modified();
	c1 -> Update();
	// ---------------------------
	TCanvas * c2 = new TCanvas("c2","c2",800,600);
	gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.14); gPad -> SetBottomMargin(0.14);
	g_dpp_v_p_et_bins[0][0] -> Draw("APL"   );
	g_dpp_v_p_et_bins[1][0] -> Draw("PLsame");
	g_dpp_v_p_et_bins[2][0] -> Draw("PLsame");
	g_dpp_v_p_et_bins[3][0] -> Draw("PLsame");
	g_dpp_v_p_et_bins[4][0] -> Draw("PLsame");
	leg -> Draw("same");
	c2 -> Modified();
	c2 -> Update();
	// ---------------------------
        TCanvas * c3 = new TCanvas("c3","c3",800,600);
        gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.14); gPad -> SetBottomMargin(0.14);
	g_dpp_v_et_p_bins[0][5] -> Draw("APL"   );
	g_dpp_v_et_p_bins[1][5] -> Draw("PLsame");
	g_dpp_v_et_p_bins[2][5] -> Draw("PLsame");
	g_dpp_v_et_p_bins[3][5] -> Draw("PLsame");
	g_dpp_v_et_p_bins[4][5] -> Draw("PLsame");
	c3 -> Modified();
        c3 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_particle_comp_c1.pdf");
	c2 -> Print("results_particle_comp_c2.pdf");
	c3 -> Print("results_particle_comp_c3.pdf");

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
	//g_l1->SetMarkerSize(2.2);
	g_l1->SetMarkerSize(1.6);
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

