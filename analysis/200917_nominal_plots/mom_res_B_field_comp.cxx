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
		"../../output/output_skimmed_pi-_det2_20x20_Beast_FastTrackingEval.root"  ,	//  0
		"../../output/output_skimmed_pi-_det2_10x10_Beast_FastTrackingEval.root"  ,	//  1
		"../../output/output_skimmed_pi-_det2_20x20_sPHENIX_FastTrackingEval.root",	//  2
		"../../output/output_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEval.root"	//  3
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
	double mom_bin[] = {1.,25.};
	const int p_1GeV  = idx_from_vector( 1.,TVT_mom_bin[0]);
	const int p_25GeV = idx_from_vector(25.,TVT_mom_bin[0]);
	int selected_bins[] = {p_1GeV,p_25GeV};
	int size_selected_bins = sizeof(selected_bins)/sizeof(*selected_bins);

	TGraphErrors ** g_dpp_v_et_selected_20um_Beast_si   = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_10um_Beast_si   = new TGraphErrors * [size_selected_bins];
	TGraphErrors ** g_dpp_v_et_selected_20um_sPHENIX_si = new TGraphErrors * [size_selected_bins];
        TGraphErrors ** g_dpp_v_et_selected_10um_sPHENIX_si = new TGraphErrors * [size_selected_bins];

	int color_30[] = {2,97};	int marker_30[] = {24,20};
	int color_15[] = {4,63};	int marker_15[] = {25,21};

    	for(int i = 0 ; i < size_selected_bins ; i++){
		g_dpp_v_et_selected_20um_Beast_si  [i] = graph_from_histo( h1_dpp_v_et_p_bins[0][selected_bins[i]] ,color_30[i],marker_30[i],0,999);
		g_dpp_v_et_selected_10um_Beast_si  [i] = graph_from_histo( h1_dpp_v_et_p_bins[1][selected_bins[i]] ,color_30[i],marker_30[i],0,999);
		g_dpp_v_et_selected_20um_sPHENIX_si[i] = graph_from_histo( h1_dpp_v_et_p_bins[2][selected_bins[i]] ,color_15[i],marker_15[i],0,999);
		g_dpp_v_et_selected_10um_sPHENIX_si[i] = graph_from_histo( h1_dpp_v_et_p_bins[3][selected_bins[i]] ,color_15[i],marker_15[i],0,999);

		g_dpp_v_et_selected_20um_Beast_si[i] -> SetTitle("#pi^{-}, 20x20 #mum pixel tracker");
		g_dpp_v_et_selected_10um_Beast_si[i] -> SetTitle("#pi^{-}, 10x10 #mum pixel tracker");

		g_dpp_v_et_selected_20um_Beast_si[i] -> SetMaximum(20);
		g_dpp_v_et_selected_20um_Beast_si[i] -> GetXaxis() -> SetRangeUser(0,4);
		g_dpp_v_et_selected_10um_Beast_si[i] -> SetMaximum(20);
                g_dpp_v_et_selected_10um_Beast_si[i] -> GetXaxis() -> SetRangeUser(0,4);
	}
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1200,900);
	gPad -> SetRightMargin(0.02); /*gPad -> SetGridx(); gPad -> SetGridy();*/ gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
	g_dpp_v_et_selected_20um_Beast_si  [1] -> Draw(   "APL");
	g_dpp_v_et_selected_20um_sPHENIX_si[1] -> Draw("samePL");
	g_dpp_v_et_selected_20um_Beast_si  [0] -> Draw("samePL");
        g_dpp_v_et_selected_20um_sPHENIX_si[0] -> Draw("samePL");
	// ------------
	TLegend * leg1 = new TLegend(0.35,0.5,0.70,0.89);
	leg1 -> SetLineColor(0);
	leg1 -> AddEntry( g_dpp_v_et_selected_20um_sPHENIX_si[0] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_1GeV ],(*TVT_mom_bin[0])[p_1GeV +1]));
	leg1 -> AddEntry( g_dpp_v_et_selected_20um_sPHENIX_si[1] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_25GeV],(*TVT_mom_bin[0])[p_25GeV+1]));
	leg1 -> AddEntry( (TObject*)0, "", "");
	leg1 -> AddEntry( g_dpp_v_et_selected_20um_Beast_si  [0] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_1GeV ],(*TVT_mom_bin[0])[p_1GeV +1]));
	leg1 -> AddEntry( g_dpp_v_et_selected_20um_Beast_si  [1] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_25GeV],(*TVT_mom_bin[0])[p_25GeV+1]));
	// ------------
	TLatex * tex_30 = new TLatex(0.5,11.3,"3.0 T");
	TLatex * tex_15 = new TLatex(0.5,17.4,"1.4 T");
	// ------------
	TLine * l1 = new TLine(0,1,4,1);	l1 -> SetLineStyle(7);	l1 -> SetLineWidth(2);
	// ------------
        leg1 -> Draw("same");
	tex_30 -> Draw("same");
	tex_15 -> Draw("same");
	l1 -> Draw("same");
	// ------------
	c1 -> Modified();
	c1 -> Update();	
	// ------------------------------------------------------------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1200,900);
        gPad -> SetRightMargin(0.02); /*gPad -> SetGridx(); gPad -> SetGridy();*/ gPad -> SetBottomMargin(0.13); gPad -> SetLeftMargin(0.13);
        g_dpp_v_et_selected_10um_Beast_si  [1] -> Draw(   "APL");
        g_dpp_v_et_selected_10um_sPHENIX_si[1] -> Draw("samePL");
        g_dpp_v_et_selected_10um_Beast_si  [0] -> Draw("samePL");
        g_dpp_v_et_selected_10um_sPHENIX_si[0] -> Draw("samePL");
	// ------------
        TLegend * leg2 = new TLegend(0.35,0.5,0.70,0.89);
        leg2 -> SetLineColor(0);
        leg2 -> AddEntry( g_dpp_v_et_selected_10um_sPHENIX_si[0] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_1GeV ],(*TVT_mom_bin[0])[p_1GeV +1]));
        leg2 -> AddEntry( g_dpp_v_et_selected_10um_sPHENIX_si[1] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_25GeV],(*TVT_mom_bin[0])[p_25GeV+1]));
        leg2 -> AddEntry( (TObject*)0, "", "");
	leg2 -> AddEntry( g_dpp_v_et_selected_10um_Beast_si  [0] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_1GeV ],(*TVT_mom_bin[0])[p_1GeV +1]));
        leg2 -> AddEntry( g_dpp_v_et_selected_10um_Beast_si  [1] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin[0])[p_25GeV],(*TVT_mom_bin[0])[p_25GeV+1]));
        // ------------
	leg2 -> Draw("same");
	tex_30 -> Draw("same");
        tex_15 -> Draw("same");
	l1 -> Draw("same");
	// ------------
	c2 -> Modified();
        c2 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_mom_res_B_field_comp.pdf(");
	c2 -> Print("results_mom_res_B_field_comp.pdf)");
	
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
	g_l1->SetMarkerSize(2.0);
	g_l1->SetLineWidth(3);

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
