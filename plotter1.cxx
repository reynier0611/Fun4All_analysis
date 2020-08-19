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

		num_eta_bin[f] = (*TVT_eta_bin[f]).NonZeros();
		num_mom_bin[f] = (*TVT_mom_bin[f]).NonZeros();

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
	// ------------------------------------------------------------------------------
	// Copying and editing histograms
	const int p_25GeV = idx_from_vector(25.,TVT_mom_bin[0]);

	TH1F * h1_dpp_v_et_p25GeV_20um_Beast_si     = (TH1F*) h1_dpp_v_et_p_bins[0][p_25GeV] -> Clone();
	TH1F * h1_dpp_v_et_p25GeV_10um_Beast_si     = (TH1F*) h1_dpp_v_et_p_bins[1][p_25GeV] -> Clone();
	TH1F * h1_dpp_v_et_p25GeV_20um_Beast_si_GEM = (TH1F*) h1_dpp_v_et_p_bins[2][p_25GeV] -> Clone();
	TH1F * h1_dpp_v_et_p25GeV_10um_Beast_si_GEM = (TH1F*) h1_dpp_v_et_p_bins[3][p_25GeV] -> Clone();

	prettyTH1F(h1_dpp_v_et_p25GeV_20um_Beast_si    ,62,21,999,999);
	prettyTH1F(h1_dpp_v_et_p25GeV_10um_Beast_si    ,94,21,999,999);
	prettyTH1F(h1_dpp_v_et_p25GeV_20um_Beast_si_GEM, 4,20,999,999);
	prettyTH1F(h1_dpp_v_et_p25GeV_10um_Beast_si_GEM, 2,20,999,999);

	TGraph * g_dpp_v_et_p25GeV_20um_Beast_si     = graph_from_histo( h1_dpp_v_et_p25GeV_20um_Beast_si     );
	TGraph * g_dpp_v_et_p25GeV_10um_Beast_si     = graph_from_histo( h1_dpp_v_et_p25GeV_10um_Beast_si     );
	TGraph * g_dpp_v_et_p25GeV_20um_Beast_si_GEM = graph_from_histo( h1_dpp_v_et_p25GeV_20um_Beast_si_GEM );
	TGraph * g_dpp_v_et_p25GeV_10um_Beast_si_GEM = graph_from_histo( h1_dpp_v_et_p25GeV_10um_Beast_si_GEM );

	// ------------------------------------------------------------------------------
	TCanvas * c1 = new TCanvas("c1","c1",1200,900);
	gPad -> SetTopMargin(0.02); gPad -> SetRightMargin(0.02);
	h1_dpp_v_et_p25GeV_20um_Beast_si     -> Draw(      );	g_dpp_v_et_p25GeV_20um_Beast_si     -> Draw("sameL");
	h1_dpp_v_et_p25GeV_10um_Beast_si     -> Draw("same");	g_dpp_v_et_p25GeV_10um_Beast_si     -> Draw("sameL");
	h1_dpp_v_et_p25GeV_20um_Beast_si_GEM -> Draw("same");	g_dpp_v_et_p25GeV_20um_Beast_si_GEM -> Draw("sameL");
	h1_dpp_v_et_p25GeV_10um_Beast_si_GEM -> Draw("same");	g_dpp_v_et_p25GeV_10um_Beast_si_GEM -> Draw("sameL");
	// ----------------
	TLegend * leg1 = new TLegend(0.3,0.4,0.7,0.9);
	leg1 -> SetLineColor(0);
	leg1 ->SetHeader(Form("%.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[p_25GeV],(*TVT_mom_bin[0])[p_25GeV+1]),"C");
	leg1 -> AddEntry( h1_dpp_v_et_p25GeV_20um_Beast_si     , "All-Si (20 #mum)"       );
	leg1 -> AddEntry( h1_dpp_v_et_p25GeV_10um_Beast_si     , "All-Si (10 #mum)"       );
	leg1 -> AddEntry( h1_dpp_v_et_p25GeV_20um_Beast_si_GEM , "All-Si (20 #mum) + GEM" );
	leg1 -> AddEntry( h1_dpp_v_et_p25GeV_10um_Beast_si_GEM , "All-Si (10 #mum) + GEM" );
	leg1 -> Draw("same");
	c1 -> Modified();
	c1 -> Update();

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
	h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks

	h1 -> SetMinimum(0.001);
}
// ============================================================================================================================================
int idx_from_vector( double value , TVectorT<double> * vec ){
	int size_vec = (*vec).NonZeros();
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
