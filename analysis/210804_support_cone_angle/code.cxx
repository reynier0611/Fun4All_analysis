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
#include "TLine.h"
#include "TVirtualPad.h"

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max, double text_size );
void pretty_axis_graph( TGraph * g, TString xtit, TString ytit, double pmin, double pmax, TString Bfield , double miny, double maxy, double text_size );
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
		/*
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_36.5_degsigma_eta_14_p_12_.root",
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_20_degsigma_eta_14_p_12_.root",
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_24_degsigma_eta_14_p_12_.root",
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_28_degsigma_eta_14_p_12_.root",
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_32_degsigma_eta_14_p_12_.root",
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_34_degsigma_eta_14_p_12_.root",
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_40_degsigma_eta_14_p_12_.root",
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_42_degsigma_eta_14_p_12_.root"
		*/
		/*
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_36.5_degsigma_eta_35_p_12_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_20_degsigma_eta_35_p_12_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_24_degsigma_eta_35_p_12_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_28_degsigma_eta_35_p_12_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_32_degsigma_eta_35_p_12_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_34_degsigma_eta_35_p_12_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_40_degsigma_eta_35_p_12_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_42_degsigma_eta_35_p_12_.root"
		*/
		/*
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_36.5_degsigma_eta_30_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_20_degsigma_eta_30_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_24_degsigma_eta_30_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_28_degsigma_eta_30_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_32_degsigma_eta_30_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_34_degsigma_eta_30_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_40_degsigma_eta_30_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_42_degsigma_eta_30_p_6_.root"
		*/
		"../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_36.5_degsigma_eta_46_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_20_degsigma_eta_46_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_24_degsigma_eta_46_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_28_degsigma_eta_46_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_32_degsigma_eta_46_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_34_degsigma_eta_46_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_40_degsigma_eta_46_p_6_.root",
                "../../output/output_mom_res_skimmed_simp_geom_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507_FastSimEval_theta_cone_42_degsigma_eta_46_p_6_.root"
	};
	TString labels[] = {
		"#theta_{C} = 36.5^{o}",
		"#theta_{C} = 20.0^{o}",
		"#theta_{C} = 24.0^{o}",
		"#theta_{C} = 28.0^{o}",
		"#theta_{C} = 32.0^{o}",
		"#theta_{C} = 34.0^{o}",
		"#theta_{C} = 40.0^{o}",
		"#theta_{C} = 42.0^{o}"
	};
	double angle_C[] = {36.5,20,24,28,32,34,40,42};
	double eta_C[sizeof(angle_C)/sizeof(*angle_C)] = {0};
	for(int i = 0 ; i < sizeof(angle_C)/sizeof(*angle_C) ; i++) eta_C[i] = -TMath::Log(TMath::Tan(angle_C[i]*3.14159/180./2.));
	double min_for_ratio = 0.55;
	double max_for_ratio = 3.60;
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
	int color [] = { 1, 2,62, 8,93,50,51,13};
	//int marker[] = {20,24,24,24,24,24,24,24};
	int marker[] = {20,1,1,1,1,1,1,1};

	TGraphErrors *** g_dpp_v_p_et_bins = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dpp_v_et_p_bins = new TGraphErrors ** [size_loaded];
	TH1F *** h1_ratio_dpp_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_ratio_dpp_v_et_p_bins = new TH1F ** [size_loaded];
	TGraphErrors *** g_ratio_dpp_v_p_et_bins = new TGraphErrors ** [size_loaded];
        TGraphErrors *** g_ratio_dpp_v_et_p_bins = new TGraphErrors ** [size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		g_dpp_v_p_et_bins[f] = new TGraphErrors * [num_eta_bin[f]];
		g_dpp_v_et_p_bins[f] = new TGraphErrors * [num_mom_bin[f]];
		h1_ratio_dpp_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
		h1_ratio_dpp_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];
		g_ratio_dpp_v_p_et_bins[f] = new TGraphErrors * [num_eta_bin[f]];
		g_ratio_dpp_v_et_p_bins[f] = new TGraphErrors * [num_mom_bin[f]];

		for(int et = 0 ; et < num_eta_bin[f] ; et++){
			g_dpp_v_p_et_bins[f][et] = graph_from_histo( h1_dpp_v_p_et_bins[f][et],color[f],marker[f],0,999,0.06);
			g_dpp_v_p_et_bins[f][et] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[et],(*TVT_eta_bin[0])[et+1]));
			g_dpp_v_p_et_bins[f][et] -> SetMaximum(et>1?3.5:1.4);
		
			h1_ratio_dpp_v_p_et_bins[f][et] = (TH1F*) h1_dpp_v_p_et_bins[f][et] -> Clone();
			h1_ratio_dpp_v_p_et_bins[f][et] -> SetName(Form("h1_ratio_dpp_v_et_p_bis_%i_%i",f,et));
			h1_ratio_dpp_v_p_et_bins[f][et] -> Divide(h1_dpp_v_p_et_bins[0][et]);
			g_ratio_dpp_v_p_et_bins[f][et] = graph_from_histo( h1_ratio_dpp_v_p_et_bins[f][et],color[f],marker[f],min_for_ratio,max_for_ratio,0.09);
			g_ratio_dpp_v_p_et_bins[f][et] -> GetYaxis() -> SetTitle("ratio to "+labels[0]);
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			g_dpp_v_et_p_bins[f][p] = graph_from_histo( h1_dpp_v_et_p_bins[f][p],color[f],marker[f],0,999,0.06);
			g_dpp_v_et_p_bins[f][p] -> SetTitle(Form("%.1f < p < %.1f GeV/#it{c}",(*TVT_mom_bin[0])[p],(*TVT_mom_bin[0])[p+1]));
		
			h1_ratio_dpp_v_et_p_bins[f][p] = (TH1F*) h1_dpp_v_et_p_bins[f][p] -> Clone();
			h1_ratio_dpp_v_et_p_bins[f][p] -> SetName(Form("h1_ratio_dpp_v_et_p_bins_%i_%i",f,p));
			h1_ratio_dpp_v_et_p_bins[f][p] -> Divide(h1_dpp_v_et_p_bins[0][p]);
			g_ratio_dpp_v_et_p_bins[f][p] = graph_from_histo( h1_ratio_dpp_v_et_p_bins[f][p],color[f],marker[f],min_for_ratio,max_for_ratio,0.09);
			g_ratio_dpp_v_et_p_bins[f][p] -> GetYaxis() -> SetTitle("ratio to "+labels[0]);
		}
	}
	// ------------------------------------------------------------------------------
	//TLegend * leg1 = new TLegend(0.35,0.5,0.75,0.89);
	TLegend * leg1 = new TLegend(0.2,0.4,0.6,0.85);
	leg1 -> SetLineColor(0);
	for(int f = 0 ; f < size_loaded ; f++) leg1 -> AddEntry(g_dpp_v_et_p_bins[f][0],labels[f]+Form(" (#eta=%.2f)",eta_C[f]));
	// ------------------------------------------------------------------------------
	TLine * l1 = new TLine(0,1,(*TVT_eta_bin[0])[num_eta_bin[0]],1);
	l1 -> SetLineStyle(2);
	
	TLine ** t_v1 = new TLine*[size_loaded];
	TLine ** t_v2 = new TLine*[size_loaded];
	for(int f = 0 ; f < size_loaded ; f++){
		t_v1[f] = new TLine(eta_C[f],0.00,eta_C[f],1.9);	t_v1[f] -> SetLineStyle(7);	t_v1[f] -> SetLineColor(color[f]);
		t_v2[f] = new TLine(eta_C[f],min_for_ratio,eta_C[f],max_for_ratio);	t_v2[f] -> SetLineStyle(7);	t_v2[f] -> SetLineColor(color[f]);
	}
	// ------------------------------------------------------------------------------
	// Plotting graphs
	// ----------------------------------------------
	TCanvas ** c4 = new TCanvas*[num_mom_bin[0]];
	TVirtualPad ** pad4t = new TVirtualPad * [num_mom_bin[0]];
	TVirtualPad ** pad4b = new TVirtualPad * [num_mom_bin[0]];
	for( int i = 0 ; i < num_mom_bin[0] ; i++ ){
		c4[i] = new TCanvas(Form("c4_%i",i),Form("c4_%i",i),800,900);
		c4[i] -> Divide(1,2);

		pad4t[i] = c4[i] -> cd(1);
		pad4t[i] -> SetPad(Form("pad4t_%i",i),"",0,0.4,1,1,kWhite,0,0);
		gPad -> SetRightMargin(0.02); gPad -> SetBottomMargin(0.0); gPad -> SetLeftMargin(0.18);
		g_dpp_v_et_p_bins[0][i] -> Draw("ALP");
		for(int f = 0 ; f < size_loaded ; f++){
			g_dpp_v_et_p_bins[f][i] -> Draw("sameL");
			//t_v1[f] -> Draw("same");
		}
		pad4b[i] = c4[i] -> cd(2);
		pad4b[i] -> SetPad(Form("pad4b_%i",i),"",0,0,1,0.4,kWhite,0,0);
		gPad -> SetRightMargin(0.02); gPad -> SetTopMargin(0.0); gPad -> SetLeftMargin(0.18); gPad -> SetBottomMargin(0.17);
                g_ratio_dpp_v_et_p_bins[1][i] -> Draw("AL");
		for(int f = 1 ; f < size_loaded ; f++){
                        g_ratio_dpp_v_et_p_bins[f][i] -> Draw("sameL");
			//t_v2[f] -> Draw("same");
		}
		l1 -> Draw("same");

	c4[i] -> cd(1);
	leg1 -> Draw();
	c4[i] -> Modified();
	c4[i] -> Update();
	}
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	TString out_pdf_file_eta = "results_code_dpp_v_eta.pdf";
	for( int i = 0 ; i < num_mom_bin[0] ; i++ ){
		if(i==0)
			c4[i] -> Print(out_pdf_file_eta+"(");
		else if(i==num_mom_bin[0]-1)
			c4[i] -> Print(out_pdf_file_eta+")");
		else
			c4[i] -> Print(out_pdf_file_eta);
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
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max , double text_size){
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
	g_l1->SetLineWidth(3);

	g_l1->GetXaxis()->SetTitle(xax);
	g_l1->GetXaxis()->SetNdivisions(108);
	g_l1->GetXaxis()->SetTitleSize(text_size);
	g_l1->GetXaxis()->SetLabelSize(text_size);
	g_l1->GetXaxis()->CenterTitle();

	g_l1->GetYaxis()->SetTitle(yax);
	g_l1->GetYaxis()->SetNdivisions(108);
	g_l1->GetYaxis()->SetTitleSize(text_size);
	g_l1->GetYaxis()->SetLabelSize(text_size);
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
void pretty_axis_graph( TGraph * g, TString xtit, TString ytit, double pmin, double pmax, TString Bfield , double miny, double maxy , double text_size ){
	g -> SetTitle(Form(Bfield+", %.1f < p < %.1f GeV/#it{c}",pmin,pmax));

	g -> GetXaxis() -> SetTitle(xtit);
	g -> GetXaxis() -> CenterTitle();
	g -> GetXaxis() -> SetNdivisions(108);
	g -> GetXaxis() -> SetLabelSize(text_size);
	g -> GetXaxis() -> SetTitleSize(text_size);

	g -> GetYaxis() -> SetTitle(ytit);
	g -> GetYaxis() -> CenterTitle();
	g -> GetYaxis() -> SetNdivisions(108);
	g -> GetYaxis() -> SetLabelSize(text_size);
	g -> GetYaxis() -> SetTitleSize(text_size);
	g -> SetMinimum( miny );
	g -> SetMaximum( maxy );

	g -> SetMarkerColor(0);
}
