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
	gStyle->SetTitleSize(0.06,"t");
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_Beast_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_sPHENIX_FastSimEvalsigma_eta_6_p_15_.root",
	    	"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_92cm_50um_Beast_FastSimEvalsigma_eta_6_p_15_.root",
               	"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_92cm_50um_sPHENIX_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_84cm_GEM_R_92cm_50um_Beast_FastSimEvalsigma_eta_6_p_15_.root",
                "../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_84cm_GEM_R_92cm_50um_sPHENIX_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_60cm_50um_Beast_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_60cm_50um_sPHENIX_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_49cm_GEM_R_60cm_50um_Beast_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_49cm_GEM_R_60cm_50um_sPHENIX_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_barrelR_75cm_Beast_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_barrelR_75cm_sPHENIX_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_DIRC_GEM_Beast_FastSimEvalsigma_eta_6_p_15_.root",
		"../../output/output_mom_res_skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_DIRC_GEM_sPHENIX_FastSimEvalsigma_eta_6_p_15_.root"
	};

	bool plotme[] = {
		true,	// 3.0T - all-si tracker only (v,b,d) = (0.05,0.55,0.24) % X0
		true,	// 1.4T - all-si tracker only (v,b,d) = (0.05,0.55,0.24) % X0
		false,	// 3.0T - all-si + GEM (sigma=50um) @ R = 92 cm
		false,	// 1.4T - all-si + GEM (sigma=50um) @ R = 92 cm
		true,	// 3.0T - all-si + DIRC + GEM (sigma=50um) @ R = 92 cm
		true,	// 1.4T - all-si + DIRC + GEM (sigma=50um) @ R = 92 cm
		false,	// 3.0T - all-si + GEM (sigma=50um) @ R = 60 cm
		false,	// 1.4T - all-si + GEM (sigma=50um) @ R = 60 cm
		false,	// 3.0T - all-si + small radius DIRC + GEM (sigma=50um) @ R = 60 cm
		false,	// 1.4T - all-si + small radius DIRC + GEM (sigma=50um) @ R = 60 cm
		false,	// 3.0T - all-si with outer two layers centered at R = 75 cm and middle two layers centered at R = 75/2 cm
		false,	// 1.4T - all-si with outer two layers centered at R = 75 cm and middle two layers centered at R = 75/2 cm
		true,	// 3.0T - all-si with outer two layers replaced with a GEM (sigma=50um) + DIRC + GEM (sigma=50um) @ R = 92 cm
		true	// 1.4T - all-si with outer two layers replaced with a GEM (sigma=50um) + DIRC + GEM (sigma=50um) @ R = 92 cm
	};

	TString labels[] = {
		"all-si only (B = 3.0 T)",
		"all-si only (B = 1.4 T)",
		"all-si + GEM (R=92 cm, #sigma=50 #mum)",
                "all-si + GEM (R=92 cm, #sigma=50 #mum)",
		"all-si + DIRC (R=84 cm) + GEM (R=92 cm, #sigma=50 #mum)",
		"all-si + DIRC (R=84 cm) + GEM (R=92 cm, #sigma=50 #mum)",
		"all-si + GEM (R=60 cm, #sigma=50 #mum)",
		"all-si + GEM (R=60 cm, #sigma=50 #mum)",
		"all-si + DIRC (R=49 cm) + GEM (R=60 cm, #sigma=50 #mum)",
		"all-si + DIRC (R=49 cm) + GEM (R=60 cm, #sigma=50 #mum)",
		"all-si (barrel R = 75 cm)",
		"all-si (barrel R = 75 cm)",
		"all-si outer two layers replaced with GEM (#sigma=50 #mum) + DIRC + GEM (#sigma=50 #mum) R = 92 cm",
		"all-si outer two layers replaced with GEM (#sigma=50 #mum) + DIRC + GEM (#sigma=50 #mum) R = 92 cm"
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
	int color[] = {1,1,2,2,62,62,8,8,94,94,50,50,210,210};
	int marker[] = {20,20,21,21,22,22,21,21,22,22,23,23,24,24};

	TGraphErrors *** g_dpp_v_p_et_bins   = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dppT_v_pT_et_bins = new TGraphErrors ** [size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		double max = f%2==0 ? 0.9 : 1.9;
		g_dpp_v_p_et_bins  [f] = new TGraphErrors * [num_eta_bin[0]];
		g_dppT_v_pT_et_bins[f] = new TGraphErrors * [num_eta_bin[0]];
		for(int p = 0 ; p < num_eta_bin[f] ; p++){
			g_dpp_v_p_et_bins  [f][p] = graph_from_histo( h1_dpp_v_p_et_bins  [f][p] ,color[f],marker[f],0.,max);
			g_dppT_v_pT_et_bins[f][p] = graph_from_histo( h1_dppT_v_pT_et_bins[f][p] ,color[f],marker[f],0.,max);
			g_dpp_v_p_et_bins  [f][p] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[p],(*TVT_eta_bin[0])[p+1]));
			g_dppT_v_pT_et_bins[f][p] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[p],(*TVT_eta_bin[0])[p+1]));
		}
	}

	// Ratios (different configurations to all-si tracker standalone)
	TH1F *** h1_dpp_v_p_et_bins_i_over_0 = new TH1F ** [size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		h1_dpp_v_p_et_bins_i_over_0[f] = new TH1F * [num_eta_bin[0]];
		int idx = f%2==0 ? 0 : 1;
		for(int p = 0 ; p < num_eta_bin[f] ; p++){
			h1_dpp_v_p_et_bins_i_over_0[f][p] = (TH1F*) h1_dpp_v_p_et_bins[f][p] -> Clone();
			h1_dpp_v_p_et_bins_i_over_0[f][p] -> SetName(Form("h1_dpp_v_p_et_bins_%i_over_0_%i",f,p));
			h1_dpp_v_p_et_bins_i_over_0[f][p] -> Divide(h1_dpp_v_p_et_bins[idx][p]);

			prettyTH1( h1_dpp_v_p_et_bins_i_over_0[f][p] , color[f] , marker[f] , 0.4 , 1.2 );
			h1_dpp_v_p_et_bins_i_over_0[f][p] -> GetYaxis() -> SetTitle("ratio to all-si tracker only");
			h1_dpp_v_p_et_bins_i_over_0[f][p] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[p],(*TVT_eta_bin[0])[p+1]));
		}
	} 

	// Ratios (B=1.4T over B=3.0T)
	TH1F *** h1_dpp_v_p_et_bins_1_4_over_3_0 = new TH1F ** [size_loaded/2];

	for(int f = 0 ; f < size_loaded/2 ; f++){
		h1_dpp_v_p_et_bins_1_4_over_3_0[f] = new TH1F * [num_eta_bin[0]];
		for(int p = 0 ; p < num_eta_bin[f] ; p++){
			h1_dpp_v_p_et_bins_1_4_over_3_0[f][p] = (TH1F*) h1_dpp_v_p_et_bins[2*f+1][p] -> Clone();
			h1_dpp_v_p_et_bins_1_4_over_3_0[f][p] -> SetName(Form("h1_dpp_v_p_et_bins_Brat_%i_%i",f,p));
			h1_dpp_v_p_et_bins_1_4_over_3_0[f][p] -> Divide(h1_dpp_v_p_et_bins[2*f][p]);

			prettyTH1( h1_dpp_v_p_et_bins_1_4_over_3_0[f][p] , color[2*f] , marker[2*f] , 0. , 3.0 );
			h1_dpp_v_p_et_bins_1_4_over_3_0[f][p] -> GetYaxis() -> SetTitle("ratio (1.4T)/(3.0T)");
			h1_dpp_v_p_et_bins_1_4_over_3_0[f][p] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[p],(*TVT_eta_bin[0])[p+1]));
		}
	}

	// ------------------------------------------------------------------------------
	// Legend
	TLegend * leg1 = new TLegend(0.3,0.18,0.95,0.4);
        leg1 -> SetLineColor(0);
        for(int j = 0 ; j < size_loaded ; j+=2){
                if(plotme[j])
			leg1 -> AddEntry(g_dpp_v_p_et_bins[j][0],labels[j]);
        }

	TLegend * leg2 = new TLegend(0.3,0.18,0.95,0.4);
        leg2 -> SetLineColor(0);
        for(int j = 1 ; j < size_loaded ; j+=2){
                if(plotme[j])
			leg2 -> AddEntry(g_dpp_v_p_et_bins[j][0],labels[j]);
        }

	TLegend * leg3 = new TLegend(0.,0.,1.,1.);
	leg3 -> SetLineColor(0);
	for(int j = 0 ; j < size_loaded ; j+=2){	
		if(plotme[j])
			leg3 -> AddEntry(g_dpp_v_p_et_bins[j][0],labels[j]);
	}

	TLegend * leg4 = new TLegend(0.,0.,1.,1.);
        leg4 -> SetLineColor(0);
        for(int j = 1 ; j < size_loaded ; j+=2){
                if(plotme[j])
			leg4 -> AddEntry(g_dpp_v_p_et_bins[j][0],labels[j]);
        }
	// ------------------------------------------------------------------------------
        // Physics Working Group Requirements
        // https://docs.google.com/spreadsheets/d/1ynU7Cu7NlwRvMtbtFdlp_B5xXkw8yBAtWJbenMf-P3U/edit#gid=368031287
        // https://wiki.bnl.gov/eicug/index.php/Yellow_Report_Physics_Common
        TF1 * f_PWG_req_eta_lt_1 = new TF1("f_PWG_req_eta_lt_1","sqrt(pow(0.05*x,2)+pow(0.5,2))",0,30);
	f_PWG_req_eta_lt_1 -> SetLineColor(11);
	f_PWG_req_eta_lt_1 -> SetLineStyle(2);
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c0 = new TCanvas("c0","c0",1400,900);
        gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.14); gPad -> SetBottomMargin(0.17);
        g_dpp_v_p_et_bins[0][0] -> Draw("APE3");
        for(int j = 2 ; j < size_loaded ; j+=2){
                if(plotme[j])
			g_dpp_v_p_et_bins[j][0] -> Draw("samePE3");
        }
        leg1 -> Draw("same");
	f_PWG_req_eta_lt_1 -> Draw("same");
        c0 -> Modified();
        c0 -> Update();
	// ----------------------------------
	TCanvas * c1 = new TCanvas("c1","c1",1400,900);
        gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.14); gPad -> SetBottomMargin(0.17);
        g_dpp_v_p_et_bins[1][0] -> Draw("APE3");
        for(int j = 3 ; j < size_loaded ; j+=2){
                if(plotme[j])
			g_dpp_v_p_et_bins[j][0] -> Draw("samePE3");
        }
        leg2 -> Draw("same");
        f_PWG_req_eta_lt_1 -> Draw("same");
        c1 -> Modified();
        c1 -> Update();
	// ----------------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1400,900);
	c2 -> Divide(3,2);
	for(int i = 0 ; i < num_eta_bin[0]-1 ; i++){
		c2 -> cd(i+1); gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);	
		g_dpp_v_p_et_bins[0][i] -> Draw("APE3");
		for(int j = 2 ; j < size_loaded ; j+=2){
			if(plotme[j])
				g_dpp_v_p_et_bins[j][i] -> Draw("samePE3");
		}
		f_PWG_req_eta_lt_1 -> Draw("same");
	}
	c2 -> cd(6);
	leg3 -> Draw("same");
	c2 -> Modified();
	c2 -> Update();
	// ----------------------------------
	TCanvas * c3 = new TCanvas("c3","c3",1400,900);
        c3 -> Divide(3,2);
        for(int i = 0 ; i < num_eta_bin[0]-1 ; i++){
                c3 -> cd(i+1); gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17); 
                g_dpp_v_p_et_bins[1][i] -> Draw("APE3");
                for(int j = 3 ; j < size_loaded ; j+=2){
                	if(plotme[j])
		        	g_dpp_v_p_et_bins[j][i] -> Draw("samePE3");
                }
		f_PWG_req_eta_lt_1 -> Draw("same");
        }
        c3 -> cd(6);
        leg4 -> Draw("same");
        c3 -> Modified();
        c3 -> Update();
	// ----------------------------------
	TCanvas * c4 = new TCanvas("c4","c4",1400,900);
        c4 -> Divide(3,2);
	for(int i = 0 ; i < num_eta_bin[0]-1 ; i++){
                c4 -> cd(i+1); gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);
                h1_dpp_v_p_et_bins_i_over_0[0][i] -> Draw();
                for(int j = 2 ; j < size_loaded ; j+=2){
                        if(plotme[j])
				h1_dpp_v_p_et_bins_i_over_0[j][i] -> Draw("same");
                }
        }
        c4 -> cd(6);
        leg3 -> Draw("same");
	c4 -> Modified();
        c4 -> Update();
	// ----------------------------------
	TCanvas * c5 = new TCanvas("c5","c5",1400,900);
        c5 -> Divide(3,2);
        for(int i = 0 ; i < num_eta_bin[0]-1 ; i++){
                c5 -> cd(i+1); gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);
                h1_dpp_v_p_et_bins_i_over_0[1][i] -> Draw();
                for(int j = 3 ; j < size_loaded ; j+=2){
			if(plotme[j])
                        	h1_dpp_v_p_et_bins_i_over_0[j][i] -> Draw("same");
                }
        }
        c5 -> cd(6);
        leg4 -> Draw("same");
	c5 -> Modified();
        c5 -> Update();
	// ----------------------------------
        TCanvas * c6 = new TCanvas("c6","c6",1400,900);
        c6 -> Divide(3,2);
	for(int i = 0 ; i < num_eta_bin[0]-1 ; i++){
                c6 -> cd(i+1); gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);
		h1_dpp_v_p_et_bins_1_4_over_3_0[0][i] -> Draw();
                for(int j = 0 ; j < size_loaded/2 ; j++){
			if(plotme[j])
				h1_dpp_v_p_et_bins_1_4_over_3_0[j][i] -> Draw("same");
		}
	}
	c6 -> cd(6);
	leg3 -> Draw("same");
	c6 -> Modified();
        c6 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	TString out_pdf_name = "results_code_newer.pdf";
	c0 -> Print(out_pdf_name+"(");
	c1 -> Print(out_pdf_name    );
	c2 -> Print(out_pdf_name    );
	c3 -> Print(out_pdf_name    );
	c4 -> Print(out_pdf_name    );
	c5 -> Print(out_pdf_name    );
	c6 -> Print(out_pdf_name+")");

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
	h1->GetYaxis()->SetTitleOffset(1.2);

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
	g_l1->SetFillColorAlpha(color,0.2);
        g_l1->SetMarkerStyle(marker);
        g_l1->SetMarkerSize(1.4);
        g_l1->SetLineWidth(2);

        g_l1->GetXaxis()->SetTitle(xax);
        g_l1->GetXaxis()->SetNdivisions(108);
        g_l1->GetXaxis()->SetTitleSize(0.07);
        g_l1->GetXaxis()->SetLabelSize(0.07);
        g_l1->GetXaxis()->CenterTitle();

        g_l1->GetYaxis()->SetTitle(yax);
        g_l1->GetYaxis()->SetNdivisions(108);
        g_l1->GetYaxis()->SetTitleSize(0.07);
        g_l1->GetYaxis()->SetLabelSize(0.07);
        g_l1->GetYaxis()->CenterTitle();
	g_l1->GetYaxis()->SetTitleOffset(1.);

        g_l1->SetTitle(tit);

        if(min!=999) g_l1 -> SetMinimum(min);
        if(max!=999) g_l1 -> SetMaximum(max);

        float minxval = h1 -> GetBinLowEdge(1);
        float maxxval = h1 -> GetBinLowEdge(nbin) + h1 -> GetBinWidth(nbin);

        g_l1->GetXaxis()->SetRangeUser( minxval , maxxval );

        return g_l1;
}

