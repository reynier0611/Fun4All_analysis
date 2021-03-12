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
		"../../output/output_mom_res_skimmed_pi-_det2_20x20_Beast_FastTrackingEvalsigma_eta_16_p_10_.root"                 ,	//  0
		"../../output/output_mom_res_skimmed_pi-_det2_10x10_Beast_FastTrackingEvalsigma_eta_16_p_10_.root"                 ,	//  1
		"../../output/output_mom_res_skimmed_pi-_det2_20x20_sPHENIX_FastTrackingEvalsigma_eta_16_p_10_.root"               ,	//  2
		"../../output/output_mom_res_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEvalsigma_eta_16_p_10_.root"               ,	//  3
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
	double max[] = {29,13,7,3,2,2,2,2,2,2,2,2,3,7,13,29};	
	// Editing histograms
	for(int i = 0 ; i < num_eta_bin[0] ; i++){ 
		h1_dpp_v_p_et_bins[0][i] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dpp_v_p_et_bins[1][i] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dpp_v_p_et_bins[2][i] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dpp_v_p_et_bins[3][i] -> SetTitle(Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1])); 

		prettyTH1( h1_dpp_v_p_et_bins[0][i] ,  2 , 21 , 999 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[1][i] , 96 , 20 , 999 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[2][i] ,  4 , 22 , 999 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[3][i] , 62 , 24 , 999 , max[i] );
        
		h1_dppT_v_pT_et_bins[0][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dppT_v_pT_et_bins[1][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dppT_v_pT_et_bins[2][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dppT_v_pT_et_bins[3][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));

                prettyTH1( h1_dppT_v_pT_et_bins[0][i] ,  2 , 21 , 999 , max[i] );
                prettyTH1( h1_dppT_v_pT_et_bins[1][i] , 96 , 20 , 999 , max[i] );
                prettyTH1( h1_dppT_v_pT_et_bins[2][i] ,  4 , 22 , 999 , max[i] );
                prettyTH1( h1_dppT_v_pT_et_bins[3][i] , 62 , 24 , 999 , max[i] );
	}
	// ------------------------------------------------------------------------------
	// Doing fits
	TF1 ** f_dpp_20um_30T  = new TF1*[num_eta_bin[0]];
	TF1 ** f_dpp_10um_30T  = new TF1*[num_eta_bin[0]];
	TF1 ** f_dpp_20um_14T  = new TF1*[num_eta_bin[0]];
	TF1 ** f_dpp_10um_14T  = new TF1*[num_eta_bin[0]];
	TF1 ** f_dppT_20um_30T = new TF1*[num_eta_bin[0]];
        TF1 ** f_dppT_10um_30T = new TF1*[num_eta_bin[0]];
        TF1 ** f_dppT_20um_14T = new TF1*[num_eta_bin[0]];
        TF1 ** f_dppT_10um_14T = new TF1*[num_eta_bin[0]];

	for(int i = 0 ; i < num_eta_bin[0] ; i++){
		f_dpp_20um_30T[i] = new TF1(Form("f_dpp_20um_30T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
		f_dpp_10um_30T[i] = new TF1(Form("f_dpp_10um_30T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
		f_dpp_20um_14T[i] = new TF1(Form("f_dpp_20um_14T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
		f_dpp_10um_14T[i] = new TF1(Form("f_dpp_10um_14T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);

		f_dpp_20um_30T[i] -> SetLineColor(2); f_dpp_20um_30T[i] -> SetLineStyle(3);
                f_dpp_10um_30T[i] -> SetLineColor(2); f_dpp_10um_30T[i] -> SetLineStyle(3);
                f_dpp_20um_14T[i] -> SetLineColor(4); f_dpp_20um_14T[i] -> SetLineStyle(3);
                f_dpp_10um_14T[i] -> SetLineColor(4); f_dpp_10um_14T[i] -> SetLineStyle(3);

		h1_dpp_v_p_et_bins[0][i] -> Fit(Form("f_dpp_20um_30T_%i",i),"R");
		h1_dpp_v_p_et_bins[1][i] -> Fit(Form("f_dpp_10um_30T_%i",i),"R");
		h1_dpp_v_p_et_bins[2][i] -> Fit(Form("f_dpp_20um_14T_%i",i),"R");
		h1_dpp_v_p_et_bins[3][i] -> Fit(Form("f_dpp_10um_14T_%i",i),"R");
		
		f_dppT_20um_30T[i] = new TF1(Form("f_dppT_20um_30T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
                f_dppT_10um_30T[i] = new TF1(Form("f_dppT_10um_30T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
                f_dppT_20um_14T[i] = new TF1(Form("f_dppT_20um_14T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);
                f_dppT_10um_14T[i] = new TF1(Form("f_dppT_10um_14T_%i",i),"sqrt(sq([0]*x)+sq([1]))",0,30);

                f_dppT_20um_30T[i] -> SetLineColor(2); f_dppT_20um_30T[i] -> SetLineStyle(3);
                f_dppT_10um_30T[i] -> SetLineColor(2); f_dppT_10um_30T[i] -> SetLineStyle(3);
                f_dppT_20um_14T[i] -> SetLineColor(4); f_dppT_20um_14T[i] -> SetLineStyle(3);
                f_dppT_10um_14T[i] -> SetLineColor(4); f_dppT_10um_14T[i] -> SetLineStyle(3);

                h1_dppT_v_pT_et_bins[0][i] -> Fit(Form("f_dppT_20um_30T_%i",i),"R");
                h1_dppT_v_pT_et_bins[1][i] -> Fit(Form("f_dppT_10um_30T_%i",i),"R");
                h1_dppT_v_pT_et_bins[2][i] -> Fit(Form("f_dppT_20um_14T_%i",i),"R");
                h1_dppT_v_pT_et_bins[3][i] -> Fit(Form("f_dppT_10um_14T_%i",i),"R");
	}
	// ------------------------------------------------------------------------------
	// Legend
	double tex_y_3_0[] = {1.81,1.81,1.81,1.81,2.69,6.27,11.55,25.86};
	double tex_y_1_4[] = {1.57,1.57,1.57,1.57,2.31,5.42,10.00,22.48};
	
	TLegend ** leg = new TLegend*[num_eta_bin[0]];
	TLatex ** tex_3_0_T_par = new TLatex * [num_eta_bin[0]];
	TLatex ** tex_1_4_T_par = new TLatex * [num_eta_bin[0]];

	for(int i = 8 ; i < num_eta_bin[0] ; i++){
		//leg[i] = new TLegend(0.20,0.7,0.96,0.89);
		leg[i] = new TLegend(0.20,0.71,0.5,0.89);
		leg[i] -> SetLineColor(0);
		//leg[i] -> SetHeader("dp/p = Ap #oplus B","C");
		leg[i] -> AddEntry(h1_dpp_v_p_et_bins[1][i]," ");//Form("3.0T, A = %.3f, B = %.3f",f_dpp_10um_30T[i]->GetParameter(0),f_dpp_10um_30T[i]->GetParameter(1)));
		leg[i] -> AddEntry(h1_dpp_v_p_et_bins[3][i]," ");//Form("1.4T, A = %.3f, B = %.3f",f_dpp_10um_14T[i]->GetParameter(0),f_dpp_10um_14T[i]->GetParameter(1)));

		tex_3_0_T_par[i] = new TLatex(3.5,tex_y_3_0[i-8],Form("#bf{3.0T, A = %.2f, B = %.1f}",f_dpp_10um_30T[i]->GetParameter(0),f_dpp_10um_30T[i]->GetParameter(1)));
		tex_1_4_T_par[i] = new TLatex(3.5,tex_y_1_4[i-8],Form("#bf{1.4T, A = %.2f, B = %.1f}",f_dpp_10um_14T[i]->GetParameter(0),f_dpp_10um_14T[i]->GetParameter(1)));

		tex_3_0_T_par[i] -> SetTextSize(0.06);
		tex_1_4_T_par[i] -> SetTextSize(0.06);

		cout << Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]) << endl;
		cout << "3.0T\t" << Form("\tsqrt(pow(%.3f,2)*x*x+pow(%.3f,2))",f_dpp_10um_30T[i]->GetParameter(0),f_dpp_10um_30T[i]->GetParameter(1)) << endl;
		cout << "1.4T\t" << Form("\tsqrt(pow(%.3f,2)*x*x+pow(%.3f,2))",f_dpp_10um_14T[i]->GetParameter(0),f_dpp_10um_14T[i]->GetParameter(1)) << endl;
	}	
	// ------------------------------------------------------------------------------
	// Pad limits
	double lowx[] = {0.00,0.25,0.50,0.75,0.00,0.25,0.50,0.75};
        double lowy[] = {0.50,0.50,0.50,0.50,0.00,0.00,0.00,0.00};
        double higx[] = {0.25,0.50,0.75,1.00,0.25,0.50,0.75,1.00};
        double higy[] = {1.00,1.00,1.00,1.00,0.50,0.50,0.50,0.50};
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1400,900);
	c1 -> Divide(4,2);
	TVirtualPad ** pad1 = new TVirtualPad * [num_eta_bin[0]];
	for(int i = 8 ; i < num_eta_bin[0] ; i++){
		//c1 -> cd(i-7);
		pad1[i-8] = c1 -> cd(i-7);
                pad1[i-8] -> SetPad(Form("pad1_%i",i-8),Form("pad1_%i",i-8),lowx[i-8],lowy[i-8],higx[i-8],higy[i-8],kWhite, 0, 0);
		gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);
		h1_dpp_v_p_et_bins[1][i] -> GetYaxis() -> SetDecimals(1);
		h1_dpp_v_p_et_bins[1][i] -> SetMinimum(0);
		h1_dpp_v_p_et_bins[1][i] -> Draw(      ); // 10um pixel, 3.0T
		h1_dpp_v_p_et_bins[3][i] -> Draw("same"); // 10um pixel, 1.4T
		leg[i] -> Draw("same");
		tex_3_0_T_par[i] -> Draw("same");
		tex_1_4_T_par[i] -> Draw("same");
	}
	// ------------------------------------------------------------------------------
        // Physics Working Group Requirements
        // https://docs.google.com/spreadsheets/d/1ynU7Cu7NlwRvMtbtFdlp_B5xXkw8yBAtWJbenMf-P3U/edit#gid=368031287
        // https://wiki.bnl.gov/eicug/index.php/Yellow_Report_Physics_Common
        TF1 * f_PWG_req_m3_5_m3_0 = new TF1("f_PWG_req_m3_5_m3_0","sqrt(pow(0.10*x,2)+pow(0.5,2))",0,30); f_PWG_req_m3_5_m3_0 -> SetLineColor(11); f_PWG_req_m3_5_m3_0 -> SetLineStyle(2);	c1 -> cd(7); f_PWG_req_m3_5_m3_0 -> Draw("same");
        TF1 * f_PWG_req_m3_0_m2_5 = new TF1("f_PWG_req_m3_0_m2_5","sqrt(pow(0.10*x,2)+pow(0.5,2))",0,30); f_PWG_req_m3_0_m2_5 -> SetLineColor(11); f_PWG_req_m3_0_m2_5 -> SetLineStyle(2);	c1 -> cd(6); f_PWG_req_m3_0_m2_5 -> Draw("same");
        TF1 * f_PWG_req_m2_5_m2_0 = new TF1("f_PWG_req_m2_5_m2_0","sqrt(pow(0.10*x,2)+pow(0.5,2))",0,21); f_PWG_req_m2_5_m2_0 -> SetLineColor(11); f_PWG_req_m2_5_m2_0 -> SetLineStyle(2);	c1 -> cd(5); f_PWG_req_m2_5_m2_0 -> Draw("same");
        TF1 * f_PWG_req_m2_0_m1_5 = new TF1("f_PWG_req_m2_0_m1_5","sqrt(pow(0.05*x,2)+pow(0.5,2))",0,30); f_PWG_req_m2_0_m1_5 -> SetLineColor(11); f_PWG_req_m2_0_m1_5 -> SetLineStyle(2);	c1 -> cd(4); f_PWG_req_m2_0_m1_5 -> Draw("same");
        TF1 * f_PWG_req_m1_5_m1_0 = new TF1("f_PWG_req_m1_5_m1_0","sqrt(pow(0.05*x,2)+pow(0.5,2))",0,30); f_PWG_req_m1_5_m1_0 -> SetLineColor(11); f_PWG_req_m1_5_m1_0 -> SetLineStyle(2);	c1 -> cd(3); f_PWG_req_m1_5_m1_0 -> Draw("same");
        TF1 * f_PWG_req_m1_0_m0_5 = new TF1("f_PWG_req_m1_0_m0_5","sqrt(pow(0.05*x,2)+pow(0.5,2))",0,30); f_PWG_req_m1_0_m0_5 -> SetLineColor(11); f_PWG_req_m1_0_m0_5 -> SetLineStyle(2);	c1 -> cd(2); f_PWG_req_m1_0_m0_5 -> Draw("same");
        TF1 * f_PWG_req_m0_5_p0_0 = new TF1("f_PWG_req_m0_5_p0_0","sqrt(pow(0.05*x,2)+pow(0.5,2))",0,30); f_PWG_req_m0_5_p0_0 -> SetLineColor(11); f_PWG_req_m0_5_p0_0 -> SetLineStyle(2);	c1 -> cd(1); f_PWG_req_m0_5_p0_0 -> Draw("same");
        TF1 * f_PWG_req_p0_0_p0_5 = new TF1("f_PWG_req_p0_0_p0_5","sqrt(pow(0.05*x,2)+pow(0.5,2))",0,30); f_PWG_req_p0_0_p0_5 -> SetLineColor(11); 						c1 -> cd(1); f_PWG_req_p0_0_p0_5 -> Draw("same");
        TF1 * f_PWG_req_p0_5_p1_0 = new TF1("f_PWG_req_p0_5_p1_0","sqrt(pow(0.05*x,2)+pow(0.5,2))",0,30); f_PWG_req_p0_5_p1_0 -> SetLineColor(11); 						c1 -> cd(2); f_PWG_req_p0_5_p1_0 -> Draw("same");
        TF1 * f_PWG_req_p1_0_p1_5 = new TF1("f_PWG_req_p1_0_p1_5","sqrt(pow(0.05*x,2)+pow(1.0,2))",0,23); f_PWG_req_p1_0_p1_5 -> SetLineColor(11); 						c1 -> cd(3); f_PWG_req_p1_0_p1_5 -> Draw("same");
        TF1 * f_PWG_req_p1_5_p2_0 = new TF1("f_PWG_req_p1_5_p2_0","sqrt(pow(0.05*x,2)+pow(1.0,2))",0,23); f_PWG_req_p1_5_p2_0 -> SetLineColor(11); 						c1 -> cd(4); f_PWG_req_p1_5_p2_0 -> Draw("same");
        TF1 * f_PWG_req_p2_0_p2_5 = new TF1("f_PWG_req_p2_0_p2_5","sqrt(pow(0.05*x,2)+pow(1.0,2))",0,30); f_PWG_req_p2_0_p2_5 -> SetLineColor(11); 						c1 -> cd(5); f_PWG_req_p2_0_p2_5 -> Draw("same");
        TF1 * f_PWG_req_p2_5_p3_0 = new TF1("f_PWG_req_p2_5_p3_0","sqrt(pow(0.10*x,2)+pow(2.0,2))",0,30); f_PWG_req_p2_5_p3_0 -> SetLineColor(11); 						c1 -> cd(6); f_PWG_req_p2_5_p3_0 -> Draw("same");
        TF1 * f_PWG_req_p3_0_p3_5 = new TF1("f_PWG_req_p3_0_p3_5","sqrt(pow(0.10*x,2)+pow(2.0,2))",0,30); f_PWG_req_p3_0_p3_5 -> SetLineColor(11); 						c1 -> cd(7); f_PWG_req_p3_0_p3_5 -> Draw("same");
	f_PWG_req_m3_5_m3_0 -> SetMarkerColor(11);
	f_PWG_req_p3_0_p3_5 -> SetMarkerColor(11);
	// --------------
        TLegend * leg_PWG = new TLegend(0.20,0.57,0.61,0.72);
        leg_PWG -> SetLineColor(0);
        leg_PWG -> AddEntry(f_PWG_req_m3_5_m3_0," ");
	leg_PWG -> AddEntry(f_PWG_req_p3_0_p3_5," ");
	TLatex * tex_PWG_for  = new TLatex(5,1.15,"#bf{Forward}");	tex_PWG_for  -> SetTextSize(0.06);
	TLatex * tex_PWG_back = new TLatex(5,1.37,"#bf{Backward}");	tex_PWG_back -> SetTextSize(0.06);
	c1 -> cd(1);
	leg_PWG -> Draw("same");
	tex_PWG_for  -> Draw("same");
	tex_PWG_back -> Draw("same");
	// ----------------------------------
	c1 -> Modified();
	c1 -> Update();
	// ----------------------------------
	/*
	TCanvas * c2 = new TCanvas("c2","c2",1400,900);
        c2 -> Divide(4,2);
        TVirtualPad ** pad2 = new TVirtualPad * [num_eta_bin[0]];
        for(int i = 8 ; i < num_eta_bin[0] ; i++){
                pad2[i-8] = c2 -> cd(i-7);
                pad2[i-8] -> SetPad(Form("pad2_%i",i-8),Form("pad2_%i",i-8),lowx[i-8],lowy[i-8],higx[i-8],higy[i-8],kWhite, 0, 0);
                gPad -> SetRightMargin(0.044); gPad -> SetLeftMargin(0.19); gPad -> SetBottomMargin(0.17);
                h1_dppT_v_pT_et_bins[1][i] -> Draw(      ); // 10um pixel, 3.0T
                h1_dppT_v_pT_et_bins[3][i] -> Draw("same"); // 10um pixel, 1.4T
                //leg[i] -> Draw("same");
                //tex_3_0_T_par[i] -> Draw("same");
                //tex_1_4_T_par[i] -> Draw("same");
        }
        c2 -> Modified();
        c2 -> Update();
	*/
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
