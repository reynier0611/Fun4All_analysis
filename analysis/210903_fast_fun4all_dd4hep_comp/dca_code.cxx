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

double mom_DD4HEP  [7]    =  {0};
double eta_DD4HEP  [7]    =  {0};
double dpp_DD4HEP  [7][7] = {{0}};
double E_dpp_DD4HEP[7][7] = {{0}};
double dth_DD4HEP  [7][7] = {{0}};
double E_dth_DD4HEP[7][7] = {{0}};
double dph_DD4HEP  [7][7] = {{0}};
double E_dph_DD4HEP[7][7] = {{0}};

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max );
void load_DD4HEP_Shujie( TString filename );
double mom_theta_from_string( std::string st , std::string subst );
double max_val_hist(TH1F * h);
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetErrorX(0.0001);
	int color[] = {64,92,8,98,51,50,6};
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		//"output_vtx_res_combined_fun4all_to_comp_to_dd4hep_first_eta_0_005_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507sigma_eta_7_p_8_.root",
		//"output_vtx_res_combined_fun4all_to_comp_to_dd4hep_first_eta_0_005_including_beampipe_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507sigma_eta_7_p_8_.root"
		"output_vtx_res_combined_fun4all_to_comp_to_dd4hep_first_eta_0_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507sigma_eta_7_p_8_.root",
                "output_vtx_res_combined_fun4all_to_comp_to_dd4hep_first_eta_0_0001_including_beampipe_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507sigma_eta_7_p_8_.root"
	};

	TString label[] = {
		"B = 3.0 T (ATHENA 05/07/21) (all-si)"
	};

	TString p_vals[] = {"1","2","5","10","15","20"};
	TString eta_vals[] = {"0","0.509","1.011","1.506","2.028","2.542","3.131"};

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

	TH1F *** h1_dvl_v_p_et_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dvt_v_p_et_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dvl_v_et_p_bins   = new TH1F ** [size_loaded];
	TH1F *** h1_dvt_v_et_p_bins   = new TH1F ** [size_loaded];

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

		h1_dvl_v_p_et_bins  [f] = new TH1F * [num_eta_bin[f]];
		h1_dvt_v_p_et_bins  [f] = new TH1F * [num_eta_bin[f]];
		h1_dvl_v_et_p_bins  [f] = new TH1F * [num_mom_bin[f]];
		h1_dvt_v_et_p_bins  [f] = new TH1F * [num_mom_bin[f]];

		for(int et = 0 ; et < num_eta_bin[f] ; et++){
			h1_dvl_v_p_et_bins  [f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dvl_v_p_et_bins_%i"  ,et));
			h1_dvt_v_p_et_bins  [f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dvt_v_p_et_bins_%i"  ,et));		
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			h1_dvl_v_et_p_bins  [f][p] = (TH1F*) Fin[f] -> Get(Form("h1_dvl_v_et_p_bins_%i"  ,p));
			h1_dvt_v_et_p_bins  [f][p] = (TH1F*) Fin[f] -> Get(Form("h1_dvt_v_et_p_bins_%i"  ,p));
		}
	}

	cout << "\neta bin boundaries:\n"; for(int et = 0 ; et < num_eta_bin[0]+1 ; et++) cout << (*TVT_eta_bin[0])[et] << ", "; cout << "\n";
	cout << "\np bin boundaries:\n"  ; for(int p  = 0 ; p  < num_mom_bin[0]+1 ; p ++) cout << (*TVT_mom_bin[0])[ p] << ", "; cout << "\n\n";

	// #######################################################################################################################################
	// EDIT THE CODE BELOW DEPENDING ON WHAT YOU WANT TO PLOT
	// ------------------------------------------------------------------------------
	TGraphErrors *** g_dvl_v_p_et_bins   = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dvt_v_p_et_bins   = new TGraphErrors ** [size_loaded];

	TGraphErrors *** g_dvl_v_et_p_bins   = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dvt_v_et_p_bins   = new TGraphErrors ** [size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		std::string temp_lab = (std::string)(label[f]);
		double max_dpp = temp_lab.find("1.4") != std::string::npos ? 13.:6.;

		g_dvl_v_p_et_bins  [f] = new TGraphErrors * [num_eta_bin[f]];
		g_dvt_v_p_et_bins  [f] = new TGraphErrors * [num_eta_bin[f]];

		g_dvl_v_et_p_bins  [f] = new TGraphErrors * [num_mom_bin[f]];
		g_dvt_v_et_p_bins  [f] = new TGraphErrors * [num_mom_bin[f]];

		for(int et = 0 ; et < num_eta_bin[f] ; et++){
			g_dvl_v_p_et_bins  [f][et] = graph_from_histo( h1_dvl_v_p_et_bins[f][et] , f==0?62:8 , f==0?22:20 , 0.0 , max_dpp );
			g_dvt_v_p_et_bins  [f][et] = graph_from_histo( h1_dvt_v_p_et_bins[f][et] , f==0?62:8 , f==0?22:20 , 2e-2, 100 );
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			g_dvl_v_et_p_bins  [f][p] = graph_from_histo( h1_dvl_v_et_p_bins[f][p] , 51 + p*5 , 20 , 0.0 , max_dpp );
			g_dvt_v_et_p_bins  [f][p] = graph_from_histo( h1_dvt_v_et_p_bins[f][p] , 51 + p*5 , 20 , 2e-2, 100 );
		}
	}
	// ------------------------------------------------------------------------------
	// Results from fast simulations (Ernst Sichtermann)
	TFile * F_fast_wo_beampipe = new TFile("Ernst_fast_sim/fast_sim_performance_wo_beampipe.root");	
	TFile * F_fast_w_beampipe  = new TFile("Ernst_fast_sim/fast_sim_performance_w_beampipe.root" );
	
	TGraph ** g_fast_DCAt_v_mom_wo_beampipe = new TGraph * [7];
	TGraph ** g_fast_DCAz_v_mom_wo_beampipe = new TGraph * [7];

	g_fast_DCAt_v_mom_wo_beampipe[0] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAt_v_mom_eta_0_000");
	g_fast_DCAt_v_mom_wo_beampipe[1] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAt_v_mom_eta_0_509");
	g_fast_DCAt_v_mom_wo_beampipe[2] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAt_v_mom_eta_1_011");
	g_fast_DCAt_v_mom_wo_beampipe[3] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAt_v_mom_eta_1_506");
	g_fast_DCAt_v_mom_wo_beampipe[4] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAt_v_mom_eta_2_028");
	g_fast_DCAt_v_mom_wo_beampipe[5] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAt_v_mom_eta_2_542");
	g_fast_DCAt_v_mom_wo_beampipe[6] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAt_v_mom_eta_3_131");

	g_fast_DCAz_v_mom_wo_beampipe[0] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAz_v_mom_eta_0_000");
        g_fast_DCAz_v_mom_wo_beampipe[1] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAz_v_mom_eta_0_509");
        g_fast_DCAz_v_mom_wo_beampipe[2] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAz_v_mom_eta_1_011");
        g_fast_DCAz_v_mom_wo_beampipe[3] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAz_v_mom_eta_1_506");
        g_fast_DCAz_v_mom_wo_beampipe[4] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAz_v_mom_eta_2_028");
        g_fast_DCAz_v_mom_wo_beampipe[5] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAz_v_mom_eta_2_542");
        g_fast_DCAz_v_mom_wo_beampipe[6] = (TGraph*) F_fast_wo_beampipe -> Get("g_DCAz_v_mom_eta_3_131");	

	for(int i = 0 ; i < 7 ; i++){
		g_fast_DCAt_v_mom_wo_beampipe[i] -> SetMarkerStyle(22);
		g_fast_DCAz_v_mom_wo_beampipe[i] -> SetMarkerStyle(22);
	}

	
        TGraph ** g_fast_DCAt_v_mom_w_beampipe = new TGraph * [7];
        TGraph ** g_fast_DCAz_v_mom_w_beampipe = new TGraph * [7];

        g_fast_DCAt_v_mom_w_beampipe[0] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAt_v_mom_eta_0_000");
        g_fast_DCAt_v_mom_w_beampipe[1] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAt_v_mom_eta_0_509");
        g_fast_DCAt_v_mom_w_beampipe[2] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAt_v_mom_eta_1_011");
        g_fast_DCAt_v_mom_w_beampipe[3] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAt_v_mom_eta_1_506");
        g_fast_DCAt_v_mom_w_beampipe[4] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAt_v_mom_eta_2_028");
        g_fast_DCAt_v_mom_w_beampipe[5] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAt_v_mom_eta_2_542");
        g_fast_DCAt_v_mom_w_beampipe[6] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAt_v_mom_eta_3_131");

        g_fast_DCAz_v_mom_w_beampipe[0] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAz_v_mom_eta_0_000");
        g_fast_DCAz_v_mom_w_beampipe[1] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAz_v_mom_eta_0_509");
        g_fast_DCAz_v_mom_w_beampipe[2] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAz_v_mom_eta_1_011");
        g_fast_DCAz_v_mom_w_beampipe[3] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAz_v_mom_eta_1_506");
        g_fast_DCAz_v_mom_w_beampipe[4] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAz_v_mom_eta_2_028");
        g_fast_DCAz_v_mom_w_beampipe[5] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAz_v_mom_eta_2_542");
        g_fast_DCAz_v_mom_w_beampipe[6] = (TGraph*) F_fast_w_beampipe -> Get("g_DCAz_v_mom_eta_3_131");

	// ------------------------------------------------------------------------------
	// Plotting graphs
	TLegend * leg = new TLegend(0,0.2,1,0.8);
	leg -> SetLineColor(0);
	//leg -> AddEntry( g_fast_dvl_v_p[0]       ,"fast"   );
	//leg -> AddEntry((TObject*)0, "" , "");
	leg -> AddEntry((TObject*)0, "B = ATHENA 05/07", "");
	leg -> AddEntry( g_dvl_v_p_et_bins[0][0] ,"Fun4All (w/o beampipe)");
	leg -> AddEntry( g_dvl_v_p_et_bins[1][0] ,"Fun4All (w/ beampipe)");
	//leg -> AddEntry( g_dd4hep_dvl_v_p[0]     ,"DD4HEP" );
	// ------------------------------------------------------------------------------
        TCanvas * c1 = new TCanvas("c1","c1",1300,900);
	c1 -> Divide(4,2);	
	for(int et = 0 ; et < num_eta_bin[0] ; et++){
        	c1 -> cd(et+1); gPad -> SetLeftMargin(0.25); gPad -> SetRightMargin(0.05); gPad -> SetBottomMargin(0.17);
		g_dvl_v_p_et_bins[1][et] -> SetMaximum(2.2*max_val_hist(h1_dvl_v_p_et_bins[0][et]));
		g_dvl_v_p_et_bins[1][et] -> SetMinimum(0);
		g_dvl_v_p_et_bins[1][et] -> SetTitle("#eta = "+eta_vals[et]);
		g_dvl_v_p_et_bins[1][et] -> Draw("APL");
        	for(int f = 0 ; f < size_loaded ; f++){
			g_dvl_v_p_et_bins[f][et] -> Draw("samePL");
			g_fast_DCAz_v_mom_wo_beampipe[et] -> Draw("sameP");
			g_fast_DCAz_v_mom_w_beampipe [et] -> Draw("sameP");
		}
	}
	c1 -> cd(8);
	leg -> Draw();
	// ------------
        c1 -> Modified();
        c1 -> Update();
	// ------------------------------------------------------------------------------
        TCanvas * c2 = new TCanvas("c2","c2",1300,900);
        c2 -> Divide(4,2);
        for(int et = 0 ; et < num_eta_bin[0] ; et++){
                c2 -> cd(et+1); gPad -> SetLeftMargin(0.25); gPad -> SetRightMargin(0.05); gPad -> SetBottomMargin(0.17);
                g_dvt_v_p_et_bins[1][et] -> SetMaximum(2.2*max_val_hist(h1_dvt_v_p_et_bins[0][et]));
                g_dvt_v_p_et_bins[1][et] -> SetMinimum(0);
                g_dvt_v_p_et_bins[1][et] -> SetTitle("#eta = "+eta_vals[et]);
                g_dvt_v_p_et_bins[1][et] -> Draw("APL");
		for(int f = 0 ; f < size_loaded ; f++){
			g_dvt_v_p_et_bins[f][et] -> Draw("samePL");
        		g_fast_DCAt_v_mom_wo_beampipe[et] -> Draw("sameP");
                        g_fast_DCAt_v_mom_w_beampipe [et] -> Draw("sameP");
		}
	}
        c2 -> cd(8);
        leg -> Draw();
        // ------------
        c2 -> Modified();
        c2 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_dca_z.pdf");
	c2 -> Print("results_dca_t.pdf");

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

			if     (xval[ctr]<1.5) xval[ctr] = 1;
			else if(xval[ctr]<3.0) xval[ctr] = 2;
			else if(xval[ctr]<7.0) xval[ctr] = 5;
			else if(xval[ctr]<12.) xval[ctr] = 10;
			else if(xval[ctr]<17.) xval[ctr] = 15;
			else if(xval[ctr]<25.) xval[ctr] = 20;
			else if(xval[ctr]<40.) xval[ctr] = 30;
			else if(xval[ctr]<70.) xval[ctr] = 50;

			err[ctr] = h1 -> GetBinError(i+1);

			cout << xval[ctr] << "\t" << yval[ctr] << "\t" << err[ctr] << endl;
			ctr++;
		}
	}

	// Creating a TGraphErrors object that mirrors the histogram
	TGraphErrors * g_l1 = new TGraphErrors(ctr,xval,yval,0,err);
	g_l1->SetLineColor(color);
	g_l1->SetMarkerColor(color);
	g_l1->SetMarkerStyle(marker);
	g_l1->SetMarkerSize(1.5);
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

	g_l1->SetTitle(tit);

	if(min!=999) g_l1 -> SetMinimum(min);
	if(max!=999) g_l1 -> SetMaximum(max);

	float minxval = h1 -> GetBinLowEdge(1);
	float maxxval = h1 -> GetBinLowEdge(nbin) + h1 -> GetBinWidth(nbin); 

	g_l1->GetXaxis()->SetRangeUser( minxval , maxxval );

	return g_l1;
}
// ==========================================================================
void load_DD4HEP_Shujie( TString filename ){
        ifstream f;
        std::string descriptor;
        double val1, val2, val3, val4, val5, val6, val7;
        double mom, theta, eta;
        f.open(filename);
        int ctr = 0;
        int ctr2 = 0;
        int et_idx = 0;
        while(!f.eof()){
                f >> descriptor; // String with some information
                f >> val1; f >> val2; f >> val3; f >> val4; f >> val5; f >> val6; f >> val7; // numerical values
                // --------------------------------
                // Get values from descriptor string
                mom = mom_theta_from_string  ( descriptor , "_p" ); // Get momentum value from descriptor
                theta = mom_theta_from_string( descriptor , "_th"); // Get polar-angle value from descriptor
                theta *= TMath::Pi()/180.; // Convert theta to radians
                eta = -TMath::Log(TMath::Tan(theta/2.)); // determine pseudorapitidy from polar angle
                // --------------------------------
                bool eta_already_there = false;
                for(int i=0;i<sizeof(eta_DD4HEP)/sizeof(*eta_DD4HEP);i++){
                        if(eta==eta_DD4HEP[i]) eta_already_there = true;
                }
                if(!eta_already_there)
                        eta_DD4HEP[ctr] = eta;
                // --------------------------------
                if(ctr2==0||mom!=mom_DD4HEP[ctr2-1]){
                        mom_DD4HEP[ctr2] = mom;
                        ctr2++;
                        et_idx = 0;
                }
                // --------------------------------
                dpp_DD4HEP  [et_idx][ctr2-1] = val2;
                E_dpp_DD4HEP[et_idx][ctr2-1] = val3;
                dth_DD4HEP  [et_idx][ctr2-1] = val4;
                E_dth_DD4HEP[et_idx][ctr2-1] = val5;
                dph_DD4HEP  [et_idx][ctr2-1] = val6;
                E_dph_DD4HEP[et_idx][ctr2-1] = val7;

                et_idx++;
                ctr++;
        }

        f.close();
}
// ==========================================================================
double mom_theta_from_string( std::string st , std::string subst ){
        int step = 2;
        if(subst=="_th") step = 3;
        //int idx = st.find(subst); // position of characters 'subst' in string 'st'
        std::string new_str = st.substr(st.find(subst)+step,2); // copy substring right after 'subst'
        if(new_str.substr(1)=="_") new_str = new_str.substr(0,1); // if last character is '_', remove it
        double result = stod(new_str);
        return result;
}
// ============================================================================================================================================
double max_val_hist(TH1F * h){
        int nBins = h -> GetNbinsX();
        double max_val = -99999.;
        for(int i = 0 ; i<nBins ; i++){
                double bc = h -> GetBinContent(i+1) + h -> GetBinError(i+1);
                if(bc>max_val) max_val = bc;
        }
        return max_val;
}

