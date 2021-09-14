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
		"output_vtx_res_combined_fun4all_to_comp_to_dd4hep_first_eta_0_005_AllSi_vbd_0.05_0.55_0.24_B_ATHENA_210507sigma_eta_7_p_8_.root"
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
			g_dvl_v_p_et_bins  [f][et] = graph_from_histo( h1_dvl_v_p_et_bins  [f][et] , 62 , 20 , 0.0 , max_dpp );
			g_dvt_v_p_et_bins  [f][et] = graph_from_histo( h1_dvt_v_p_et_bins  [f][et] , 62 , 20 , 2e-2, 100 );
		}

		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			g_dvl_v_et_p_bins  [f][p] = graph_from_histo( h1_dvl_v_et_p_bins  [f][p] , 51 + p*5 , 20 , 0.0 , max_dpp );
			g_dvt_v_et_p_bins  [f][p] = graph_from_histo( h1_dvt_v_et_p_bins  [f][p] , 51 + p*5 , 20 , 2e-2, 100 );
		}
	}
	// ------------------------------------------------------------------------------
	// Results from fast simulations (Ernst Sichtermann)
	/*
	double fast_p[] = {1,5,10,15,20,30,40,50};
	double fast_dpp[5][8] = {
		{0.3049,0.5241,0.6117,0.6618,0.7072,0.8027,0.9108,1.028 },
		{0.3112,0.5319,0.6290,0.6822,0.7229,0.8027,0.8904,0.9890},
		{0.3519,0.5413,0.6697,0.7339,0.7793,0.8434,0.9014,0.9609},
		{0.4082,0.4207,0.4380,0.4646,0.4975,0.5804,0.6791,0.7871},
		{0.5068,0.4755,0.4912,0.5022,0.5162,0.5491,0.5898,0.6368}
	};
	TGraph ** g_fast_dvl_v_p = new TGraph*[5];
	for(int i = 0 ; i < 5 ; i++){
		g_fast_dvl_v_p[i] = new TGraph(8,fast_p,fast_dpp[i]);
		g_fast_dvl_v_p[i] -> SetLineColor(1);
		g_fast_dvl_v_p[i] -> SetLineStyle(2);
		g_fast_dvl_v_p[i] -> SetLineWidth(3);
	}
	*/
	// ------------------------------------------------------------------------------
        // Results from dd4hep (Shujie Li)
	/*
	double dd4hep_p[] = {1,2,5,10,20,30,50};
	double dd4hep_dpp[5][7] = {
		{0.4609,0.5234,0.6429,0.7132,0.7176,0.7812,0.9487},
		{0.4609,0.5212,0.6652,0.6797,0.7388,0.8158,1.0090},
		{0.5569,0.5324,0.6663,0.8683,0.8705,0.8415,1.0100},
		{0.6116,0.5915,0.4520,0.4342,0.5458,0.5513,0.8281},
		{0.9743,0.8683,0.6585,0.6217,0.6138,0.6205,0.6752}
	};
	TGraph ** g_dd4hep_dvl_v_p = new TGraph*[5];
        for(int i = 0 ; i < 5 ; i++){
                g_dd4hep_dvl_v_p[i] = new TGraph(7,dd4hep_p,dd4hep_dpp[i]);
                g_dd4hep_dvl_v_p[i] -> SetLineColor(2);
		g_dd4hep_dvl_v_p[i] -> SetMarkerColor(2);
		g_dd4hep_dvl_v_p[i] -> SetMarkerStyle(21);
        }
	*/
	/*
	load_DD4HEP_Shujie("tracker_only_1.2sigma_211.txt");
	TGraphErrors ** g_dd4hep_dvl_v_p = new TGraphErrors*[7];
        for(int i = 0 ; i < 7 ; i++){
                g_dd4hep_dvl_v_p[i] = new TGraphErrors(7,mom_DD4HEP,dpp_DD4HEP[i],0,E_dpp_DD4HEP[i]);
                g_dd4hep_dvl_v_p[i] -> SetLineColor(2);
                g_dd4hep_dvl_v_p[i] -> SetMarkerColor(2);
                g_dd4hep_dvl_v_p[i] -> SetMarkerStyle(21);
        }
	*/
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TLegend * leg = new TLegend(0,0.2,1,0.8);
	leg -> SetLineColor(0);
	//leg -> AddEntry( g_fast_dvl_v_p[0]       ,"fast"   );
	//leg -> AddEntry((TObject*)0, "" , "");
	leg -> AddEntry((TObject*)0, "B = ATHENA 05/07", "");
	leg -> AddEntry( g_dvl_v_p_et_bins[0][0] ,"Fun4All");
	//leg -> AddEntry( g_dd4hep_dvl_v_p[0]     ,"DD4HEP" );
	// ------------------------------------------------------------------------------
	TCanvas ** c1 = new TCanvas*[size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
                c1[f] = new TCanvas(Form("c1_%i",f),Form("c1_%i",f),1300,900);
		c1[f] -> Divide(4,2);
	
		for(int et = 0 ; et < num_eta_bin[f] ; et++){
                	c1[f] -> cd(et+1); gPad -> SetLeftMargin(0.25); gPad -> SetRightMargin(0.05); gPad -> SetBottomMargin(0.17);
			g_dvl_v_p_et_bins[f][et] -> SetMaximum(1.2*max_val_hist(h1_dvl_v_p_et_bins[f][et]));
			g_dvl_v_p_et_bins[f][et] -> SetMinimum(0);
			g_dvl_v_p_et_bins[f][et] -> SetTitle("#eta = "+eta_vals[et]);
			g_dvl_v_p_et_bins[f][et] -> Draw("APL");
                }
		// ------------
		/*
		for(int i = 0 ; i < 7 ; i++){
                        c1[f] -> cd(i+1);
			if(i<5) g_fast_dvl_v_p[i] -> Draw("same");
                        g_dd4hep_dvl_v_p[i] -> Draw("samePL");
                }
		*/
		c1[f] -> cd(8);
		leg -> Draw();
		// ------------
                c1[f] -> Modified();
                c1[f] -> Update();
	}

	// ------------------------------------------------------------------------------
        TCanvas ** c2 = new TCanvas*[size_loaded];

        for(int f = 0 ; f < size_loaded ; f++){
                c2[f] = new TCanvas(Form("c2_%i",f),Form("c2_%i",f),1300,900);
                c2[f] -> Divide(4,2);

                for(int et = 0 ; et < num_eta_bin[f] ; et++){
                        c2[f] -> cd(et+1); gPad -> SetLeftMargin(0.25); gPad -> SetRightMargin(0.05); gPad -> SetBottomMargin(0.17);
                        g_dvt_v_p_et_bins[f][et] -> SetMaximum(1.2*max_val_hist(h1_dvt_v_p_et_bins[f][et]));
                        g_dvt_v_p_et_bins[f][et] -> SetMinimum(0);
                        g_dvt_v_p_et_bins[f][et] -> SetTitle("#eta = "+eta_vals[et]);
                        g_dvt_v_p_et_bins[f][et] -> Draw("APL");
                }
                // ------------
                /*
                for(int i = 0 ; i < 7 ; i++){
                        c2[f] -> cd(i+1);
                        if(i<5) g_fast_dvl_v_p[i] -> Draw("same");
                        g_dd4hep_dvl_v_p[i] -> Draw("samePL");
                }
                */
                c2[f] -> cd(8);
                leg -> Draw();
                // ------------
                c2[f] -> Modified();
                c2[f] -> Update();
        }
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	if(size_loaded==1){
		c1[0] -> Print("results_dca_z.pdf");
		c2[0] -> Print("results_dca_t.pdf");
	}
	else{
		cout << "IMPLEMENT THIS! Bailing out!" << endl;
		exit(0);
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

			//cout << xval[ctr] << "\t" << yval[ctr] << "\t" << err[ctr] << endl;
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

