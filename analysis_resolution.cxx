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

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	if(argc!=4){
		cout << "Run this code as:\n\033[1;32m./analysis_resolution A B filename.root\033[0m\n";
		cout << "where:\nA = 1 -> Widths from table will be used\n  = 2 -> Widths from table \033[1;31mwon't\033[0m be used\n";
		cout << "B = 1 -> Table will be updated\n  = 2 -> Table \033[1;31mwon't\033[0m be updated\n";
		exit(0);
	}

	bool use_widths = true;
	bool update_tab = false;

	if (!fs::is_directory("data"  ) || !fs::exists("data"  )) fs::create_directory("data"  ); // Create directory if it does not exist
	if (!fs::is_directory("tables") || !fs::exists("tables")) fs::create_directory("tables");
	if (!fs::is_directory("output") || !fs::exists("output")) fs::create_directory("output");
	if (!fs::is_directory("fits"  ) || !fs::exists("fits"  )) fs::create_directory("fits"  );

	cout << "\033[1;31m********************************************************************\nUSEFUL INFO:\033[0m\nWill be loading data from file: '" << argv[3] << "' assumed to be in directory 'data'" << endl;

	if     (atoi(argv[1])==1){use_widths = true ;	cout << "Will be using widths from table\n" ;}
	else if(atoi(argv[1])==2){use_widths = false;	cout << "Won't be using widths from table\n";}
	else{cout << "Something wrong with your election of input parameter 'A'. Bailing out!\n"; exit(0);}

	if     (atoi(argv[2])==1){update_tab = true ;   cout << "Table will be updated\n" ;}
	else if(atoi(argv[2])==2){update_tab = false;   cout << "Table won't be updated\n";}
	else{cout << "Something wrong with your election of input parameter 'B'. Bailing out!\n"; exit(0);}

	// -------------------------
	// Binning
	float eta_bin[] = {-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75};
	float mom_bin[] = {0.,1.,2.,3.,4.,5.,10.,15.,20.,25.,30.};

	const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	const int size_mom_bin = sizeof(mom_bin)/sizeof(*mom_bin);

	TVectorT<double> TVT_eta_bin(size_eta_bin);	for(int i = 0 ; i < size_eta_bin ; i++) TVT_eta_bin[i] = eta_bin[i];
	TVectorT<double> TVT_mom_bin(size_mom_bin);	for(int i = 0 ; i < size_mom_bin ; i++) TVT_mom_bin[i] = mom_bin[i];
	// -------------------------
	// useful strings
	string raw_fname = argv[3];
	TString infile = "data/" + raw_fname;
	TString outfile = "output/output_" + raw_fname;
	raw_fname.resize(raw_fname.size()-5);
	TString tab_name = "tables/tab_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".txt";
	TString out_pdf = "output_fits_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
	TString out_pdf2 = "results_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
	// -------------------------------------------------------------
	// Some settings
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);	
	// -------------------------------------------------------------
	// Loading all the needed info from the root file
	TFile * F = new TFile(infile);
	TTree * T = (TTree*) F -> Get("tracks");
	float gpx, gpy, gpz, px, py, pz;
	T -> SetBranchAddress("gpx"  ,&gpx  );
	T -> SetBranchAddress("gpy"  ,&gpy  );
	T -> SetBranchAddress("gpz"  ,&gpz  );
	T -> SetBranchAddress("px"   ,&px   );
	T -> SetBranchAddress("py"   ,&py   );
	T -> SetBranchAddress("pz"   ,&pz   ); 
	int nEntries = T -> GetEntries();
	// -------------------------------------------------------------
	fstream tab;
	float approx_sig_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	TString temp_str;
	if(use_widths){
		tab.open(tab_name);
		if(!tab){cout << "Could not find file '" << tab_name << "'" << endl; use_widths = false; update_tab = true;}
		else{
			cout << "Loading parameters from file '" << tab_name << "'" << endl;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dpp[et][p];}}	//tab >> temp_str;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dth[et][p];}}	//tab >> temp_str;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dph[et][p];}}
		}
		tab.close();
	}

	float approx_sig_dpp_3_0[size_eta_bin-1][size_mom_bin-1] = {{0}}; float approx_sig_dpp_1_1[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dth_3_0[size_eta_bin-1][size_mom_bin-1] = {{0}}; float approx_sig_dth_1_1[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dph_3_0[size_eta_bin-1][size_mom_bin-1] = {{0}}; float approx_sig_dph_1_1[size_eta_bin-1][size_mom_bin-1] = {{0}};

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			approx_sig_dpp_3_0[et][p] = 3.0*approx_sig_dpp[et][p]; approx_sig_dpp_1_1[et][p] = 1.1*approx_sig_dpp[et][p];
			approx_sig_dth_3_0[et][p] = 3.0*approx_sig_dth[et][p]; approx_sig_dth_1_1[et][p] = 1.1*approx_sig_dth[et][p];
			approx_sig_dph_3_0[et][p] = 3.0*approx_sig_dph[et][p]; approx_sig_dph_1_1[et][p] = 1.1*approx_sig_dph[et][p];
		}
	}
	// -------------------------------------------------------------
	// Defining histograms
	TH1F *** h1_dpp_p_et_bins = new TH1F**[size_eta_bin-1];	// delta p / p vs. p in eta bins
	TH1F *** h1_dth_p_et_bins = new TH1F**[size_eta_bin-1];	// delta theta vs. p in eta bins
	TH1F *** h1_dph_p_et_bins = new TH1F**[size_eta_bin-1];	// delta phi   vs. p in eta bins
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dth_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dph_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dp/p;Counts"         ,80,-approx_sig_dpp_3_0[et][p],approx_sig_dpp_3_0[et][p]);
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",80,-approx_sig_dth_3_0[et][p],approx_sig_dth_3_0[et][p]);
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts"  ,80,-approx_sig_dph_3_0[et][p],approx_sig_dph_3_0[et][p]);
			}
			else{
				h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dp/p;Counts"         ,80,-0.08  ,0.08  );
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",80,-0.0014,0.0014);
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts"  ,80,-0.04  ,0.04  );
			}

			h1_dpp_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dth_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dph_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
		}
	}
	// -------------------------------------------------------------	
	//int color[] = {1,2,62,8,95,52,6,28,209,92,15};
	//int marker[] = {20,21,23,24,25,26,27,28,29,30};

	TH1F ** h1_dpp_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dth_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dph_v_p_et_bins = new TH1F*[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_v_p_et_bins[et] = new TH1F(Form("h1_dpp_v_p_et_bins_%i",et),";p [GeV/c];dp/p [%]"      ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dpp_v_p_et_bins[et] , 50+et /*color[et]*/ , 20 /*marker[et]*/ , 0. , 10. );
		h1_dth_v_p_et_bins[et] = new TH1F(Form("h1_dth_v_p_et_bins_%i",et),";p [GeV/c];d#theta [mrad]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dth_v_p_et_bins[et] , 50+et /*color[et]*/ , 20 /*marker[et]*/ , 0. , 1.  );
		h1_dph_v_p_et_bins[et] = new TH1F(Form("h1_dph_v_p_et_bins_%i",et),";p [GeV/c];d#phi [mrad]"  ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dph_v_p_et_bins[et] , 50+et /*color[et]*/ , 20 /*marker[et]*/ , 0. , 25. );
	}

	TH1F ** h1_dpp_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dth_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dph_v_et_p_bins = new TH1F*[size_mom_bin-1];

	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dpp_v_et_p_bins[p] = new TH1F(Form("h1_dpp_v_et_p_bins_%i",p),";#eta;dp/p [%]"      ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dpp_v_et_p_bins[p] , 50+p /*color[p]*/ , 20 /*marker[p]*/ , 0. , 10. );
		h1_dth_v_et_p_bins[p] = new TH1F(Form("h1_dth_v_et_p_bins_%i",p),";#eta;d#theta [mrad]",size_eta_bin-1,eta_bin);	prettyTH1F( h1_dth_v_et_p_bins[p] , 50+p /*color[p]*/ , 20 /*marker[p]*/ , 0. , 1.  );
		h1_dph_v_et_p_bins[p] = new TH1F(Form("h1_dph_v_et_p_bins_%i",p),";#eta;d#phi [mrad]"  ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dph_v_et_p_bins[p] , 50+p /*color[p]*/ , 20 /*marker[p]*/ , 0. , 25. );
	}
	// -------------------------------------------------------------
	// Declaring other useful variables and functions
	float width_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};

	TF1 *** f_gaus_dpp = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dth = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dph = new TF1**[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		f_gaus_dpp[et] = new TF1*[size_mom_bin-1];
		f_gaus_dth[et] = new TF1*[size_mom_bin-1];
		f_gaus_dph[et] = new TF1*[size_mom_bin-1];

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-approx_sig_dpp_1_1[et][p],approx_sig_dpp_1_1[et][p]);
				f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-approx_sig_dth_1_1[et][p],approx_sig_dth_1_1[et][p]);
				f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-approx_sig_dph_1_1[et][p],approx_sig_dph_1_1[et][p]);
			}
			else{
				f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-0.007 ,0.007 );
				f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-0.0007,0.0007);
				f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-0.02  ,0.02  );
			}
		}
	}
	cout << "\033[1;31m********************************************************************\033[0m\n";
	// -------------------------------------------------------------
	// Loop over entries of the tree	
	for(int ev = 0 ; ev < nEntries ; ev++){
		T -> GetEntry(ev);
		if(ev%1000000==0) cout << "Looping over entry " << ev << " out of " << nEntries << endl;

		// Calculating some variables
		float gtheta = TMath::ACos(gpz/sqrt(gpx*gpx+gpy*gpy+gpz*gpz));	
		float theta = TMath::ACos(pz/sqrt(px*px+py*py+pz*pz));
		float dth = theta - gtheta;

		float geta = -TMath::Log(TMath::Tan(gtheta/2.));

		float p_reco = sqrt(px*px+py*py+pz*pz);
		float p_truth = sqrt(gpx*gpx+gpy*gpy+gpz*gpz);
		float dp_p = (p_reco-p_truth)/p_truth;

		float gphi = TMath::ATan(gpy/gpx);
		float phi = TMath::ATan(py/px);
		float dph = phi - gphi;

		// Filling histograms
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			if( geta >  eta_bin[et] &&  geta <= eta_bin[et+1] ){
				//if( abs(geta) >  eta_bin[et] &&  abs(geta) <= eta_bin[et+1] ){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					if( p_truth > mom_bin[p] && p_truth <= mom_bin[p+1] ){
						h1_dpp_p_et_bins[et][p] -> Fill( dp_p );
						h1_dth_p_et_bins[et][p] -> Fill( dth  );
						h1_dph_p_et_bins[et][p] -> Fill( dph  );
					}	
				}
			}
			}
		}
		cout << "\033[1;31m********************************************************************\033[0m\n";
		// -------------------------------------------------------------
		// Doing fits
		TCanvas ** c_fits_p  = new TCanvas*[size_eta_bin-1];
		TCanvas ** c_fits_th = new TCanvas*[size_eta_bin-1];
		TCanvas ** c_fits_ph = new TCanvas*[size_eta_bin-1];

		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			c_fits_p [et] = new TCanvas(Form("c_fits_p_%i" ,et),Form("dp/p  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_p [et] -> Divide(5,2);
			c_fits_th[et] = new TCanvas(Form("c_fits_th_%i",et),Form("dtheta, %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_th[et] -> Divide(5,2);
			c_fits_ph[et] = new TCanvas(Form("c_fits_ph_%i",et),Form("dphi  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_ph[et] -> Divide(5,2);

			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				// Momentum resolutions
				c_fits_p [et] -> cd(p+1);
				h1_dpp_p_et_bins[et][p] -> Draw();	h1_dpp_p_et_bins[et][p] -> Fit(Form("f_gaus_dpp_%i_%i",et,p),"RQ");
				width_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(2);
				error_dpp[et][p] = (f_gaus_dpp[et][p] -> GetParError(2))*(f_gaus_dpp[et][p] -> GetChisquare())/(f_gaus_dpp[et][p] -> GetNDF());

				// Theta resolution
				c_fits_th[et] -> cd(p+1);
				h1_dth_p_et_bins[et][p] -> Draw();	h1_dth_p_et_bins[et][p] -> Fit(Form("f_gaus_dth_%i_%i",et,p),"RQ");
				width_dth[et][p] = f_gaus_dth[et][p] -> GetParameter(2);
				error_dth[et][p] = (f_gaus_dth[et][p] -> GetParError(2))*(f_gaus_dth[et][p] -> GetChisquare())/(f_gaus_dth[et][p] -> GetNDF());

				// Phi resolution
				c_fits_ph[et] -> cd(p+1);
				h1_dph_p_et_bins[et][p] -> Draw();	h1_dph_p_et_bins[et][p] -> Fit(Form("f_gaus_dph_%i_%i",et,p),"RQ");
				width_dph[et][p] = f_gaus_dph[et][p] -> GetParameter(2);
				error_dph[et][p] = (f_gaus_dph[et][p] -> GetParError(2))*(f_gaus_dph[et][p] -> GetChisquare())/(f_gaus_dph[et][p] -> GetNDF());

				// ----

				h1_dpp_v_p_et_bins[et] -> SetBinContent(p +1,width_dpp[et][p]*100. );
				h1_dth_v_p_et_bins[et] -> SetBinContent(p +1,width_dth[et][p]*1000.);
				h1_dph_v_p_et_bins[et] -> SetBinContent(p +1,width_dph[et][p]*1000.);
				h1_dpp_v_et_p_bins[ p] -> SetBinContent(et+1,width_dpp[et][p]*100. );
				h1_dth_v_et_p_bins[ p] -> SetBinContent(et+1,width_dth[et][p]*1000.);
				h1_dph_v_et_p_bins[ p] -> SetBinContent(et+1,width_dph[et][p]*1000.);

				h1_dpp_v_p_et_bins[et] -> SetBinError  (p +1,error_dpp[et][p]*100. );
				h1_dth_v_p_et_bins[et] -> SetBinError  (p +1,error_dth[et][p]*1000.);
				h1_dph_v_p_et_bins[et] -> SetBinError  (p +1,error_dph[et][p]*1000.);
				h1_dpp_v_et_p_bins[ p] -> SetBinError  (et+1,error_dpp[et][p]*100. );
				h1_dth_v_et_p_bins[ p] -> SetBinError  (et+1,error_dth[et][p]*1000.);
				h1_dph_v_et_p_bins[ p] -> SetBinError  (et+1,error_dph[et][p]*1000.);
			}
			c_fits_p [et] -> Modified();	c_fits_p [et] -> Update();
			c_fits_th[et] -> Modified();	c_fits_th[et] -> Update();
			c_fits_ph[et] -> Modified();	c_fits_ph[et] -> Update();
		}
		// -------------------------------------------------------------
		// Updating table with width values
		ofstream updated_tab;
		if(update_tab){
			updated_tab.open(tab_name);
			for(int et = 0 ; et < size_eta_bin-1 ; et++){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					updated_tab << width_dpp[et][p];
					if(p == size_mom_bin-2) updated_tab << "\n";
					else updated_tab << "\t";
				}
			}
			updated_tab << "\n";
			for(int et = 0 ; et < size_eta_bin-1 ; et++){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					updated_tab << width_dth[et][p];
					if(p == size_mom_bin-2) updated_tab << "\n";
					else updated_tab << "\t";
				}
			}
			updated_tab << "\n";
			for(int et = 0 ; et < size_eta_bin-1 ; et++){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					updated_tab << width_dph[et][p];
					if(p == size_mom_bin-2) updated_tab << "\n";
					else updated_tab << "\t";
				}
			}
			updated_tab.close();
		}

		h1_dpp_v_p_et_bins[0] -> SetMinimum(0.3 );	h1_dpp_v_p_et_bins[0] -> GetYaxis() -> SetMoreLogLabels();	h1_dpp_v_p_et_bins[0]->GetYaxis()->SetTitleOffset(2.3);
		h1_dth_v_p_et_bins[0] -> SetMinimum(0.03);      h1_dth_v_p_et_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dth_v_p_et_bins[0]->GetYaxis()->SetTitleOffset(2.3);
		h1_dph_v_p_et_bins[0] -> SetMinimum(0.06);      h1_dph_v_p_et_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dph_v_p_et_bins[0]->GetYaxis()->SetTitleOffset(2.3);
		h1_dpp_v_et_p_bins[0] -> SetMinimum(0.3 );      h1_dpp_v_et_p_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dpp_v_et_p_bins[0]->GetYaxis()->SetTitleOffset(2.3);
		h1_dth_v_et_p_bins[0] -> SetMinimum(0.03);      h1_dth_v_et_p_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dth_v_et_p_bins[0]->GetYaxis()->SetTitleOffset(2.3);
		h1_dph_v_et_p_bins[0] -> SetMinimum(0.06);      h1_dph_v_et_p_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dph_v_et_p_bins[0]->GetYaxis()->SetTitleOffset(2.3);

		h1_dpp_v_p_et_bins[0] -> SetMaximum(20);
		h1_dth_v_p_et_bins[0] -> SetMaximum(3 );
		h1_dph_v_p_et_bins[0] -> SetMaximum(30);
		h1_dpp_v_et_p_bins[0] -> SetMaximum(10.1);
		h1_dth_v_et_p_bins[0] -> SetMaximum(4   );
		h1_dph_v_et_p_bins[0] -> SetMaximum(40  );

		// -------------------------------------------------------------
		// Plotting histograms
		TCanvas * c1 = new TCanvas("c1","c1",1300,900);
		c1 -> Divide(3,2);

		c1 -> cd(1); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
		h1_dpp_v_p_et_bins[0] -> Draw();
		for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dpp_v_p_et_bins[et] -> Draw("same");
		c1 -> cd(2); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
		h1_dth_v_p_et_bins[0] -> Draw();
		for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dth_v_p_et_bins[et] -> Draw("same");
		c1 -> cd(3); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
		h1_dph_v_p_et_bins[0] -> Draw();
		for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dph_v_p_et_bins[et] -> Draw("same");
		c1 -> cd(4); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
		h1_dpp_v_et_p_bins[0] -> Draw();
		for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dpp_v_et_p_bins[p] -> Draw("same");
		c1 -> cd(5); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
		h1_dth_v_et_p_bins[0] -> Draw();
		for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dth_v_et_p_bins[p] -> Draw("same");
		c1 -> cd(6); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
		h1_dph_v_et_p_bins[0] -> Draw();
		for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dph_v_et_p_bins[p] -> Draw("same");
		TLegend * leg1 = new TLegend(0.50,0.6,0.95,0.95);
		leg1 -> SetLineColor(0);
		for(int et = 0 ; et < size_eta_bin-1 ; et++) leg1 -> AddEntry(h1_dph_v_p_et_bins[et],Form("%.1f < |#eta| < %.1f",eta_bin[et],eta_bin[et+1]));
		c1 -> cd(2); 
		leg1 -> Draw("same");
		TLegend * leg2 = new TLegend(0.20,0.5,0.65,0.95);
		leg2 -> SetLineColor(0);
		for(int p = 0 ; p < size_mom_bin-1 ; p++) leg2 -> AddEntry(h1_dph_v_et_p_bins[p],Form("%.1f < p < %.1f GeV/c",mom_bin[p],mom_bin[p+1]));
		c1 -> cd(4);
		leg2 -> Draw("same");
		c1 -> Modified();
		c1 -> Update();

		// -------------------------------------------------------------
		// Saving fits to pdf
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			TString fname = out_pdf;
			if(et == 0) fname+="(";
			else if(et == size_eta_bin-2) fname+=")";
			c_fits_p [et] -> Print("fits/dpp_"+fname);
			c_fits_th[et] -> Print("fits/dth_"+fname);
			c_fits_ph[et] -> Print("fits/dph_"+fname);
		}

		// -------------------------------------------------------------
		// Saving histograms
		TFile * Fout = new TFile(outfile,"recreate");
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			h1_dpp_v_p_et_bins[et] -> Write(Form("h1_dpp_v_p_et_bins_%i",et));
			h1_dth_v_p_et_bins[et] -> Write(Form("h1_dth_v_p_et_bins_%i",et));
			h1_dph_v_p_et_bins[et] -> Write(Form("h1_dph_v_p_et_bins_%i",et));
		}
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			h1_dpp_v_et_p_bins[p] -> Write(Form("h1_dpp_v_et_p_bins_%i",p));
			h1_dth_v_et_p_bins[p] -> Write(Form("h1_dth_v_et_p_bins_%i",p));
			h1_dph_v_et_p_bins[p] -> Write(Form("h1_dph_v_et_p_bins_%i",p));
		}
		c1 -> Write("c1");
		TVT_eta_bin.Write("TVT_eta_bin");
		TVT_mom_bin.Write("TVT_mom_bin");
		Fout -> Close();

		// -------------------------------------------------------------
		// Saving plots
		c1 -> Print(out_pdf2);

		myapp -> Run();
		return 0;
	}
	// ============================================================================================================================================
	void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max ){
		h1 -> SetLineWidth(2);
		h1 -> SetLineColor(color);
		h1 -> SetMarkerStyle(marker);
		h1 -> SetMarkerColor(color);

		h1 -> SetMinimum(min);
		h1 -> SetMaximum(max);

		h1 -> GetXaxis() -> CenterTitle();
		h1 -> GetXaxis() -> SetNdivisions(107); // to draw less tick marks
		h1 -> GetYaxis() -> CenterTitle();
		h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks

		h1 -> SetMinimum(0.001);
	}
