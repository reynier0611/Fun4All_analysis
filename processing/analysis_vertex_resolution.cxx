#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <Riostream.h>

// Root includes
#include <TROOT.h>
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

	if(argc!=6){
		cout << "Run this code as:\n\033[1;32m./analysis_vertex_resolution A B C filename.root config.txt\033[0m\n";
		cout << "where:\nA = 1 -> Widths from table will be used\n  = 2 -> Widths from table \033[1;31mwon't\033[0m be used\n";
		cout << "B = 1 -> Table will be updated\n  = 2 -> Table \033[1;31mwon't\033[0m be updated\n";
		cout << "C = 1 -> Run code and quit\n  = 2 -> Run code and show plots\n";
		cout << "The files 'filename.root' and 'config.txt' are expected to be in the directories '../data' and 'config' respectively\n";
		exit(0);
	}

	bool use_widths = true;
	bool update_tab = false;
	bool keep_plots = false;

	if (!fs::is_directory("../data"  ) || !fs::exists("../data"  )) fs::create_directory("../data"  ); // Create directory if it does not exist
	if (!fs::is_directory("tables"   ) || !fs::exists("tables"   )) fs::create_directory("tables"   );
	if (!fs::is_directory("../output") || !fs::exists("../output")) fs::create_directory("../output");
	if (!fs::is_directory("fits"     ) || !fs::exists("fits"     )) fs::create_directory("fits"     );
	if (!fs::is_directory("results"  ) || !fs::exists("results"  )) fs::create_directory("results"  );

	cout << "\033[1;31m********************************************************************\nUSEFUL INFO:\033[0m\nWill be loading data from file: '" << argv[4] << "' assumed to be in directory 'data'" << endl;

	cout << "Will be loading configuration from file: '" << argv[5] << "' assumed to be in directory 'config'" << endl;

	if     (atoi(argv[1])==1){use_widths = true ;	cout << "Will be using widths from table\n" ;}
	else if(atoi(argv[1])==2){use_widths = false;	cout << "Won't be using widths from table\n";}
	else{cout << "Something wrong with your election of input parameter 'A'. Bailing out!\n"; exit(0);}

	if     (atoi(argv[2])==1){update_tab = true ;   cout << "Table will be updated\n" ;}
	else if(atoi(argv[2])==2){update_tab = false;   cout << "Table won't be updated\n";}
	else{cout << "Something wrong with your election of input parameter 'B'. Bailing out!\n"; exit(0);}

	if     (atoi(argv[3])==1){keep_plots = false;   cout << "Will run and quit. Examine the output files for resulting plots\n"; gROOT->SetBatch(kTRUE);}
        else if(atoi(argv[3])==2){keep_plots = true ;   cout << "Will run and show the plots\n" ;}
        else{cout << "Something wrong with your election of input parameter 'C'. Bailing out!\n"; exit(0);}

	// -------------------------
	// Binning
	float eta_bin[100] = {0};
        float mom_bin[100] = {0};
        int ctr_eta = 0;
        int ctr_mom = 0;
        fstream f_conf;
        double val = 0;
        string sval;
	string config_filename = argv[5];
        f_conf.open("config/"+config_filename);
        if(!f_conf){cout << "Could not find config file: '" << config_filename << "'. Bailing out!" <<  endl; exit(0);}
	bool reached_end_of_line = false;
        while(!(f_conf.eof())){
                f_conf >> sval;
                if(sval=="\\") reached_end_of_line = true;
                else{
                        val = stod(sval);
                        if(!reached_end_of_line){
                                eta_bin[ctr_eta] = val;
                                ctr_eta++;
                        }
                        else{
                                mom_bin[ctr_mom] = val;
                                ctr_mom++;
                        }
                }
        }
        ctr_mom--;
	f_conf.close();
	const int size_eta_bin = ctr_eta;
	const int size_mom_bin = ctr_mom;	

	TVectorT<double> TVT_eta_bin(size_eta_bin);	for(int i = 0 ; i < size_eta_bin ; i++) TVT_eta_bin[i] = eta_bin[i];
	TVectorT<double> TVT_mom_bin(size_mom_bin);	for(int i = 0 ; i < size_mom_bin ; i++) TVT_mom_bin[i] = mom_bin[i];
	// -------------------------
	// useful strings
	string raw_fname = argv[4];
	TString infile = "../data/" + raw_fname;
	raw_fname.resize(raw_fname.size()-5);
	TString outfile = "../output/output_vtx_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".root";
	TString tab_name = "tables/tab_vtx_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".txt";
	TString out_pdf = "output_fits_vtx_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
	TString out_pdf2 = "results/results_vtx_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
	// -------------------------------------------------------------
	// Some settings
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);	
	// -------------------------------------------------------------
	// Loading all the needed info from the root file
	TFile * F = new TFile(infile);
	TTree * T = (TTree*) F -> Get("tracks");
	float gpx, gpy, gpz, px, py, pz, pcaz, dca2d, gvz;
	T -> SetBranchAddress("gpx"  ,&gpx  );
	T -> SetBranchAddress("gpy"  ,&gpy  );
	T -> SetBranchAddress("gpz"  ,&gpz  );
	T -> SetBranchAddress("px"   ,&px   );
	T -> SetBranchAddress("py"   ,&py   );
	T -> SetBranchAddress("pz"   ,&pz   ); 
	T -> SetBranchAddress("pcaz" ,&pcaz );
    	T -> SetBranchAddress("dca2d",&dca2d);
	T -> SetBranchAddress("gvz"  ,&gvz  );
	int nEntries = T -> GetEntries();
	// -------------------------------------------------------------
	fstream tab;
	float approx_sig_dvl_p [100][100] = {{0}};
	float approx_sig_dvt_p [100][100] = {{0}};
	float approx_sig_dvl_pT[100][100] = {{0}};
	float approx_sig_dvt_pT[100][100] = {{0}};
	TString temp_str;
	if(use_widths){
		tab.open(tab_name);
		if(!tab){cout << "Could not find file '" << tab_name << "'" << endl; use_widths = false; update_tab = true;}
		else{
			cout << "Loading parameters from file '" << tab_name << "'" << endl;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dvl_p [et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dvt_p [et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dvl_pT[et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dvt_pT[et][p];}}
		}
		tab.close();
	}

	float approx_sig_dvl_p_3_0 [100][100] = {{0}};	float approx_sig_dvl_p_1_1 [100][100] = {{0}};
	float approx_sig_dvt_p_3_0 [100][100] = {{0}};	float approx_sig_dvt_p_1_1 [100][100] = {{0}};
	float approx_sig_dvl_pT_3_0[100][100] = {{0}};	float approx_sig_dvl_pT_1_1[100][100] = {{0}};
	float approx_sig_dvt_pT_3_0[100][100] = {{0}};	float approx_sig_dvt_pT_1_1[100][100] = {{0}};

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			approx_sig_dvl_p_3_0 [et][p] = 3.0*approx_sig_dvl_p [et][p];	approx_sig_dvl_p_1_1 [et][p] = 1.1*approx_sig_dvl_p [et][p];
			approx_sig_dvt_p_3_0 [et][p] = 3.0*approx_sig_dvt_p [et][p];	approx_sig_dvt_p_1_1 [et][p] = 1.1*approx_sig_dvt_p [et][p];
			approx_sig_dvl_pT_3_0[et][p] = 3.0*approx_sig_dvl_pT[et][p];	approx_sig_dvl_pT_1_1[et][p] = 1.1*approx_sig_dvl_pT[et][p];
			approx_sig_dvt_pT_3_0[et][p] = 3.0*approx_sig_dvt_pT[et][p];	approx_sig_dvt_pT_1_1[et][p] = 1.1*approx_sig_dvt_pT[et][p];
		}
	}
	// -------------------------------------------------------------
	// Defining histograms
	TH1F *** h1_dvl_p_et_bins  = new TH1F**[size_eta_bin-1];
	TH1F *** h1_dvt_p_et_bins  = new TH1F**[size_eta_bin-1];
	TH1F *** h1_dvl_pT_et_bins = new TH1F**[size_eta_bin-1];
        TH1F *** h1_dvt_pT_et_bins = new TH1F**[size_eta_bin-1];
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dvl_p_et_bins [et] = new TH1F*[size_mom_bin-1];
		h1_dvt_p_et_bins [et] = new TH1F*[size_mom_bin-1];
		h1_dvl_pT_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dvt_pT_et_bins[et] = new TH1F*[size_mom_bin-1];
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				h1_dvl_p_et_bins [et][p] = new TH1F(Form("h1_dvl_p_et_bins_%i_%i" ,et,p),";dca_{z} [cm];Counts",60,-approx_sig_dvl_p_3_0 [et][p],approx_sig_dvl_p_3_0 [et][p]);
				h1_dvt_p_et_bins [et][p] = new TH1F(Form("h1_dvt_p_et_bins_%i_%i" ,et,p),";dca_{T} [cm];Counts",60,-approx_sig_dvt_p_3_0 [et][p],approx_sig_dvt_p_3_0 [et][p]);
				h1_dvl_pT_et_bins[et][p] = new TH1F(Form("h1_dvl_pT_et_bins_%i_%i",et,p),";dca_{z} [cm];Counts",60,-approx_sig_dvl_pT_3_0[et][p],approx_sig_dvl_pT_3_0[et][p]);
				h1_dvt_pT_et_bins[et][p] = new TH1F(Form("h1_dvt_pT_et_bins_%i_%i",et,p),";dca_{T} [cm];Counts",60,-approx_sig_dvt_pT_3_0[et][p],approx_sig_dvt_pT_3_0[et][p]);
			}
			else{
				h1_dvl_p_et_bins [et][p] = new TH1F(Form("h1_dvl_p_et_bins_%i_%i" ,et,p),";dca_{z} [cm];Counts",60,-0.015,0.015);
				h1_dvt_p_et_bins [et][p] = new TH1F(Form("h1_dvt_p_et_bins_%i_%i" ,et,p),";dca_{T} [cm];Counts",60,-0.015,0.015);
				h1_dvl_pT_et_bins[et][p] = new TH1F(Form("h1_dvl_pT_et_bins_%i_%i",et,p),";dca_{z} [cm];Counts",60,-0.015,0.015);
                                h1_dvt_pT_et_bins[et][p] = new TH1F(Form("h1_dvt_pT_et_bins_%i_%i",et,p),";dca_{T} [cm];Counts",60,-0.015,0.015);
			}

			h1_dvl_p_et_bins [et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dvt_p_et_bins [et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dvl_pT_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p_{T} < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dvt_pT_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p_{T} < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
		}
	}
	// -------------------------------------------------------------
	TH1F ** h1_dvl_v_p_et_bins  = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dvt_v_p_et_bins  = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dvl_v_pT_et_bins = new TH1F*[size_eta_bin-1];
        TH1F ** h1_dvt_v_pT_et_bins = new TH1F*[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dvl_v_p_et_bins [et] = new TH1F(Form("h1_dvl_v_p_et_bins_%i" ,et),";p [GeV/c];#sigma(DCA_{z}) [#mum]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dvl_v_p_et_bins [et] , 51+7*et , 20 , 0. , 10. );
		h1_dvt_v_p_et_bins [et] = new TH1F(Form("h1_dvt_v_p_et_bins_%i" ,et),";p [GeV/c];#sigma(DCA_{T}) [#mum]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dvt_v_p_et_bins [et] , 51+7*et , 20 , 0. , 1.  );
		h1_dvl_v_pT_et_bins[et] = new TH1F(Form("h1_dvl_v_pT_et_bins_%i",et),";p_{T} [GeV/c];#sigma(DCA_{z}) [#mum]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dvl_v_pT_et_bins[et] , 51+7*et , 20 , 0. , 10. );
		h1_dvt_v_pT_et_bins[et] = new TH1F(Form("h1_dvt_v_pT_et_bins_%i",et),";p_{T} [GeV/c];#sigma(DCA_{T}) [#mum]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dvt_v_pT_et_bins[et] , 51+7*et , 20 , 0. , 1.  );
	}

	TH1F ** h1_dvl_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dvt_v_et_p_bins = new TH1F*[size_mom_bin-1];

	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dvl_v_et_p_bins[p] = new TH1F(Form("h1_dvl_v_et_p_bins_%i",p),";#eta;#sigma(DCA_{z}) [#mum]",size_eta_bin-1,eta_bin);	prettyTH1F( h1_dvl_v_et_p_bins[p] , 51+3*p , 20 , 0. , 10. );
		h1_dvt_v_et_p_bins[p] = new TH1F(Form("h1_dvt_v_et_p_bins_%i",p),";#eta;#sigma(DCA_{T}) [#mum]",size_eta_bin-1,eta_bin);	prettyTH1F( h1_dvt_v_et_p_bins[p] , 51+3*p , 20 , 0. , 1.  );
	}
	// -------------------------------------------------------------
	// Declaring other useful variables and functions
	float width_dvl_p[100][100] = {{0}};	float width_dvl_pT[100][100] = {{0}};
	float error_dvl_p[100][100] = {{0}};	float error_dvl_pT[100][100] = {{0}};
	float width_dvt_p[100][100] = {{0}};	float width_dvt_pT[100][100] = {{0}};
	float error_dvt_p[100][100] = {{0}};	float error_dvt_pT[100][100] = {{0}};

	TF1 *** f_gaus_dvl_p = new TF1**[size_eta_bin-1];	TF1 *** f_gaus_dvl_pT = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dvt_p = new TF1**[size_eta_bin-1];	TF1 *** f_gaus_dvt_pT = new TF1**[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		f_gaus_dvl_p[et] = new TF1*[size_mom_bin-1];	f_gaus_dvl_pT[et] = new TF1*[size_mom_bin-1];
		f_gaus_dvt_p[et] = new TF1*[size_mom_bin-1];	f_gaus_dvt_pT[et] = new TF1*[size_mom_bin-1];

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				f_gaus_dvl_p [et][p] = new TF1(Form("f_gaus_dvl_p_%i_%i" ,et,p),"gaus",-approx_sig_dvl_p_1_1 [et][p],approx_sig_dvl_p_1_1 [et][p]);
				f_gaus_dvt_p [et][p] = new TF1(Form("f_gaus_dvt_p_%i_%i" ,et,p),"gaus",-approx_sig_dvt_p_1_1 [et][p],approx_sig_dvt_p_1_1 [et][p]);
				f_gaus_dvl_pT[et][p] = new TF1(Form("f_gaus_dvl_pT_%i_%i",et,p),"gaus",-approx_sig_dvl_pT_1_1[et][p],approx_sig_dvl_pT_1_1[et][p]);
				f_gaus_dvt_pT[et][p] = new TF1(Form("f_gaus_dvt_pT_%i_%i",et,p),"gaus",-approx_sig_dvt_pT_1_1[et][p],approx_sig_dvt_pT_1_1[et][p]);
			}
			else{
				f_gaus_dvl_p [et][p] = new TF1(Form("f_gaus_dvl_p_%i_%i" ,et,p),"gaus",-0.01,0.01);
				f_gaus_dvt_p [et][p] = new TF1(Form("f_gaus_dvt_p_%i_%i" ,et,p),"gaus",-0.01,0.01);
				f_gaus_dvl_pT[et][p] = new TF1(Form("f_gaus_dvl_pT_%i_%i",et,p),"gaus",-0.01,0.01);
				f_gaus_dvt_pT[et][p] = new TF1(Form("f_gaus_dvt_pT_%i_%i",et,p),"gaus",-0.01,0.01);
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
		float geta = -TMath::Log(TMath::Tan(gtheta/2.));	
		float p_truth = sqrt(gpx*gpx+gpy*gpy+gpz*gpz);
		float pT_truth = sqrt(gpx*gpx+gpy*gpy);

		float dvt = dca2d;
                float dvl = pcaz - gvz;

		// Filling histograms
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			if( geta >  eta_bin[et] &&  geta <= eta_bin[et+1] ){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					if( p_truth > mom_bin[p] && p_truth <= mom_bin[p+1] ){
						h1_dvl_p_et_bins[et][p] -> Fill( dvl );
						h1_dvt_p_et_bins[et][p] -> Fill( dvt );
					}
					// transverse momentum
					if( pT_truth > mom_bin[p] && pT_truth <= mom_bin[p+1]){
						h1_dvl_pT_et_bins[et][p] -> Fill( dvl );
						h1_dvt_pT_et_bins[et][p] -> Fill( dvt );
					}	
				}
			}
		}
	}
	cout << "\033[1;31m********************************************************************\033[0m\n";
	// -------------------------------------------------------------
	// Doing fits
	TCanvas ** c_fits_dvl_p  = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_dvt_p  = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_dvl_pT = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_dvt_pT = new TCanvas*[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		c_fits_dvl_p [et] = new TCanvas(Form("c_fits_dvl_p_%i" ,et),Form("d vtx_{z} p,  %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_dvl_p [et] -> Divide(4,4);
		c_fits_dvt_p [et] = new TCanvas(Form("c_fits_dvt_p_%i" ,et),Form("d vtx_{T} p,  %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_dvt_p [et] -> Divide(4,4);
		c_fits_dvl_pT[et] = new TCanvas(Form("c_fits_dvl_pT_%i",et),Form("d vtx_{z} pT, %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_dvl_pT[et] -> Divide(4,4);
		c_fits_dvt_pT[et] = new TCanvas(Form("c_fits_dvt_pT_%i",et),Form("d vtx_{T} pT, %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_dvt_pT[et] -> Divide(4,4);

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			// Longitudinal (z) vertex resolution (p)
			c_fits_dvl_p[et] -> cd(p+1);
			h1_dvl_p_et_bins[et][p] -> Draw();	h1_dvl_p_et_bins[et][p] -> Fit(Form("f_gaus_dvl_p_%i_%i",et,p),"RQ");
			width_dvl_p[et][p] = f_gaus_dvl_p[et][p] -> GetParameter(2);
			error_dvl_p[et][p] = (f_gaus_dvl_p[et][p] -> GetParError(2))*(f_gaus_dvl_p[et][p] -> GetChisquare())/(f_gaus_dvl_p[et][p] -> GetNDF());

			// Transverse (rho) vertex resolution (p)
			c_fits_dvt_p[et] -> cd(p+1);
			h1_dvt_p_et_bins[et][p] -> Draw();	h1_dvt_p_et_bins[et][p] -> Fit(Form("f_gaus_dvt_p_%i_%i",et,p),"RQ");
			width_dvt_p[et][p] = f_gaus_dvt_p[et][p] -> GetParameter(2);
			error_dvt_p[et][p] = (f_gaus_dvt_p[et][p] -> GetParError(2))*(f_gaus_dvt_p[et][p] -> GetChisquare())/(f_gaus_dvt_p[et][p] -> GetNDF());

			// Longitudinal (z) vertex resolution (pT)
                        c_fits_dvl_pT[et] -> cd(p+1);
                        h1_dvl_pT_et_bins[et][p] -> Draw();	h1_dvl_pT_et_bins[et][p] -> Fit(Form("f_gaus_dvl_pT_%i_%i",et,p),"RQ");
                        width_dvl_pT[et][p] = f_gaus_dvl_pT[et][p] -> GetParameter(2);
                        error_dvl_pT[et][p] = (f_gaus_dvl_pT[et][p] -> GetParError(2))*(f_gaus_dvl_pT[et][p] -> GetChisquare())/(f_gaus_dvl_pT[et][p] -> GetNDF());

                        // Transverse (rho) vertex resolution (pT)
                        c_fits_dvt_pT[et] -> cd(p+1);
                        h1_dvt_pT_et_bins[et][p] -> Draw();	h1_dvt_pT_et_bins[et][p] -> Fit(Form("f_gaus_dvt_pT_%i_%i",et,p),"RQ");
                        width_dvt_pT[et][p] = f_gaus_dvt_pT[et][p] -> GetParameter(2);
                        error_dvt_pT[et][p] = (f_gaus_dvt_pT[et][p] -> GetParError(2))*(f_gaus_dvt_pT[et][p] -> GetChisquare())/(f_gaus_dvt_pT[et][p] -> GetNDF());

			// ----
			if(h1_dvl_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dvl_v_p_et_bins[et] -> SetBinContent(p +1,width_dvl_p[et][p]*10000.);
				h1_dvl_v_p_et_bins[et] -> SetBinError  (p +1,error_dvl_p[et][p]*10000.);
			}
			if(h1_dvt_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dvt_v_p_et_bins[et] -> SetBinContent(p +1,width_dvt_p[et][p]*10000.);
				h1_dvt_v_p_et_bins[et] -> SetBinError  (p +1,error_dvt_p[et][p]*10000.);
			}
			if(h1_dvl_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dvl_v_et_p_bins[ p] -> SetBinContent(et+1,width_dvl_p[et][p]*10000.);
				h1_dvl_v_et_p_bins[ p] -> SetBinError  (et+1,error_dvl_p[et][p]*10000.);
			}
			if(h1_dvt_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dvt_v_et_p_bins[ p] -> SetBinContent(et+1,width_dvt_p[et][p]*10000.);
				h1_dvt_v_et_p_bins[ p] -> SetBinError  (et+1,error_dvt_p[et][p]*10000.);
			}
			// Transverse momentum
			if(h1_dvl_pT_et_bins[et][p]->GetMaximum()>50.){
                                h1_dvl_v_pT_et_bins[et] -> SetBinContent(p+1,width_dvl_pT[et][p]*10000.);
                                h1_dvl_v_pT_et_bins[et] -> SetBinError  (p+1,error_dvl_pT[et][p]*10000.);
                        }
                        if(h1_dvt_pT_et_bins[et][p]->GetMaximum()>50.){
                                h1_dvt_v_pT_et_bins[et] -> SetBinContent(p+1,width_dvt_pT[et][p]*10000.);
                                h1_dvt_v_pT_et_bins[et] -> SetBinError  (p+1,error_dvt_pT[et][p]*10000.);
                        }
		}
		c_fits_dvl_p [et] -> Modified();	c_fits_dvl_p [et] -> Update();
		c_fits_dvt_p [et] -> Modified();	c_fits_dvt_p [et] -> Update();
		c_fits_dvl_pT[et] -> Modified();	c_fits_dvl_pT[et] -> Update();
		c_fits_dvt_pT[et] -> Modified();	c_fits_dvt_pT[et] -> Update();
	}
	// -------------------------------------------------------------
	// Updating table with width values
	ofstream updated_tab;
	if(update_tab){
		updated_tab.open(tab_name);
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dvl_p[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dvt_p[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
                for(int et = 0 ; et < size_eta_bin-1 ; et++){
                        for(int p = 0 ; p < size_mom_bin-1 ; p++){
                                updated_tab << width_dvl_pT[et][p];
                                if(p == size_mom_bin-2) updated_tab << "\n";
                                else updated_tab << "\t";
                        }
                }
                updated_tab << "\n";
                for(int et = 0 ; et < size_eta_bin-1 ; et++){
                        for(int p = 0 ; p < size_mom_bin-1 ; p++){
                                updated_tab << width_dvt_pT[et][p];
                                if(p == size_mom_bin-2) updated_tab << "\n";
                                else updated_tab << "\t";
                        }
                }
		updated_tab.close();
	}

	h1_dvl_v_p_et_bins[0] -> SetMinimum(2);	h1_dvl_v_p_et_bins[0] -> GetYaxis() -> SetMoreLogLabels();	h1_dvl_v_p_et_bins[0]->GetYaxis()->SetTitleOffset(1.5);
	h1_dvt_v_p_et_bins[0] -> SetMinimum(2);	h1_dvt_v_p_et_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dvt_v_p_et_bins[0]->GetYaxis()->SetTitleOffset(1.5);
	h1_dvl_v_et_p_bins[0] -> SetMinimum(2);	h1_dvl_v_et_p_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dvl_v_et_p_bins[0]->GetYaxis()->SetTitleOffset(1.5);
	h1_dvt_v_et_p_bins[0] -> SetMinimum(2);	h1_dvt_v_et_p_bins[0] -> GetYaxis() -> SetMoreLogLabels();      h1_dvt_v_et_p_bins[0]->GetYaxis()->SetTitleOffset(1.5);

	h1_dvl_v_p_et_bins[0] -> SetMaximum(70000);
	h1_dvt_v_p_et_bins[0] -> SetMaximum(70000);
	h1_dvl_v_et_p_bins[0] -> SetMaximum(70000);
	h1_dvt_v_et_p_bins[0] -> SetMaximum(70000);

	// -------------------------------------------------------------
	// Plotting histograms
	TCanvas * c1 = new TCanvas("c1","c1",1300,900);
	c1 -> Divide(2,2);

	c1 -> cd(1); gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.14); gPad -> SetLogy();
	h1_dvl_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dvl_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(2); gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.14); gPad -> SetLogy();
	h1_dvt_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dvt_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(3); gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.14); gPad -> SetLogy();
	h1_dvl_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dvl_v_et_p_bins[p] -> Draw("same");
	c1 -> cd(4); gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.14); gPad -> SetLogy();
	h1_dvt_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dvt_v_et_p_bins[p] -> Draw("same");
	TLegend * leg1 = new TLegend(0.20,0.75,0.90,0.95);
	leg1 -> SetLineColor(0);
	leg1 -> SetNColumns(3);
	for(int et = 0 ; et < size_eta_bin-1 ; et++) leg1 -> AddEntry(h1_dvt_v_p_et_bins[et],Form("%.1f < |#eta| < %.1f",eta_bin[et],eta_bin[et+1]));
	c1 -> cd(2);
	leg1 -> Draw("same");
	TLegend * leg2 = new TLegend(0.17,0.70,0.96,0.95);
	leg2 -> SetLineColor(0);
	leg2 -> SetNColumns(3);
	for(int p = 0 ; p < size_mom_bin-1 ; p++) leg2 -> AddEntry(h1_dvt_v_et_p_bins[p],Form("%.1f < p < %.1f GeV/c",mom_bin[p],mom_bin[p+1]));
	c1 -> cd(4);
	leg2 -> Draw("same");
	c1 -> Modified();
	c1 -> Update();

	TCanvas * c2 = new TCanvas("c2","c2",1300,900);
        c2 -> Divide(2,1);
        c2 -> cd(1); gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.14); gPad -> SetLogy();
        h1_dvl_v_pT_et_bins[0] -> Draw();
        for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dvl_v_pT_et_bins[et] -> Draw("same");
        c2 -> cd(2); gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.14); gPad -> SetLogy();
        h1_dvt_v_pT_et_bins[0] -> Draw();
        for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dvt_v_pT_et_bins[et] -> Draw("same"); 
        leg1 -> Draw("same");
        c2 -> Modified();
        c2 -> Update();

	// -------------------------------------------------------------
	// Saving fits to pdf
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		TString fname = out_pdf;
		if(et == 0) fname+="(";
		else if(et == size_eta_bin-2) fname+=")";
		c_fits_dvl_p [et] -> Print("fits/dvl_p_" +fname);
		c_fits_dvt_p [et] -> Print("fits/dvt_p_" +fname);
		c_fits_dvl_pT[et] -> Print("fits/dvl_pT_"+fname);
		c_fits_dvt_pT[et] -> Print("fits/dvt_pT_"+fname);
	}

	// -------------------------------------------------------------
	// Saving histograms
	TFile * Fout = new TFile(outfile,"recreate");
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dvl_v_p_et_bins [et] -> Write(Form("h1_dvl_v_p_et_bins_%i" ,et));
		h1_dvt_v_p_et_bins [et] -> Write(Form("h1_dvt_v_p_et_bins_%i" ,et));
		h1_dvl_v_pT_et_bins[et] -> Write(Form("h1_dvl_v_pT_et_bins_%i",et));
		h1_dvt_v_pT_et_bins[et] -> Write(Form("h1_dvt_v_pT_et_bins_%i",et));
	}
	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dvl_v_et_p_bins[p] -> Write(Form("h1_dvl_v_et_p_bins_%i",p));
		h1_dvt_v_et_p_bins[p] -> Write(Form("h1_dvt_v_et_p_bins_%i",p));
	}
	c1 -> Write("c1");
	c2 -> Write("c2");
	TVT_eta_bin.Write("TVT_eta_bin");
	TVT_mom_bin.Write("TVT_mom_bin");
	Fout -> Close();

	// -------------------------------------------------------------
	// Saving plots
	c1 -> Print(out_pdf2+"(");
	c2 -> Print(out_pdf2+")");

	if(keep_plots)
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
	h1 -> GetXaxis() -> SetTitleSize(0.05);
	h1 -> GetXaxis() -> SetLabelSize(0.05);
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks
	h1 -> GetYaxis() -> SetTitleSize(0.05);
	h1 -> GetYaxis() -> SetLabelSize(0.05);

	h1 -> SetMinimum(0.001);
}
