#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <filesystem>

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
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
void better_yaxis(TH1F ** h1_array,int array_size);
TF1 * double_gaus(TH1F *h1,float min1, float max1, float min2, float max2,TString type ,int et, int p);
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	if(argc!=5){
		cout << "Run this code as:\n\033[1;32m./analysis_momentum_resolution A B C filename.root\033[0m\n";
		cout << "where:\nA = 1 -> Widths from table will be used\n  = 2 -> Widths from table \033[1;31mwon't\033[0m be used\n";
		cout << "B = 1 -> Table will be updated\n  = 2 -> Table \033[1;31mwon't\033[0m be updated\n";
		cout << "C = 1 -> Run code and quit\n  = 2 -> Run code and show plots\n";
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

	if     (atoi(argv[1])==1){use_widths = true ;	cout << "Will be using widths from table\n" ;}
	else if(atoi(argv[1])==2){use_widths = false;	cout << "Won't be using widths from table\n";}
	else{cout << "Something wrong with your election of input parameter 'A'. Bailing out!\n"; exit(0);}

	if     (atoi(argv[2])==1){update_tab = true ;   cout << "Table will be updated\n" ;}
	else if(atoi(argv[2])==2){update_tab = false;   cout << "Table won't be updated\n";}
	else{cout << "Something wrong with your election of input parameter 'B'. Bailing out!\n"; exit(0);}

	if     (atoi(argv[3])==1){keep_plots = false;	cout << "Will run and quit. Examine the output files for resulting plots\n";gROOT->SetBatch(kTRUE);}
	else if(atoi(argv[3])==2){keep_plots = true ;	cout << "Will run and show the plots\n" ; }
	else{cout << "Something wrong with your election of input parameter 'C'. Bailing out!\n"; exit(0);}

	// -------------------------
	// Binning
	float eta_bin[] = {-2.0,-1.0,0.0,1.0,2.0,4.0};
	float mom_bin[] = {4.,5.,6.,8.,10.,12.,20.};

	const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	const int size_mom_bin = sizeof(mom_bin)/sizeof(*mom_bin);

	TVectorT<double> TVT_eta_bin(size_eta_bin);	for(int i = 0 ; i < size_eta_bin ; i++) TVT_eta_bin[i] = eta_bin[i];
	TVectorT<double> TVT_mom_bin(size_mom_bin);	for(int i = 0 ; i < size_mom_bin ; i++) TVT_mom_bin[i] = mom_bin[i];
	// -------------------------
	// useful strings
	string raw_fname = argv[4];
	TString infile = "../data/" + raw_fname;
	raw_fname.resize(raw_fname.size()-5);
	TString outfile = "../output/output_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".root";
	TString tab_name = "tables/tab_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".txt";
	TString out_pdf = "output_fits_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
	TString out_pdf2 = "results/results_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
	// -------------------------------------------------------------
	// Some settings
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);	
	// -------------------------------------------------------------
	// Loading all the needed info from the root file
	TFile * F = new TFile(infile);

	TTreeReader Tree("T",F);
	int nEntries = Tree.GetEntries();
	TTreeReaderValue<int> njets(Tree,"njets");
	TTreeReaderArray<Int_t> NConst(Tree,"nComponent");
	TTreeReaderArray<Float_t> E(Tree,"e");
	TTreeReaderArray<Float_t> Eta(Tree,"eta");
	TTreeReaderArray<Float_t> Phi(Tree,"phi");
	TTreeReaderArray<Float_t> Pt(Tree,"pt");
	
	TTreeReaderArray<Float_t> gE(Tree,"matched_charged_truthE");
	TTreeReaderArray<Float_t> gEta(Tree,"matched_charged_truthEta");
	TTreeReaderArray<Float_t> gPhi(Tree,"matched_charged_truthPhi");
	TTreeReaderArray<Float_t> gPt(Tree,"matched_charged_truthPt");
	
	
	// -------------------------------------------------------------
	fstream tab;
	float approx_sig_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};

	float approx_mean_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_mean_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_mean_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	
	TString temp_str;
	if(use_widths){
		tab.open(tab_name);
		if(!tab){cout << "Could not find file '" << tab_name << "'" << endl; use_widths = false; update_tab = true;}
		else{
			cout << "Loading parameters from file '" << tab_name << "'" << endl;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dpp[et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dth[et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dph[et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_mean_dpp[et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_mean_dth[et][p];}}
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_mean_dph[et][p];}}
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
				h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dE/E;Counts",40,
							  approx_mean_dpp[et][p]-approx_sig_dpp_3_0[et][p],approx_mean_dpp[et][p]+approx_sig_dpp_3_0[et][p]);
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",40,
							  approx_mean_dth[et][p]-approx_sig_dth_3_0[et][p],approx_mean_dth[et][p]+approx_sig_dth_3_0[et][p]);
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts",40,
							  approx_mean_dph[et][p]-approx_sig_dph_3_0[et][p],approx_mean_dph[et][p]+approx_sig_dph_3_0[et][p]);
			}
			else{
			  h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dE/E;Counts"         ,40,-0.5  ,1.5  );
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",40,-0.2,0.2);
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts"  ,40,-0.2  ,0.2  );
}

			h1_dpp_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < E < %.1f GeV",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dth_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < E < %.1f GeV",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dph_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < E < %.1f GeV",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
		}
	}
	// -------------------------------------------------------------	
	TH1F ** h1_dpp_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dth_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dph_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH2F * mom_response = new TH2F("mom_response","P_{Charge}^{Reco} vs. P_{Charge}^{Truth}",50,0,50,0,50,0);
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_v_p_et_bins[et] = new TH1F(Form("h1_dpp_v_p_et_bins_%i",et),";E [GeV];dE/E [%]"      ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dpp_v_p_et_bins[et] , 50+et , 20 , 1. , 100. );
		h1_dth_v_p_et_bins[et] = new TH1F(Form("h1_dth_v_p_et_bins_%i",et),";E [GeV];d#theta [rad]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dth_v_p_et_bins[et] , 50+et , 20 , 0.001 , 1.  );
		h1_dph_v_p_et_bins[et] = new TH1F(Form("h1_dph_v_p_et_bins_%i",et),";E [GeV];d#phi [rad]"  ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dph_v_p_et_bins[et] , 50+et , 20 , 0.001 , 1. );
	}

	TH1F ** h1_dpp_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dth_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dph_v_et_p_bins = new TH1F*[size_mom_bin-1];

	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dpp_v_et_p_bins[p] = new TH1F(Form("h1_dpp_v_et_p_bins_%i",p),";#eta;dE/E [%]"      ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dpp_v_et_p_bins[p] , 50+p , 20 , 1. , 100. );
		h1_dth_v_et_p_bins[p] = new TH1F(Form("h1_dth_v_et_p_bins_%i",p),";#eta;d#theta [rad]",size_eta_bin-1,eta_bin);	prettyTH1F( h1_dth_v_et_p_bins[p] , 50+p , 20 , 0.001 , 1.  );
		h1_dph_v_et_p_bins[p] = new TH1F(Form("h1_dph_v_et_p_bins_%i",p),";#eta;d#phi [rad]"  ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dph_v_et_p_bins[p] , 50+p , 20 , 0.001 , 1000. );
	}
	cout << "\033[1;31m********************************************************************\033[0m\n";
	// -------------------------------------------------------------
	// Loop over entries of the tree	
	int ev = 0;
	while ( Tree.Next() ){
	  if (ev%10000==0) fprintf(stderr,"%d: Entry %i out of %d\n",__LINE__,ev,nEntries);
	  if (ev==500000) break;
	  for (int n = 0; n < *njets; ++n) {

	    if (NConst[n] < 3) continue;
	    if (isnan(gE[n])) continue;

	    ROOT::Math::PtEtaPhiEVector Lorentz(Pt[n],Eta[n],Phi[n],E[n]);
	    ROOT::Math::PtEtaPhiEVector gLorentz(gPt[n],gEta[n],gPhi[n],gE[n]);
		  
	    float dth = Lorentz.Theta() - gLorentz.Theta();
	    
	    float geta = gLorentz.Eta();
	    
	    float E_reco = Lorentz.E();
	    float E_truth = gLorentz.E();
	    float dE_E = (E_reco-E_truth)/E_truth;
	    //float dE_E = (E_reco)/E_truth;		
	    mom_response->Fill(gLorentz.E(),Lorentz.E());
	    
	    float dph = Lorentz.Phi() - gLorentz.Phi();

	    //cout<<"Delta Theta = "<< dth <<" dE_E = " << dE_E << " Delta Phi = " << dph << endl;
	    
	    // Filling histograms
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			if( geta >  eta_bin[et] &&  geta <= eta_bin[et+1] ){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					if( E_truth > mom_bin[p] && E_truth <= mom_bin[p+1] ){
					  //if(abs(dph) < 0.05) //better dE/E with this cut, but not physical
					    h1_dpp_p_et_bins[et][p] -> Fill( dE_E );
					  h1_dth_p_et_bins[et][p] -> Fill( dth  );
					  h1_dph_p_et_bins[et][p] -> Fill( dph  );
					}	
				}
			}
		}
	  }
	  ++ev;
	}
	cout << "\033[1;31m********************************************************************\033[0m\n";
	// -------------------------------------------------------------
	// Declaring other useful variables and functions
	float width_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float mean_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float mean_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float mean_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	
	TF1 *** f_gaus_dpp = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dth = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dph = new TF1**[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		f_gaus_dpp[et] = new TF1*[size_mom_bin-1];
		f_gaus_dth[et] = new TF1*[size_mom_bin-1];
		f_gaus_dph[et] = new TF1*[size_mom_bin-1];

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
		  if(use_widths){
		    
		    double par[6];
		    TF1 *G1 = new TF1 (Form("G1_%i_%i",et,p),"landau",-0.2,0.1);
		    TF1 *G2 = new TF1 (Form("G2_%i_%i",et,p),"gaus",approx_mean_dpp[et][p]-approx_sig_dpp_3_0[et][p],
				       approx_mean_dpp[et][p]+approx_sig_dpp_3_0[et][p]);
		    
		    h1_dpp_p_et_bins[et][p]->Fit(G1,"R");
		    h1_dpp_p_et_bins[et][p]->Fit(G2,"R");
		    
		    G1->GetParameters(&par[0]);
		    G2->GetParameters(&par[3]);
		    
		    f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus(0)+gaus(3)",
						approx_mean_dpp[et][p]-approx_sig_dpp_3_0[et][p],approx_mean_dpp[et][p]+approx_sig_dpp_3_0[et][p]);

		    f_gaus_dpp[et][p]->SetParameters(par);
		    
		    f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",
						approx_mean_dth[et][p]-approx_sig_dth_1_1[et][p],approx_mean_dth[et][p]+approx_sig_dth_1_1[et][p]);
		    f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",
						approx_mean_dph[et][p]-approx_sig_dph_1_1[et][p],approx_mean_dph[et][p]+approx_sig_dph_1_1[et][p]);
		  }
		  else{

		    //double_gaus args: histo, min1, max1, min2, max2, eta bin, p bin
		    f_gaus_dpp[et][p] = double_gaus(h1_dpp_p_et_bins[et][p],-0.2,0.1,-0.5,1.5,TString("dpp"),et,p);
		    //f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"landau",-0.5 ,0l.5 );
			  //f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"[0]*TMath::Landau(x,[1],[2])",-0.5,1.5);
		    f_gaus_dth[et][p] = double_gaus(h1_dth_p_et_bins[et][p],-0.04,0.04,-0.2,0.2,TString("dth"),et,p);
		    f_gaus_dph[et][p] = double_gaus(h1_dph_p_et_bins[et][p],-0.04,0.04,-0.2,0.2,TString("dph"),et,p);
			  // f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-0.2,0.2);
			  // f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-0.2  ,0.2  );
		  }
		}
	}

	// -------------------------------------------------------------
	// Doing fits
	TCanvas ** c_fits_p  = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_th = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_ph = new TCanvas*[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		c_fits_p [et] = new TCanvas(Form("c_fits_p_%i" ,et),Form("dE/E  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_p [et] -> Divide(5,2);
		c_fits_th[et] = new TCanvas(Form("c_fits_th_%i",et),Form("dtheta, %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_th[et] -> Divide(5,2);
		c_fits_ph[et] = new TCanvas(Form("c_fits_ph_%i",et),Form("dphi  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_ph[et] -> Divide(5,2);

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			// Momentum resolutions
			c_fits_p [et] -> cd(p+1);
			h1_dpp_p_et_bins[et][p] -> Draw();	//h1_dpp_p_et_bins[et][p] -> Fit(Form("f_gaus_dpp_%i_%i",et,p),"RQ");
			h1_dpp_p_et_bins[et][p] -> Fit(f_gaus_dpp[et][p],"RQ");
			width_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(2);
			error_dpp[et][p] = (f_gaus_dpp[et][p] -> GetParError(2));//*(f_gaus_dpp[et][p] -> GetChisquare())/(f_gaus_dpp[et][p] -> GetNDF());
			mean_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(1);
			
			// Theta resolution
			c_fits_th[et] -> cd(p+1);
			h1_dth_p_et_bins[et][p] -> Draw();	h1_dth_p_et_bins[et][p] -> Fit(Form("f_gaus_dth_%i_%i",et,p),"RQ");
			width_dth[et][p] = f_gaus_dth[et][p] -> GetParameter(2);
			error_dth[et][p] = (f_gaus_dth[et][p] -> GetParError(2));//*(f_gaus_dth[et][p] -> GetChisquare())/(f_gaus_dth[et][p] -> GetNDF());
			mean_dth[et][p] = f_gaus_dth[et][p] -> GetParameter(1);
			
			// Phi resolution
			c_fits_ph[et] -> cd(p+1);
			h1_dph_p_et_bins[et][p] -> Draw();	h1_dph_p_et_bins[et][p] -> Fit(Form("f_gaus_dph_%i_%i",et,p),"RQ");
			width_dph[et][p] = f_gaus_dph[et][p] -> GetParameter(2);
			error_dph[et][p] = (f_gaus_dph[et][p] -> GetParError(2));//*(f_gaus_dph[et][p] -> GetChisquare())/(f_gaus_dph[et][p] -> GetNDF());
			mean_dph[et][p] = f_gaus_dph[et][p] -> GetParameter(1);
			cout<<__LINE__<<": dph width = "<<width_dph[et][p]<<endl;
			cout<<__LINE__<<": dpp width = "<<width_dpp[et][p]<<endl;
			cout<<__LINE__<<": dth width = "<<width_dth[et][p]<<endl;
			// ----
			if(h1_dpp_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dpp_v_p_et_bins[et] -> SetBinContent(p +1,width_dpp[et][p]*100. );
				h1_dpp_v_p_et_bins[et] -> SetBinError  (p +1,error_dpp[et][p]*100. );
				cout<<"Bin set to "<<h1_dpp_v_p_et_bins[et] -> GetBinContent(p +1)<<endl;;
			}
			if(h1_dth_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dth_v_p_et_bins[et] -> SetBinContent(p +1,width_dth[et][p]);
				h1_dth_v_p_et_bins[et] -> SetBinError  (p +1,error_dth[et][p]);
			}
			if(h1_dph_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dph_v_p_et_bins[et] -> SetBinContent(p +1,width_dph[et][p]);
				h1_dph_v_p_et_bins[et] -> SetBinError  (p +1,error_dph[et][p]);
			}
			if(h1_dpp_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dpp_v_et_p_bins[ p] -> SetBinContent(et+1,width_dpp[et][p]*100. );
				h1_dpp_v_et_p_bins[ p] -> SetBinError  (et+1,error_dpp[et][p]*100. );
			}
			if(h1_dth_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dth_v_et_p_bins[ p] -> SetBinContent(et+1,width_dth[et][p]);
				h1_dth_v_et_p_bins[ p] -> SetBinError  (et+1,error_dth[et][p]);
			}
			if(h1_dph_p_et_bins[et][p]->GetMaximum()>50.){
				h1_dph_v_et_p_bins[ p] -> SetBinContent(et+1,width_dph[et][p]);
				h1_dph_v_et_p_bins[ p] -> SetBinError  (et+1,error_dph[et][p]);
			}

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
		//widths
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
		//Means
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << mean_dpp[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << mean_dth[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << mean_dph[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab.close();
	}
	// -------------------------------------------------------------
        // Manipulate Histogram Y-Axis 
        better_yaxis(h1_dpp_v_p_et_bins,size_eta_bin-1);
	better_yaxis(h1_dth_v_p_et_bins,size_eta_bin-1);
	better_yaxis(h1_dph_v_p_et_bins,size_eta_bin-1);
        better_yaxis(h1_dpp_v_et_p_bins,size_mom_bin-1);
        better_yaxis(h1_dth_v_et_p_bins,size_mom_bin-1);
        better_yaxis(h1_dph_v_et_p_bins,size_mom_bin-1);
	
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
	for(int p = 0 ; p < size_mom_bin-1 ; p++) leg2 -> AddEntry(h1_dph_v_et_p_bins[p],Form("%.1f < E < %.1f GeV",mom_bin[p],mom_bin[p+1]));
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
	mom_response->Write();
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

	// h1 -> SetMinimum(min);
	// h1 -> SetMaximum(max);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetNdivisions(107); // to draw less tick marks
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks

	//h1 -> SetMinimum(0.001);
}

void better_yaxis(TH1F ** h1_array,int array_size){

  //-----------------------------------------------
  //Determine the Y-axis min and max
  float min = 999;
  float max = -999;

  for(int i = 0 ; i < array_size; i++){
    int min_bin = h1_array[i]->GetMinimumBin();
    int max_bin = h1_array[i]->GetMaximumBin();
    
    float temp_min = h1_array[i]->GetBinContent(min_bin)
                     - h1_array[i]->GetBinError(min_bin);

    if (temp_min < 0)
      temp_min = h1_array[i]->GetBinContent(min_bin);

    float temp_max = h1_array[i]->GetBinContent(max_bin)
                     + h1_array[i]->GetBinError(max_bin);

    if ((temp_min > 0) && (min > temp_min))
      min = temp_min;
    if (max < temp_max){
      max = temp_max;
    }
   }
  cout<<"Minimum = "<<min<<", Maximum = "<<max<<endl;
  max = 1.2*max;
  min = 0.8*min;
  //-----------------------------------------------
  //Set the Y-axis Range
  for(int i = 0 ; i < array_size; i++){
    h1_array[i]->GetYaxis()->SetRangeUser(min,max);
  }
  //-----------------------------------------------
  //Labels
  h1_array[0]->GetYaxis()->SetTitleOffset(2.3);
  h1_array[0]->GetYaxis()->SetMoreLogLabels();

  return;
}

TF1 * double_gaus(TH1F *h1,float min1, float max1, float min2, float max2,TString type ,int et,int p){
  //First range is for narrow peak. Second range is full fit range desired.
  Double_t par[6];
  TF1 *G1 = new TF1 (Form("G1_%i_%i",et,p),"gaus",min1,max1);
  TF1 *G2 = new TF1 (Form("G2_%i_%i",et,p),"gaus",min2,max2);
  h1->Fit(G1,"R"); G1->GetParameters(&par[0]);
  h1->Fit(G2,"R"); G2->GetParameters(&par[3]);
  TF1 *double_gaus = new TF1(Form("f_gaus_%s_%i_%i",type.Data(),et,p),"gaus(0)+gaus(3)",min2,max2);
  double_gaus->SetParameters(par);
  return double_gaus;
}
