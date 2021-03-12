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
#include "TLatex.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVectorT.h"
#include "TGraph.h"
#include "Math/LorentzVector.h"

using namespace std;

// Forward-declaring functions
void prettyTH1D( TH1D * h1 , int color , int marker , float min , float max );
void prettyTGraph( TGraph * g , int color , int marker , float min , float max , TString xtit , TString ytit , TString tit );
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	// -------------------------------------------------------------
	// Some settings
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);	
	// -------------------------------------------------------------
	// Loading all the needed info from the root file
	TString fname[] = {
		//"out_flat_pT_pi-_det2_10x10_Beast_FastTrackingEval.root",		// Realistic geometry,  (v,b,d) = (0.30,0.30,0.30) % X0, Beast
		//"out_flat_pT_pi-_det2_10x10_sPHENIX_FastTrackingEval.root",		// Realistic geometry,  (v,b,d) = (0.30,0.30,0.30) % X0, sPHENIX
		"combined_simp_geom_v2_vbd_0.3_0.3_0.3_Beast_FastSimEval.root",		// Simplified geometry, (v,b,d) = (0.30,0.30,0.30) % X0, Beast
		"combined_simp_geom_v2_vbd_0.3_0.3_0.3_sPHENIX_FastSimEval.root",	// Simplified geometry, (v,b,d) = (0.30,0.30,0.30) % X0, sPHENIX
		"combined_simp_geom_v2_vbd_0.05_0.55_0.24_Beast_FastSimEval.root",	// Simplified geometry, (v,b,d) = (0.05,0.55,0.24) % X0, Beast
		"combined_simp_geom_v2_vbd_0.05_0.55_0.24_sPHENIX_FastSimEval.root",	// Simplified geometry, (v,b,d) = (0.05,0.55,0.24) % X0, sPHENIX
		"combined_simp_geom_3vtx_vbd_0.05_0.55_0.24_Beast_FastSimEval.root",	// Simplified geometry, (v,b,d) = (0.05,0.55,0.24) % X0, Beast  , 3 vtx layers
		"combined_simp_geom_3vtx_vbd_0.05_0.55_0.24_sPHENIX_FastSimEval.root",	// Simplified geometry, (v,b,d) = (0.05,0.55,0.24) % X0, sPHENIX, 3 vtx layers
	};
	//"combined_simp_geom_3vtx_vbd_0.3_0.3_0.3_Beast_FastSimEval.root"
	//"combined_simp_geom_3vtx_vbd_0.3_0.3_0.3_sPHENIX_FastSimEval.root"

	const int size_fname = sizeof(fname)/sizeof(*fname);
	TString B_field[] = {
		"(0.30,0.30,0.30), 3.0 T",
		"(0.30,0.30,0.30), 1.4 T",
		"(0.05,0.55,0.24), 3.0 T",
		"(0.05,0.55,0.24), 1.4 T",
		"(0.05,0.55,0.24), 3.0 T (3 vtx)",
                "(0.05,0.55,0.24), 1.4 T (3 vtx)"
	};

	TFile ** F = new TFile * [size_fname];
	TTree ** T = new TTree * [size_fname];

	float gpx[size_fname], gpy[size_fname], gpz[size_fname];
	int charge[size_fname], nhits[size_fname], trackID[size_fname], nEntries[size_fname];

	for(int i = 0 ; i < size_fname ; i++){
		F[i] = new TFile("data/"+fname[i]);
		T[i] = (TTree*) F[i] -> Get("tracks");
		T[i] -> SetBranchAddress("gpx"    ,&gpx    [i]);
		T[i] -> SetBranchAddress("gpy"    ,&gpy    [i]);
		T[i] -> SetBranchAddress("gpz"    ,&gpz    [i]);
		T[i] -> SetBranchAddress("charge" ,&charge [i]);
		T[i] -> SetBranchAddress("nhits"  ,&nhits  [i]);
		T[i] -> SetBranchAddress("trackID",&trackID[i]);
		nEntries[i] = T[i] -> GetEntries();
	}

	// -------------------------------------------------------------
	// Defining histograms for resolution distributions per p and eta bins
	TH2D ** h2_pT_eta_num_finer  = new TH2D * [size_fname];
	TH2D ** h2_pT_eta_num_coarse = new TH2D * [size_fname];
	TH2D ** h2_pT_eta_den_finer  = new TH2D * [size_fname];
	TH2D ** h2_pT_eta_den_coarse = new TH2D * [size_fname];
	TH2D ** h2_pT_eta_rat_finer  = new TH2D * [size_fname];
	TH2D ** h2_pT_eta_rat_coarse = new TH2D * [size_fname];

	//double eta_bin[] = {-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5};
	double eta_bin[] = {-3.5,-3,-2.5,-2,-1.5,-1,-0.5,-0.3,0.3,0.5,1,1.5,2,2.5,3,3.5};
	const int n_eta_bins = sizeof(eta_bin)/sizeof(*eta_bin) - 1;
	cout << "eta binning used: ";
	for(int i = 0 ; i < n_eta_bins ; i++) cout << eta_bin[i] << ", ";
	cout << eta_bin[n_eta_bins] << endl;

	for(int file = 0 ; file < size_fname ; file++){
		h2_pT_eta_num_finer [file] = new TH2D(Form("h2_pT_eta_num_finer_%i" ,file),B_field[file]+", only reconstructed;#eta;p_{T} [GeV/#it{c}]",100,-4,4,100,0,3);
		h2_pT_eta_den_finer [file] = new TH2D(Form("h2_pT_eta_den_finer_%i" ,file),B_field[file]+", all generated;#eta;p_{T} [GeV/#it{c}]"     ,100,-4,4,100,0,3);

		h2_pT_eta_num_coarse[file] = new TH2D(Form("h2_pT_eta_num_coarse_%i",file),B_field[file]+", only reconstructed;#eta;p_{T} [GeV/#it{c}]",  n_eta_bins,eta_bin,20,0,2.5);
		h2_pT_eta_den_coarse[file] = new TH2D(Form("h2_pT_eta_den_coarse_%i",file),B_field[file]+", all generated;#eta;p_{T} [GeV/#it{c}]"     ,  n_eta_bins,eta_bin,20,0,2.5);
	}

	// -------------------------------------------------------------
	// Loop over entries of the tree
	for(int file = 0 ; file < size_fname ; file++){
		cout << endl << "Looping over file: '" << fname[file] << "'\n";
		//nEntries[file] = 1000000;
		for(int ev = 0 ; ev < nEntries[file] ; ev++){
			T[file] -> GetEntry(ev);
			if(ev%1000000==0) cout << "Looping over entry " << ev << " out of " << nEntries[file] << endl;

			// Calculating some variables
			float gtheta = TMath::ACos(gpz[file]/sqrt(gpx[file]*gpx[file]+gpy[file]*gpy[file]+gpz[file]*gpz[file]));	
			float geta = -TMath::Log(TMath::Tan(gtheta/2.));
			float pt_truth = sqrt(gpx[file]*gpx[file]+gpy[file]*gpy[file]);

			// Filling histograms for tracks that pass reconstruction
			if(
					charge [file] != -9999 &&
					nhits  [file] != -9999 &&
					trackID[file] != -9999
			  ){
				h2_pT_eta_num_finer [file] -> Fill( geta , pt_truth );
				h2_pT_eta_num_coarse[file] -> Fill( geta , pt_truth );
			}

			// Filling histograms for all trackes, whether reconstructed or not
			h2_pT_eta_den_finer [file] -> Fill( geta , pt_truth );
			h2_pT_eta_den_coarse[file] -> Fill( geta , pt_truth );
		}
	}
	// -------------------------------------------------------------
	// Doing some histogram manipulation
	// ------------------------
	// Taking ratios to determine pseudo-efficiency
	for(int file = 0 ; file < size_fname ; file++){
		h2_pT_eta_rat_finer [file] = (TH2D*) h2_pT_eta_num_finer [file] -> Clone();
		h2_pT_eta_rat_coarse[file] = (TH2D*) h2_pT_eta_num_coarse[file] -> Clone();
		h2_pT_eta_rat_finer [file] -> SetNameTitle(Form("h2_pT_eta_rat_finer_%i" ,file),B_field[file]+", ratio");
		h2_pT_eta_rat_coarse[file] -> SetNameTitle(Form("h2_pT_eta_rat_coarse_%i",file),B_field[file]+", ratio");
		h2_pT_eta_rat_finer [file] -> Divide( h2_pT_eta_den_finer [file] );
		h2_pT_eta_rat_coarse[file] -> Divide( h2_pT_eta_den_coarse[file] );
	}
	// ------------------------
	// Doing 1D projections	and sigmoid fits
	//int color[] = {51,62,8,1,92,96,2,50};
	TH1D *** h1_eff_pT = new TH1D ** [size_fname];
	TF1 *** f_sigmoid = new TF1 ** [size_fname];
	for(int file = 0 ; file < size_fname ; file++){
		h1_eff_pT[file] = new TH1D * [n_eta_bins];
		f_sigmoid[file] = new TF1 * [n_eta_bins];
		for(int et = 0 ; et < n_eta_bins ; et++){
			int clr = 51 + et*50/n_eta_bins;
			cout << clr << endl;
			h1_eff_pT[file][et] = h2_pT_eta_rat_coarse[file] -> ProjectionY(Form("h1_eff_pT_%i_%i",file,et),et+1,et+1,"e");		
			prettyTH1D( h1_eff_pT[file][et] , clr , 20 , 0 , 1.1 );
			f_sigmoid[file][et] = new TF1(Form("f_sigm_%i_%i",file,et), "[0]/(1+ TMath::Exp(-[1]*(x-[2])))", 0, 3);
			f_sigmoid[file][et] -> SetParameters(1,10,.5);
			f_sigmoid[file][et] -> SetLineColor(clr);
			h1_eff_pT[file][et] -> Fit(Form("f_sigm_%i_%i",file,et));
		}
	}
	// ------------------------
	// Creating new functions with the results
	double pT_at_50percent[size_fname][n_eta_bins] = {{0}};
	double max_height     [size_fname][n_eta_bins] = {{0}};
	double avg_eta[n_eta_bins] = {0};
	for(int et = 0 ; et < n_eta_bins ; et++) avg_eta[et] = (eta_bin[et+1]+eta_bin[et])/2.;

	int color[] = {2,62,96,4,93,8,50,51};

	TGraph ** g_pTthresh_v_eta = new TGraph * [size_fname];
	for(int file = 0 ; file < size_fname ; file++){
		for(int et = 0 ; et < n_eta_bins ; et++){
			double par0 = f_sigmoid[file][et] -> GetParameter(0);
			double par1 = f_sigmoid[file][et] -> GetParameter(1);
			double par2 = f_sigmoid[file][et] -> GetParameter(2);
			pT_at_50percent[file][et] = par2 - TMath::Log(par0/0.5-1)/par1;
			pT_at_50percent[file][et] *= 1000.; // Convert from GeV to MeV
			max_height     [file][et] = f_sigmoid[file][et] -> Eval(3);
		}
		g_pTthresh_v_eta[file] = new TGraph(n_eta_bins,avg_eta, pT_at_50percent[file]);
		double max_y = 1300;
		if(file%2!=0) max_y = 800;
		prettyTGraph( g_pTthresh_v_eta[file] , color[file] , 20+file , 0 , max_y , "#eta" , "p_{T} [MeV/#it{c}]" , "" );
	}
	// ------------------------
	// Ratios
	double threshold_rat_B[size_fname/2][n_eta_bins] = {{0}};
	TGraph ** g_pTthresh_rat_v_eta = new TGraph * [size_fname/2];
	int ctr = 0;
	for(int file = 0 ; file < size_fname ; file+=2){
		for(int et = 0 ; et < n_eta_bins ; et++){
			threshold_rat_B[ctr][et] = pT_at_50percent[file][et]/pT_at_50percent[file+1][et];	
		}
		g_pTthresh_rat_v_eta[ctr] = new TGraph(n_eta_bins,avg_eta,threshold_rat_B[ctr]);
		prettyTGraph( g_pTthresh_rat_v_eta[ctr] , color[file] , 20+file , 0 , 4 , "#eta" , "p_{T} [MeV/#it{c}]" , "" );
		ctr++;
	}
	
	// -------------------------------------------------------------
        // Comparison to results from hybrid detector
        double hyb_eta    [] = {-3,-2.25,-1.75,-1.25,0,1.25,1.75,2.25,3};
        double hyb_pT_3_0T[] = {.15,.22,.16,.3,.4,.3,.16,.22,.15};
        double hyb_pT_1_4T[] = {.1,.13,.07,.15,.2,.15,.07,.13,.1};
        for(int i = 0 ; i < sizeof(hyb_eta)/sizeof(*hyb_eta) ; i++){
                hyb_pT_3_0T[i] *= 1000.; // Convert from GeV to MeV
                hyb_pT_1_4T[i] *= 1000.; // Convert from GeV to MeV
        }
        TGraph * g_hyb_pT_3_0T = new TGraph(9,hyb_eta,hyb_pT_3_0T);
        TGraph * g_hyb_pT_1_4T = new TGraph(9,hyb_eta,hyb_pT_1_4T);
        g_hyb_pT_3_0T -> SetLineColor(1);       g_hyb_pT_3_0T -> SetLineStyle(2);
        g_hyb_pT_1_4T -> SetLineColor(2);       g_hyb_pT_1_4T -> SetLineStyle(2);
	
	// -------------------------------------------------------------
	// Plotting histograms
	TCanvas * c0 = new TCanvas("c0","c0",1300,900);
        gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        h2_pT_eta_rat_finer[0] -> Draw("colz");
        c0 -> Modified();
        c0 -> Update();
        // ------------------------
	TCanvas * c1 = new TCanvas("c1","c1",1300,900);
	c1 -> Divide(3,size_fname);
	for(int file = 0 ; file < size_fname ; file++){
		c1 -> cd(1+3*file); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		h2_pT_eta_den_finer[file] -> Draw("colz");
		c1 -> cd(2+3*file); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		h2_pT_eta_num_finer[file] -> Draw("colz");
		c1 -> cd(3+3*file); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		h2_pT_eta_rat_finer[file] -> Draw("colz");
	}
	c1 -> Modified();
	c1 -> Update();
	// ------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1300,900);
	c2 -> Divide(3,size_fname);
	for(int file = 0 ; file < size_fname ; file++){
		c2 -> cd(1+3*file); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		h2_pT_eta_den_coarse[file] -> Draw("colz");
		c2 -> cd(2+3*file); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		h2_pT_eta_num_coarse[file] -> Draw("colz");
		c2 -> cd(3+3*file); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		h2_pT_eta_rat_coarse[file] -> Draw("colz");
	}
	c2 -> Modified();
	c2 -> Update();
	// ------------------------
	TCanvas * c3 = new TCanvas("c3","c3",1300,900);
	c3 -> Divide(size_fname/2,2);
	for(int file = 0 ; file < size_fname ; file++){
		c3 -> cd(file+1);
		h1_eff_pT[file][0] -> Draw();
		for(int et = 0 ; et < n_eta_bins ; et++){
			h1_eff_pT[file][et] -> Draw("same");
		}
	}
	// -----
	TLegend * leg1 = new TLegend(0.3,0.15,0.85,0.5);
	leg1 -> SetLineColor(0);
	leg1 -> SetNColumns(2);
	leg1 -> SetHeader("#eta bin","C");
	for(int et = 0 ; et < n_eta_bins ; et++){
		leg1 -> AddEntry( h1_eff_pT[0][et] , Form("(%.1f, %.1f)",eta_bin[et],eta_bin[et+1]) );
	}
	leg1 -> Draw("same");
	// -----
	c3 -> Modified();
	c3 -> Update();
	// ------------------------
	TCanvas * c4 = new TCanvas("c4","c4",1300,900);
	gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.03);
	g_pTthresh_v_eta[0] -> Draw("APL");
	for(int file = 0 ; file < size_fname ; file++){
		g_pTthresh_v_eta[file] -> Draw("samePL");
	}
	// -----
	TLegend * leg2 = new TLegend(0.7,0.7,0.95,0.95);
	leg2 -> SetLineColor(0);
	for(int file = 0 ; file < size_fname ; file++) leg2 -> AddEntry( g_pTthresh_v_eta[file] , B_field[file]);
	leg2 -> Draw("same");
	// -----
	g_hyb_pT_3_0T -> Draw("same");
	g_hyb_pT_1_4T -> Draw("same");
	TLatex * tex = new TLatex(3.07,111.86,"hybrid");
	tex -> SetTextSize(0.04);
	tex -> Draw("same");
	// -----
	c4 -> Modified();
	c4 -> Update();
	// ------------------------
        TCanvas * c5 = new TCanvas("c5","c5",1300,900);	// Plotting odd files only (which should correspond to 1.4 T)
        gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.03);
        g_pTthresh_v_eta[1] -> Draw("APL");
        for(int file = 1 ; file < size_fname ; file+=2){
                g_pTthresh_v_eta[file] -> Draw("samePL");
        }
        // -----
        TLegend * leg3 = new TLegend(0.7,0.7,0.95,0.95);
        leg3 -> SetLineColor(0);
        for(int file = 1 ; file < size_fname ; file+=2) leg3 -> AddEntry( g_pTthresh_v_eta[file] , B_field[file]);
        leg3 -> Draw("same");
        // -----
        // Comparison to results from hybrid detector
        g_hyb_pT_1_4T -> Draw("same");
        tex -> Draw("same");
        // -----
        c5 -> Modified();
        c5 -> Update();
	// ------------------------
        TCanvas * c6 = new TCanvas("c6","c6",1300,900); // Plotting even files only (which should correspond to 3.0 T)
        gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.03);
        g_pTthresh_v_eta[0] -> Draw("APL");
        for(int file = 0 ; file < size_fname ; file+=2){
                g_pTthresh_v_eta[file] -> Draw("samePL");
        }
        // -----
        TLegend * leg4 = new TLegend(0.7,0.7,0.95,0.95);
        leg4 -> SetLineColor(0);
        for(int file = 0 ; file < size_fname ; file+=2) leg4 -> AddEntry( g_pTthresh_v_eta[file] , B_field[file]);
        leg4 -> Draw("same");
        // -----
        // Comparison to results from hybrid detector
        g_hyb_pT_3_0T -> Draw("same");
        tex -> Draw("same");
        // -----
        c6 -> Modified();
        c6 -> Update();
	// ------------------------
        TCanvas * c7 = new TCanvas("c7","c7",1300,900); 
        gPad -> SetRightMargin(0.03); gPad -> SetTopMargin(0.03);
	g_pTthresh_rat_v_eta[0] -> Draw("AL");
	for(int i = 0 ; i < ctr ; i++){
		g_pTthresh_rat_v_eta[i] -> Draw("sameL");
	}
	// -----
        c7 -> Modified();
        c7 -> Update();

	// -------------------------------------------------------------
	// Saving results as pdfs
	TString res_fname = "results.pdf";
	c1 -> Print( res_fname + "(" );
	c2 -> Print( res_fname );
	c3 -> Print( res_fname );
	c4 -> Print( res_fname );
	c5 -> Print( res_fname );
	c6 -> Print( res_fname + ")" );

	// -------------------------------------------------------------
	myapp -> Run();
}
// ============================================================================================================================================
void prettyTH1D( TH1D * h1 , int color , int marker , float min , float max ){
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
// ============================================================================================================================================
void prettyTGraph( TGraph * g , int color , int marker , float min , float max , TString xtit , TString ytit , TString tit ){
	g -> SetLineWidth(2);
	g -> SetLineColor(color);
	g -> SetMarkerStyle(marker);
	g -> SetMarkerColor(color);

	g -> SetMinimum(min);
	g -> SetMaximum(max);

	g -> GetXaxis() -> CenterTitle();
	g -> GetXaxis() -> SetNdivisions(108); // to draw less tick marks
	g -> GetXaxis() -> SetTitle(xtit);
	g -> GetYaxis() -> CenterTitle();
	g -> GetYaxis() -> SetNdivisions(108); // to draw less tick marks
	g -> GetYaxis() -> SetTitle(ytit);

	g -> SetTitle(tit);

	g -> SetMinimum(0.001);
}
