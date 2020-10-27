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
#include "TLine.h"

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max );
double sq(double x){return{x*x};}
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	TString pix_size = "10";
	gStyle->SetErrorX(0.0001);
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames_vtx[] = {
		"../../output/output_vtx_res_skimmed_combined_vtx_new_pi-_det2_"+pix_size+"x"+pix_size+"_Beast_FastTrackingEvalsigma_eta_4_p_16_.root"
		//"../../output/output_vtx_res_skimmed_combined_vtx_new_pi-_det2_10x10_sPHENIX_FastTrackingEval.root"
	};
	TString fnames_mom[] = {	
		"../../output/output_skimmed_pi-_det2_"+pix_size+"x"+pix_size+"_Beast_FastTrackingEval.root"   ,    //  0 
                "../../output/output_skimmed_pi-_det2_"+pix_size+"x"+pix_size+"_sPHENIX_FastTrackingEval.root"      //  1
	};
	// #######################################################################################################################################
	// YOU SHOULDN'T NEED TO MODIFY ANYTHING IN THE BLOCK OF CODE BELOW AND UNTIL AFTER THE NEXT LINE WITH ###...
	const int size_loaded_vtx = sizeof(fnames_vtx)/sizeof(*fnames_vtx);
	const int size_loaded_mom = sizeof(fnames_mom)/sizeof(*fnames_mom);
	// ------------------------------------------------------------------------------
	// Preparing variables that will later on be filled from root files
	TVectorT<double> ** TVT_eta_bin_vtx = new TVectorT<double>*[size_loaded_vtx];
	TVectorT<double> ** TVT_mom_bin_vtx = new TVectorT<double>*[size_loaded_vtx];
	TVectorT<double> ** TVT_eta_bin_mom = new TVectorT<double>*[size_loaded_mom];
	TVectorT<double> ** TVT_mom_bin_mom = new TVectorT<double>*[size_loaded_mom];
	int num_eta_bin_vtx[size_loaded_vtx] = {0};
	int num_mom_bin_vtx[size_loaded_vtx] = {0};
	int num_eta_bin_mom[size_loaded_mom] = {0};
	int num_mom_bin_mom[size_loaded_mom] = {0};
	// ------------------------------------------------------------------------------
	// Loading root files and info therein
	TFile ** Fin_vtx = new TFile*[size_loaded_vtx];
	TFile ** Fin_mom = new TFile*[size_loaded_mom];

	TH1F *** h1_dvl_v_p_et_bins = new TH1F ** [size_loaded_vtx];
	TH1F *** h1_dvt_v_p_et_bins = new TH1F ** [size_loaded_vtx];
	TH1F *** h1_dvl_v_et_p_bins = new TH1F ** [size_loaded_vtx];
	TH1F *** h1_dvt_v_et_p_bins = new TH1F ** [size_loaded_vtx];

	TH1F *** h1_dpp_v_p_et_bins = new TH1F ** [size_loaded_mom]; 
        TH1F *** h1_dpp_v_et_p_bins = new TH1F ** [size_loaded_mom];

	for(int f = 0 ; f < size_loaded_vtx ; f++){
		ifstream fin;
		fin.open(fnames_vtx[f]);
		if(!fin){ cout << "\033[1;31mCouldn't find input file '" << fnames_vtx[f] << "'. Bailing out!\033[0m" << endl; exit(0);}
		fin.close();

		Fin_vtx[f] = new TFile(fnames_vtx[f]);

		TVT_eta_bin_vtx[f] = (TVectorT<double> *) Fin_vtx[f] -> Get("TVT_eta_bin");
		TVT_mom_bin_vtx[f] = (TVectorT<double> *) Fin_vtx[f] -> Get("TVT_mom_bin");

		num_eta_bin_vtx[f] = (*TVT_eta_bin_vtx[f]).GetNoElements()-1;
		num_mom_bin_vtx[f] = (*TVT_mom_bin_vtx[f]).GetNoElements()-1;

		h1_dvl_v_p_et_bins[f] = new TH1F * [num_eta_bin_vtx[f]];
		h1_dvt_v_p_et_bins[f] = new TH1F * [num_eta_bin_vtx[f]];
		h1_dvl_v_et_p_bins[f] = new TH1F * [num_mom_bin_vtx[f]];
		h1_dvt_v_et_p_bins[f] = new TH1F * [num_mom_bin_vtx[f]];

		for(int et = 0 ; et < num_eta_bin_vtx[f] ; et++){
			h1_dvl_v_p_et_bins[f][et] = (TH1F*) Fin_vtx[f] -> Get(Form("h1_dvl_v_pT_et_bins_%i",et));
			h1_dvt_v_p_et_bins[f][et] = (TH1F*) Fin_vtx[f] -> Get(Form("h1_dvt_v_pT_et_bins_%i",et));
		}
		for(int p = 0 ; p < num_mom_bin_vtx[f] ; p++){
			h1_dvl_v_et_p_bins[f][p ] = (TH1F*) Fin_vtx[f] -> Get(Form("h1_dvl_v_et_p_bins_%i",p ));
			h1_dvt_v_et_p_bins[f][p ] = (TH1F*) Fin_vtx[f] -> Get(Form("h1_dvt_v_et_p_bins_%i",p ));
		}
	}
	cout << "\nvertex resolution file:\n";
	cout << "\neta bin boundaries:\n"; for(int et = 0 ; et < num_eta_bin_vtx[0]+1 ; et++) cout << (*TVT_eta_bin_vtx[0])[et] << ", "; cout << "\n";
	cout << "\np bin boundaries:\n"  ; for(int p  = 0 ; p  < num_mom_bin_vtx[0]+1 ; p ++) cout << (*TVT_mom_bin_vtx[0])[ p] << ", "; cout << "\n\n";

	for(int f = 0 ; f < size_loaded_mom ; f++){
		
                ifstream fin;
                fin.open(fnames_mom[f]);
                if(!fin){ cout << "\033[1;31mCouldn't find input file '" << fnames_mom[f] << "'. Bailing out!\033[0m" << endl; exit(0);}
                fin.close();

                Fin_mom[f] = new TFile(fnames_mom[f]);

                TVT_eta_bin_mom[f] = (TVectorT<double> *) Fin_mom[f] -> Get("TVT_eta_bin");
                TVT_mom_bin_mom[f] = (TVectorT<double> *) Fin_mom[f] -> Get("TVT_mom_bin");

                num_eta_bin_mom[f] = (*TVT_eta_bin_mom[f]).GetNoElements()-1;
                num_mom_bin_mom[f] = (*TVT_mom_bin_mom[f]).GetNoElements()-1;

                h1_dpp_v_p_et_bins[f] = new TH1F * [num_eta_bin_mom[f]];
                h1_dpp_v_et_p_bins[f] = new TH1F * [num_mom_bin_mom[f]];
		
                for(int et = 0 ; et < num_eta_bin_mom[f] ; et++){
                        h1_dpp_v_p_et_bins[f][et] = (TH1F*) Fin_mom[f] -> Get(Form("h1_dpp_v_p_et_bins_%i",et));
                }
                for(int p = 0 ; p < num_mom_bin_mom[f] ; p++){
                        h1_dpp_v_et_p_bins[f][p ] = (TH1F*) Fin_mom[f] -> Get(Form("h1_dpp_v_et_p_bins_%i",p ));
                }
	}
	cout << "\nmomentum resolution file:\n";
        cout << "\neta bin boundaries:\n"; for(int et = 0 ; et < num_eta_bin_mom[0]+1 ; et++) cout << (*TVT_eta_bin_mom[0])[et] << ", "; cout << "\n";
        cout << "\np bin boundaries:\n"  ; for(int p  = 0 ; p  < num_mom_bin_mom[0]+1 ; p ++) cout << (*TVT_mom_bin_mom[0])[ p] << ", "; cout << "\n\n";	

	// #######################################################################################################################################
        // EDIT THE CODE BELOW DEPENDING ON WHAT YOU WANT TO PLOT
	// ------------------------------------------------------------------------------
	double vtx_bin[] = {0.,1.,2.};
	double mom_bin[] = {1.,25.};
        const int eta_0_vtx = idx_from_vector(0.,TVT_eta_bin_vtx[0]);
        const int eta_1_vtx = idx_from_vector(1.,TVT_eta_bin_vtx[0]);
	const int eta_2_vtx = idx_from_vector(2.,TVT_eta_bin_vtx[0]);
        const int p_1GeV_mom  = idx_from_vector( 1.,TVT_mom_bin_mom[0]);
	const int p_25GeV_mom = idx_from_vector(25.,TVT_mom_bin_mom[0]);
	int selected_bins_vtx[] = {eta_0_vtx,eta_1_vtx,eta_2_vtx};
	int selected_bins_mom[] = {p_1GeV_mom,p_25GeV_mom};
        int size_selected_bins_vtx = sizeof(selected_bins_vtx)/sizeof(*selected_bins_vtx);
	int size_selected_bins_mom = sizeof(selected_bins_mom)/sizeof(*selected_bins_mom);

	TGraphErrors ** g_dvl_v_p_Beast = new TGraphErrors * [size_selected_bins_vtx];
        TGraphErrors ** g_dvl_v_p_Babar = new TGraphErrors * [size_selected_bins_vtx];
	
	TGraphErrors ** g_dvt_v_p_Beast = new TGraphErrors * [size_selected_bins_vtx];
        TGraphErrors ** g_dvt_v_p_Babar = new TGraphErrors * [size_selected_bins_vtx];

	TGraphErrors ** g_dpp_v_eta_Beast = new TGraphErrors * [size_selected_bins_mom];
        TGraphErrors ** g_dpp_v_eta_Babar = new TGraphErrors * [size_selected_bins_mom];

	int vtx_color[] = {62,2,8};
	int vtx_marker[] = {20,21,22};
	
    	for(int i = 0 ; i < size_selected_bins_vtx ; i++){	
		g_dvl_v_p_Beast[i] = graph_from_histo(h1_dvl_v_p_et_bins[0][selected_bins_vtx[i]], vtx_color[i],vtx_marker[i],4,100);	g_dvl_v_p_Beast[i] -> GetYaxis() -> SetMoreLogLabels();
		//g_dvl_v_p_Babar[i] = graph_from_histo(h1_dvl_v_p_et_bins[1][selected_bins_vtx[i]], vtx_color[i],21,2,100);	g_dvl_v_p_Babar[i] -> GetYaxis() -> SetMoreLogLabels();
	
		g_dvt_v_p_Beast[i] = graph_from_histo(h1_dvt_v_p_et_bins[0][selected_bins_vtx[i]], vtx_color[i],vtx_marker[i],4,100);  g_dvt_v_p_Beast[i] -> GetYaxis() -> SetMoreLogLabels();
                //g_dvt_v_p_Babar[i] = graph_from_histo(h1_dvt_v_p_et_bins[1][selected_bins_vtx[i]],  2,21,2,100);  g_dvt_v_p_Babar[i] -> GetYaxis() -> SetMoreLogLabels();

		//g_dvl_v_p_Beast[i] -> SetTitle(Form("3.0 T, %.1f < |#eta| < %.1f, "+mat_bud+"%% X/X_{0}, "+pix_size+"#mum pixel",(*TVT_eta_bin_vtx[0])[selected_bins_vtx[i]],(*TVT_eta_bin_vtx[0])[selected_bins_vtx[i]+1])); 

		g_dvl_v_p_Beast[i] -> GetXaxis() -> SetRangeUser(0,7);
		g_dvt_v_p_Beast[i] -> GetXaxis() -> SetRangeUser(0,7);
		g_dvl_v_p_Beast[i] -> SetMaximum(1000);
		g_dvt_v_p_Beast[i] -> SetMaximum(1000);
	}
	
	int mom_marker[] = {24,20,25,21};

	for(int i = 0 ; i < size_selected_bins_mom ; i++){      
         	g_dpp_v_eta_Beast[i] = graph_from_histo(h1_dpp_v_et_p_bins[0][selected_bins_mom[i]],  2,mom_marker[i  ],.3,.7);        g_dpp_v_eta_Beast[i] -> GetYaxis() -> SetMoreLogLabels();
                g_dpp_v_eta_Babar[i] = graph_from_histo(h1_dpp_v_et_p_bins[1][selected_bins_mom[i]], 62,mom_marker[i+2],.3,.9);        g_dpp_v_eta_Babar[i] -> GetYaxis() -> SetMoreLogLabels();

		g_dpp_v_eta_Beast[i] -> SetMaximum(12); 
		g_dpp_v_eta_Beast[i] -> GetXaxis() -> SetRangeUser(0,3.7);
        }
	
	// ------------------------------------------------------------------------------
	TLegend * leg1 = new TLegend(0.20,0.4,0.75,0.95);
        leg1 -> SetLineColor(0);
        leg1 -> AddEntry((TObject*)0, "B = 3.0 T", "B = 3.0 T");
	leg1 -> AddEntry( g_dpp_v_eta_Beast[0] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin_mom[0])[selected_bins_mom[0]],(*TVT_mom_bin_mom[0])[selected_bins_mom[0]+1]) );
	leg1 -> AddEntry( g_dpp_v_eta_Beast[1] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin_mom[0])[selected_bins_mom[1]],(*TVT_mom_bin_mom[1])[selected_bins_mom[1]+1]) );
        leg1 -> AddEntry( (TObject*)0, "" , "");
	leg1 -> AddEntry( (TObject*)0, "B = 1.4 T", "B = 1.4 T");
	leg1 -> AddEntry( g_dpp_v_eta_Babar[0] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin_mom[0])[selected_bins_mom[0]],(*TVT_mom_bin_mom[0])[selected_bins_mom[0]+1]) );
	leg1 -> AddEntry( g_dpp_v_eta_Babar[1] , Form("%.0f < p < %.0f GeV/#it{c}",(*TVT_mom_bin_mom[0])[selected_bins_mom[1]],(*TVT_mom_bin_mom[1])[selected_bins_mom[1]+1]) );
	// ------------
	TLegend * leg2 = new TLegend(0.65,0.5,0.95,0.95);
        leg2 -> SetLineColor(0);
	leg2 -> AddEntry((TObject*)0, "B = 3.0 T", "B = 3.0 T");
	leg2 -> AddEntry( g_dvl_v_p_Beast[0] , Form("%.0f < #eta < %.0f",(*TVT_eta_bin_vtx[0])[selected_bins_vtx[0]],(*TVT_eta_bin_vtx[0])[selected_bins_vtx[0]+1]) );
	leg2 -> AddEntry( g_dvl_v_p_Beast[1] , Form("%.0f < #eta < %.0f",(*TVT_eta_bin_vtx[0])[selected_bins_vtx[1]],(*TVT_eta_bin_vtx[0])[selected_bins_vtx[1]+1]) );
	leg2 -> AddEntry( g_dvl_v_p_Beast[2] , Form("%.0f < #eta < %.0f",(*TVT_eta_bin_vtx[0])[selected_bins_vtx[2]],(*TVT_eta_bin_vtx[0])[selected_bins_vtx[2]+1]) );
	// ------------
	// ------------------------------------------------------------------------------
	// Plotting graphs
	TCanvas * c1 = new TCanvas("c1","c1",1400,400);
	c1 -> Divide(3,1);
	// ------------
        c1 -> cd(1);
        gPad -> SetRightMargin(0.01); gPad -> SetBottomMargin(0.16); gPad -> SetLeftMargin(0.19); gPad -> SetTopMargin(0.03);
	g_dpp_v_eta_Beast[0] -> Draw(   "APL");
        g_dpp_v_eta_Beast[1] -> Draw("samePL");
	g_dpp_v_eta_Babar[0] -> Draw("samePL");
	g_dpp_v_eta_Babar[1] -> Draw("samePL");
	leg1 -> Draw("same");
	// ------------
	c1 -> cd(2);
	gPad -> SetRightMargin(0.01); gPad -> SetBottomMargin(0.16); gPad -> SetLeftMargin(0.19); gPad -> SetLogy(); gPad -> SetTopMargin(0.03);
	g_dvl_v_p_Beast[0] -> Draw(   "APL");
	g_dvl_v_p_Beast[1] -> Draw("samePL");
	g_dvl_v_p_Beast[2] -> Draw("samePL");
	//g_dvl_v_p_Babar[0] -> Draw("samePL");
	// ------------
	c1 -> cd(3);
        gPad -> SetRightMargin(0.01); gPad -> SetBottomMargin(0.16); gPad -> SetLeftMargin(0.19); gPad -> SetLogy(); gPad -> SetTopMargin(0.03);
        g_dvt_v_p_Beast[0] -> Draw(   "APL");
	g_dvt_v_p_Beast[1] -> Draw("samePL");
	g_dvt_v_p_Beast[2] -> Draw("samePL");
        //g_dvt_v_p_Babar[0] -> Draw("samePL");
	leg2 -> Draw("same");
	// ------------------------------------------------------------------------------
	// Parametrizations
	TF1 * f_LBNL_0_lt_eta_lt_1 = new TF1("f_LBNL_0_lt_eta_lt_1","sqrt(sq(20./x)+sq(5.0))",0,7);	f_LBNL_0_lt_eta_lt_1 -> SetLineColor(vtx_color[0]);	f_LBNL_0_lt_eta_lt_1 -> SetLineStyle(3);
	TF1 * f_LBNL_1_lt_eta_lt_2 = new TF1("f_LBNL_1_lt_eta_lt_2","sqrt(sq(25./x)+sq(10.))",0,7);	f_LBNL_1_lt_eta_lt_2 -> SetLineColor(vtx_color[1]);	f_LBNL_1_lt_eta_lt_2 -> SetLineStyle(3);
	TF1 * f_LBNL_2_lt_eta_lt_3 = new TF1("f_LBNL_2_lt_eta_lt_3","sqrt(sq(30./x)+sq(10.))",0,7);	f_LBNL_2_lt_eta_lt_3 -> SetLineColor(vtx_color[2]);	f_LBNL_2_lt_eta_lt_3 -> SetLineStyle(3);
	TF1 * f_LANL_0_lt_eta_lt_1 = new TF1("f_LANL_0_lt_eta_lt_1","sqrt(sq(25./x)+sq(0.0))",0,7);	f_LANL_0_lt_eta_lt_1 -> SetLineColor(vtx_color[0]);	f_LANL_0_lt_eta_lt_1 -> SetLineStyle(7);
	TF1 * f_LANL_1_lt_eta_lt_2 = new TF1("f_LANL_1_lt_eta_lt_2","sqrt(sq(30./x)+sq(20.))",0,7);	f_LANL_1_lt_eta_lt_2 -> SetLineColor(vtx_color[1]);	f_LANL_1_lt_eta_lt_2 -> SetLineStyle(7);
	TF1 * f_LANL_2_lt_eta_lt_3 = new TF1("f_LANL_2_lt_eta_lt_3","sqrt(sq(30./x)+sq(40.))",0,7);	f_LANL_2_lt_eta_lt_3 -> SetLineColor(vtx_color[2]);	f_LANL_2_lt_eta_lt_3 -> SetLineStyle(7);
	/*
	f_LBNL_0_lt_eta_lt_1 -> Draw("same");
	f_LBNL_1_lt_eta_lt_2 -> Draw("same");
	f_LBNL_2_lt_eta_lt_3 -> Draw("same");
	f_LANL_0_lt_eta_lt_1 -> Draw("same");
	f_LANL_1_lt_eta_lt_2 -> Draw("same");
	f_LANL_2_lt_eta_lt_3 -> Draw("same");
	*/
	// ------------------------------------------------------------------------------
	c1 -> Modified();
	c1 -> Update();
	// ------------------------------------------------------------------------------
	// Saving results to pdf files
	c1 -> Print("results_plotter.pdf");
	
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
			err[ctr] = h1 -> GetBinError(i+1);
			ctr++;
		}
	}

	// Creating a TGraphErrors object that mirrors the histogram
	TGraphErrors * g_l1 = new TGraphErrors(ctr,xval,yval,0,err);
	g_l1->SetLineColor(color);
	g_l1->SetMarkerColor(color);
	g_l1->SetMarkerStyle(marker);
	g_l1->SetMarkerSize(1.1);
	g_l1->SetLineWidth(2);

	g_l1->GetXaxis()->SetTitle(xax);
	g_l1->GetXaxis()->SetNdivisions(108);
	g_l1->GetXaxis()->SetTitleSize(0.06);
	g_l1->GetXaxis()->SetLabelSize(0.06);
	g_l1->GetXaxis()->SetTitleOffset(1.20);
	g_l1->GetXaxis()->CenterTitle();

	g_l1->GetYaxis()->SetTitle(yax);
	g_l1->GetYaxis()->SetNdivisions(108);
	g_l1->GetYaxis()->SetTitleSize(0.06);
	g_l1->GetYaxis()->SetLabelSize(0.06);
	g_l1->GetYaxis()->CenterTitle();
	g_l1->GetYaxis()->SetTitleOffset(1.55);
	g_l1->SetTitle(tit);

	if(min!=999) g_l1 -> SetMinimum(min);
        if(max!=999) g_l1 -> SetMaximum(max);

	float minxval = h1 -> GetBinLowEdge(1);
	float maxxval = h1 -> GetBinLowEdge(nbin) + h1 -> GetBinWidth(nbin); 

	g_l1->GetXaxis()->SetRangeUser( minxval , maxxval );

	return g_l1;
}
