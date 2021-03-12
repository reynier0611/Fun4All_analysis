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
TGraphErrors * graph_from_histo( TH1F * h1 , int color , int marker , float min , float max );
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetErrorX(0.0001);
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
		"../../output/output_mom_res_skimmed_pi-_det2_10x10_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",			//  0 - Beast, 10um, All-Si
		"../../output/output_mom_res_skimmed_pi-_det2_both_GEMs_RICH_10x10_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",       //  1 - Beast, 10um, All-Si + 50um res GEM
		"../../output/output_mom_res_skimmed_pi-_det2_10umGEM_RICH_10x10_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",         //  2 - Beast, 10um, All-Si + 10um Si disk
		"../../output/output_mom_res_skimmed_pi-_det2_LoResGEM_RICH_10x10_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",        //  3 - Beast, 10um, All-Si + low-res GEM

		"../../output/output_mom_res_skimmed_pi-_det2_20x20_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",                      //  4 - Beast, 20um, All-Si
		"../../output/output_mom_res_skimmed_pi-_det2_both_GEMs_RICH_20x20_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",       //  5 - Beast, 20um, All-Si + 50um res GEM
		"../../output/output_mom_res_skimmed_pi-_det2_10umGEM_RICH_20x20_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",         //  6 - Beast, 20um, All-Si + 10um Si disk
		"../../output/output_mom_res_skimmed_pi-_det2_LoResGEM_RICH_20x20_Beast_FastTrackingEvalsigma_eta_16_p_6_.root",        //  7 - Beast, 20um, All-Si + low-res GEM

		"../../output/output_mom_res_skimmed_pi-_det2_10x10_sPHENIX_FastTrackingEvalsigma_eta_16_p_6_.root",                    //  8 - BaBar, 10um, All-Si
		"../../output/output_mom_res_skimmed_pi-_det2_both_GEMs_RICH_10x10_sPHENIX_FastTrackingEvalsigma_eta_16_p_6_.root",     //  9 - BaBar, 10um, All-Si + 50um res GEM

		"../../output/output_mom_res_skimmed_pi-_det2_20x20_sPHENIX_FastTrackingEvalsigma_eta_16_p_6_.root",                    // 10 - BaBar, 20um, All-Si
		"../../output/output_mom_res_skimmed_pi-_det2_both_GEMs_RICH_20x20_sPHENIX_FastTrackingEvalsigma_eta_16_p_6_.root"      // 11 - BaBar, 20um, All-Si + 50um res GEM
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
	int color [] = { 1,51,94, 2, 1,51,94, 2, 1,51, 1,51};
	int marker[] = {21,22,23,24,21,22,23,24,21,22,21,22};
	TGraphErrors *** g_dpp_v_et_p_bins_BW = new TGraphErrors ** [size_loaded];
	TGraphErrors *** g_dpp_v_et_p_bins_FW = new TGraphErrors ** [size_loaded];
	for(int f = 0 ; f < size_loaded ; f++){
		g_dpp_v_et_p_bins_BW[f] = new TGraphErrors * [num_mom_bin[0]];
		g_dpp_v_et_p_bins_FW[f] = new TGraphErrors * [num_mom_bin[0]];
		for(int p = 0 ; p < num_mom_bin[f] ; p++){
			g_dpp_v_et_p_bins_BW[f][p] = graph_from_histo( h1_dpp_v_et_p_bins[f][p] ,color[f],marker[f],0.1,11);
			g_dpp_v_et_p_bins_FW[f][p] = graph_from_histo( h1_dpp_v_et_p_bins[f][p] ,color[f],marker[f],0.1,11);

			g_dpp_v_et_p_bins_BW[f][p] -> GetXaxis() -> SetRangeUser(-3.9,-1.6);
			g_dpp_v_et_p_bins_FW[f][p] -> GetXaxis() -> SetRangeUser( 1.6, 3.9);
		}
	}
	// ------------------------------------------------------------------------------
	// Plotting graphs
	double lowx[] = {0.00,0.38,0.69,0.00,0.38,0.69};
	double lowy[] = {0.55,0.55,0.55,0.00,0.00,0.00};
	double higx[] = {0.38,0.69,1.00,0.38,0.69,1.00};
	double higy[] = {1.00,1.00,1.00,0.55,0.55,0.55};
	// --------------------------
	TLegend * leg1 = new TLegend(0.25,0.47,0.99,0.86);
	leg1 -> SetLineColor(0);
	leg1 -> AddEntry(g_dpp_v_et_p_bins_BW[0][0],"All-Si");
	leg1 -> AddEntry(g_dpp_v_et_p_bins_BW[1][0],"All-Si + #sigma(50 #mum) GEM");
	leg1 -> AddEntry(g_dpp_v_et_p_bins_BW[2][0],"All-Si + 10 #mum pixel Si disk");
	//leg1 -> AddEntry(g_dpp_v_et_p_bins_BW[3][0],"All-Si + #sigma_{#hat{#phi}}(70 #mum) GEM");
	// --------------------------
        TLegend * leg2 = new TLegend(0.25,0.47,0.99,0.86);
        leg2 -> SetLineColor(0);
        leg2 -> AddEntry(g_dpp_v_et_p_bins_FW[0][0],"All-Si");
        leg2 -> AddEntry(g_dpp_v_et_p_bins_FW[1][0],"All-Si + #sigma(50 #mum) GEM");
        leg2 -> AddEntry(g_dpp_v_et_p_bins_FW[2][0],"All-Si + 10 #mum pixel Si disk");
	//leg2 -> AddEntry(g_dpp_v_et_p_bins_FW[3][0],"All-Si + #sigma_{#hat{#phi}}(70 #mum) GEM");
	// --------------------------
	TLatex ** tex_tit_back = new TLatex * [num_mom_bin[0]];
	TLatex ** tex_tit_forw = new TLatex * [num_mom_bin[0]];
	for(int p  = 0 ; p  < num_mom_bin[0] ; p ++){
		tex_tit_back[p] = new TLatex(-3.4,10,Form("#bf{%.0f < p < %.0f GeV/#it{c}}",(*TVT_mom_bin[0])[p],(*TVT_mom_bin[0])[p+1]));
		tex_tit_forw[p] = new TLatex( 2.0,10,Form("#bf{%.0f < p < %.0f GeV/#it{c}}",(*TVT_mom_bin[0])[p],(*TVT_mom_bin[0])[p+1]));
		if(p < 3){
			tex_tit_back[p] -> SetTextSize(0.09);
			tex_tit_forw[p] -> SetTextSize(0.09);
		}
		else{
			tex_tit_back[p] -> SetTextSize(0.08);
			tex_tit_forw[p] -> SetTextSize(0.08);
		}
	}
	TLatex * tex_z_back = new TLatex(-3.,3.85,"#bf{#splitline{auxiliary tracking station}{at z = -180 cm}}");
	tex_z_back -> SetTextSize(0.07);
	TLatex * tex_z_forw = new TLatex( 2.,3.85,"#bf{#splitline{auxiliary tracking station}{at z = 300 cm}}");
        tex_z_forw -> SetTextSize(0.07);
	// --------------------------
	TCanvas * c1 = new TCanvas("c1","c1",1200,700);
	c1 -> Divide(3,2);
	TVirtualPad ** pad1 = new TVirtualPad * [num_mom_bin[0]];
	for(int p  = 0 ; p  < num_mom_bin[0] ; p ++){
		c1 -> cd(p+1);
		pad1[p] = c1 -> cd(p+1);
		pad1[p] -> SetPad(Form("pad1_%i",p),Form("pad1_%i",p),lowx[p],lowy[p],higx[p],higy[p],kWhite, 0, 0);

		if(p!=2&&p!=5) pad1[p] -> SetRightMargin(0.0);
		else pad1[p] -> SetRightMargin(0.02);
		if(p>=3) pad1[p] -> SetTopMargin(0);
		else pad1[p] -> SetTopMargin(0.02);
		if(p<3) pad1[p] -> SetBottomMargin(0.0);
		else pad1[p] -> SetBottomMargin(0.18);
		if(p!=0&&p!=3) pad1[p] -> SetLeftMargin(0.0);
		else pad1[p] -> SetLeftMargin(0.20);

		if(p==0){
			g_dpp_v_et_p_bins_BW[0][p]->GetYaxis()->SetTitleSize(0.09);
        		g_dpp_v_et_p_bins_BW[0][p]->GetYaxis()->SetLabelSize(0.09);
		}

		g_dpp_v_et_p_bins_BW[0][p] -> Draw("APL");	// All-Si alone
		g_dpp_v_et_p_bins_BW[1][p] -> Draw("samePL");	// 50um-res GEM
		g_dpp_v_et_p_bins_BW[2][p] -> Draw("samePL");	// Si-disk	
		//g_dpp_v_et_p_bins_BW[3][p] -> Draw("samePL");   // Low-res GEM

		tex_tit_back[p] -> Draw("same");
	}
	c1 -> cd(1);
	leg1 -> Draw("same");
	tex_z_back -> Draw("same");
	c1 -> Modified();
	c1 -> Update();
	// --------------------------
	TCanvas * c2 = new TCanvas("c2","c2",1200,700);
        c2 -> Divide(3,2);
        TVirtualPad ** pad2 = new TVirtualPad * [num_mom_bin[0]];
        for(int p  = 0 ; p  < num_mom_bin[0] ; p ++){
                c2 -> cd(p+1);
                pad2[p] = c2 -> cd(p+1);
                pad2[p] -> SetPad(Form("pad2_%i",p),Form("pad2_%i",p),lowx[p],lowy[p],higx[p],higy[p],kWhite, 0, 0);

                if(p!=2&&p!=5) pad2[p] -> SetRightMargin(0.0);
                else pad2[p] -> SetRightMargin(0.02);
                if(p>=3) pad2[p] -> SetTopMargin(0);
                else pad2[p] -> SetTopMargin(0.02);
                if(p<3) pad2[p] -> SetBottomMargin(0.0);
                else pad2[p] -> SetBottomMargin(0.18);
                if(p!=0&&p!=3) pad2[p] -> SetLeftMargin(0.0);
                else pad2[p] -> SetLeftMargin(0.20);

                if(p==0){
                        g_dpp_v_et_p_bins_FW[0][p]->GetYaxis()->SetTitleSize(0.09);
                        g_dpp_v_et_p_bins_FW[0][p]->GetYaxis()->SetLabelSize(0.09);
                }

                g_dpp_v_et_p_bins_FW[0][p] -> Draw("APL");      // All-Si alone
                g_dpp_v_et_p_bins_FW[1][p] -> Draw("samePL");   // 50um-res GEM
                g_dpp_v_et_p_bins_FW[2][p] -> Draw("samePL");   // Si-disk
                //g_dpp_v_et_p_bins_FW[3][p] -> Draw("samePL"); // Low-res GEM

                tex_tit_forw[p] -> Draw("same");
        }
        c2 -> cd(1);
        leg2 -> Draw("same");
        tex_z_forw -> Draw("same");
        c2 -> Modified();
        c2 -> Update();
	// --------------------------
	c1 -> Print("results_plots_for_YR.pdf(");
	c2 -> Print("results_plots_for_YR.pdf)");

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
	g_l1->SetMarkerSize(1.4);
	g_l1->SetLineWidth(2);

	g_l1->GetXaxis()->SetTitle(xax);
	g_l1->GetXaxis()->SetNdivisions(108);
	g_l1->GetXaxis()->SetTitleSize(0.08);
	g_l1->GetXaxis()->SetLabelSize(0.08);
	g_l1->GetXaxis()->CenterTitle();

	g_l1->GetYaxis()->SetTitle(yax);
	g_l1->GetYaxis()->SetNdivisions(108);
	g_l1->GetYaxis()->SetTitleSize(0.08);
	g_l1->GetYaxis()->SetLabelSize(0.08);
	g_l1->GetYaxis()->CenterTitle();

	g_l1->SetTitle(tit);

	if(min!=999) g_l1 -> SetMinimum(min);
	if(max!=999) g_l1 -> SetMaximum(max);

	float minxval = h1 -> GetBinLowEdge(1);
	float maxxval = h1 -> GetBinLowEdge(nbin) + h1 -> GetBinWidth(nbin); 

	g_l1->GetXaxis()->SetRangeUser( minxval , maxxval );

	return g_l1;
}
