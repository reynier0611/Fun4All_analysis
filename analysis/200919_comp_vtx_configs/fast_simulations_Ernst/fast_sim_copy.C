{

	double x[] = {0.3,0.4,0.6,0.8,1.,2.,5.,10.};
	double y_111[] = {44.1379,33.6990,24.1065,19.2163,16.3949,9.62382,4.63949,3.04075};
	double y_101[] = {41.7868,31.1598,21.4733,16.9592,14.1379,8.11912,4.07523,3.13479};
	double y_110[] = {44.7021,34.0752,23.7304,18.5579,15.1724,8.21316,4.07523,3.04075};

	TGraph * g_111 = new TGraph(8,x,y_111);
	TGraph * g_101 = new TGraph(8,x,y_101);
	TGraph * g_110 = new TGraph(8,x,y_110);

	g_111 -> GetXaxis() -> SetTitle("Momentum [GeV/#it{c}]");
	g_111 -> GetXaxis() -> SetNdivisions(107);
	g_111 -> GetXaxis() -> CenterTitle();
	g_111 -> GetXaxis() -> SetLabelSize(0.05);
	g_111 -> GetXaxis() -> SetTitleSize(0.05);
	g_111 -> GetYaxis() -> SetTitle("Longitudinal Pointing Resolution [#mum]");
	g_111 -> GetYaxis() -> SetNdivisions(107);
	g_111 -> GetYaxis() -> CenterTitle();
	g_111 -> GetYaxis() -> SetLabelSize(0.05);
	g_111 -> GetYaxis() -> SetTitleSize(0.05);
	g_111 -> SetTitle("3.0T, #eta = 0.5, 10#mum pixel, vtx #rightarrow 0.05% X/X_{0}");

	g_111 -> SetMarkerStyle(20);
	g_101 -> SetMarkerStyle(21);
	g_110 -> SetMarkerStyle(25);

	g_111 -> SetMarkerColor( 1); g_111 -> SetLineColor( 1);
	g_101 -> SetMarkerColor( 2); g_101 -> SetLineColor( 2);
	g_110 -> SetMarkerColor( 4); g_110 -> SetLineColor( 4);

	g_111 -> GetXaxis() -> SetRangeUser(0,8);

	TCanvas * c1 = new TCanvas("c1");
	gPad -> SetLeftMargin(0.12); gPad -> SetBottomMargin(0.12); gPad -> SetRightMargin(0.03);
	g_111 -> Draw("APL");
	g_101 -> Draw("samePL");
	g_110 -> Draw("samePL");

	TLegend * leg = new TLegend(0.7,0.6,0.85,0.85);
	leg -> SetLineColor(0);
	leg -> AddEntry(g_111,"1,1,1");
	leg -> AddEntry(g_110,"1,1,0");
	leg -> AddEntry(g_101,"1,0,1");
	leg -> Draw("same");

	TLatex * tex1 = new TLatex(1,43,"#bf{fast simulations}");	tex1 -> Draw("same");
	TLatex * tex2 = new TLatex(1,40,"#bf{E. Sichtermann}"  );	tex2 -> Draw("same");

	c1 -> Print("results_fast_sim_copy.pdf");
}
