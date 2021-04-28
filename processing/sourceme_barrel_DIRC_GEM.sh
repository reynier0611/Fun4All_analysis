cd ../data/
#cp barrel_DIRC_GEM/*.root .
cp barrel_DIRC_GEM/skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_barrelR_75cm_Beast_FastSimEval.root .
cp barrel_DIRC_GEM/skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_barrelR_75cm_sPHENIX_FastSimEval.root .
cd ../processing/
make clean
make

#for file in "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_Beast_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_49cm_GEM_R_60cm_50um_Beast_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_49cm_GEM_R_60cm_50um_sPHENIX_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_84cm_GEM_R_92cm_50um_Beast_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_DIRC_R_84cm_GEM_R_92cm_50um_sPHENIX_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_60cm_50um_Beast_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_60cm_50um_sPHENIX_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_92cm_50um_Beast_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_R_92cm_50um_sPHENIX_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_sPHENIX_FastSimEval.root"
for file in "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_barrelR_75cm_Beast_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_barrelR_75cm_sPHENIX_FastSimEval.root"
do
	# Do this three times to adjust tables
	for i in {1..3}
	do
		./analysis_momentum_resolution 1 1 1 $file "barrel_DIRC_GEM.txt"
	done

	# Do one last time without updating table
	./analysis_momentum_resolution 1 2 1 $file "barrel_DIRC_GEM.txt"

	cd ../data/
	rm $file
	cd ../processing/
done

say "Finished running"
