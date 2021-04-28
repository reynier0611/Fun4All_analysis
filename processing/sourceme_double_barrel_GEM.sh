cd ../data/
cp barrel_DIRC_GEM/new/*.root .
cd ../processing/
make clean
make

for file in "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_DIRC_GEM_Beast_FastSimEval.root" "skimmed_out_simp_geom_vbd_0.05_0.55_0.24_10um_pix_GEM_DIRC_GEM_sPHENIX_FastSimEval.root"
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
