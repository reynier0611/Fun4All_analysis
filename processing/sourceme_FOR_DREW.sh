cd ../data/
cp FOR_DREW/*.root .
cd ../processing/
make clean
make

for file in "out_simp_geom_vbd_0.05_0.55_0.24_Beast_FastSimEval.root" "out_simp_geom_vbd_0.05_0.55_0.24_sPHENIX_FastSimEval.root" "out_simp_geom_vbd_0.3_0.3_0.3_Beast_FastSimEval.root" "out_simp_geom_vbd_0.3_0.3_0.3_sPHENIX_FastSimEval.root"
do
	# Do this three times to adjust tables
	for i in {1..3}
	do
		./analysis_momentum_resolution 1 1 1 $file "positive_eta.txt"
	done

	# Do one last time without updating table
	./analysis_momentum_resolution 1 2 1 $file "positive_eta.txt"

	cd ../data/
	rm $file
	cd ../processing/
done

say "Finished running"
