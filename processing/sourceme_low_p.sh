make clean
make

# Loop over magnetic field setting
for B in {"Beast","sPHENIX"}
do
	# Loop over pixel size
	for pix in {20,10}
	do
		# Do this three times to adjust tables
		for i in {1..3}
		do
			#./analysis_momentum_resolution 1 1 1 "skimmed_combined_p_0_1_GeV_pi-_det2_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
			./analysis_momentum_resolution 1 1 1 "new_pi-_det2_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
		done

		# Do one last time without updating table
		#./analysis_momentum_resolution 1 2 1 "skimmed_combined_p_0_1_GeV_pi-_det2_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
		./analysis_momentum_resolution 1 2 1 "new_pi-_det2_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
	done
done

say "All done!"
