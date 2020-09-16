make clean
make

rm tables/tab_mom_res_*.txt

# Loop over magnetic field setting
for B in {"Beast","sPHENIX"}
do
	# Loop over pixel size
	for pix in {20,10}
	do
		# Loop over detector configuration
		for det in {"","_both_GEMs","_both_GEMs_RICH"}
		do
			# Do this three times to adjust tables
			for i in {1..3}
			do
				./analysis_resolution 1 1 "skimmed_pi-_det2"$det"_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
			done

			# Do one last time without updating table
			./analysis_resolution 1 2 "skimmed_pi-_det2"$det"_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
		done
	done
done

say "All done!"
