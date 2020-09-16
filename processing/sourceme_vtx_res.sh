make clean
make

rm tables/tab_vtx_res_*.txt

# Loop over magnetic field setting
for B in {"Beast","sPHENIX"}
do
	# Loop over pixel size
	for pix in {20,10}
	do
		# Loop over detector configuration
		for det in {1,2,3}
		do
			# Do this three times to adjust tables
			for i in {1..3}
			do
				./analysis_vertex_resolution 1 1 "skimmed_combined_vtx_new_pi-_det"$det"_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
			done

			# Do one last time without updating table
			./analysis_vertex_resolution 1 1 "skimmed_combined_vtx_new_pi-_det"$det"_"$pix"x"$pix"_"$B"_FastTrackingEval.root"
		done
	done
done

say "All done!"
