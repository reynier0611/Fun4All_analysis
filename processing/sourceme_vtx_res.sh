make clean
make

# Loop over magnetic field setting
for B in {"Beast","sPHENIX"}
do
	# Do this three times to adjust tables
	for i in {1..3}
	do
		./analysis_vertex_resolution 1 1 1 "skimmed_combined_vtx_new_pi-_det2_10x10_"$B"_FastTrackingEval.root"
	done

	# Do one last time without updating table
	./analysis_vertex_resolution 1 2 1 "skimmed_combined_vtx_new_pi-_det2_10x10_"$B"_FastTrackingEval.root"

done

say "All done!"
