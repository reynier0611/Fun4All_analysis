make clean
make

# Loop over detector configuration
for det in {"111","110","101","011","100","010","001"}
do
	for BF in {"1.4","3.1"}
	do
		# Do this three times to adjust tables
		for i in {1..3}
		do
			./copy_momentum 1 1 1 "skimmed_combined_vtx_"$det"_B_"$BF"T_vtx_0.05XX0_10um_FastSimEval.root"
			./copy_vertex 1 1 1 "skimmed_combined_vtx_"$det"_B_"$BF"T_vtx_0.05XX0_10um_FastSimEval.root"
		done

		# Do one last time without updating table
		./copy_momentum 1 2 1 "skimmed_combined_vtx_"$det"_B_"$BF"T_vtx_0.05XX0_10um_FastSimEval.root"
		./copy_vertex 1 2 1 "skimmed_combined_vtx_"$det"_B_"$BF"T_vtx_0.05XX0_10um_FastSimEval.root"
	done
done

say "All done!"
