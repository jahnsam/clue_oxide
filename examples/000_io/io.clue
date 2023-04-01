
input_structure_file = "../../assets/TEMPO_wat_gly_70A.pdb";
radius = 12e-10; // m.
detected_spin_position = centroid_over_serials([28,29]);
//write_structure_pdb = out.pdb;

pulse_sequence = Hahn;
number_timepoints = [101];
time_increments = [1e-7];
cluster_method = r2CCE;
max_cluster_size = 2;
magnetic_field = 1.2;
