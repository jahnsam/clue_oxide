
//input_structure_file = "../../assets/TEMPO_wat_gly_70A.pdb";
input_structure_file = "../../assets/TEMPO.pdb";
radius = 18; // m.
detected_spin_position = centroid_over_serials([28,29]);
//write_structure_pdb = "";

neighbor_cutoff_dipole_dipole = 1; // Hz.
pulse_sequence = hahn;
number_timepoints = [101];
time_increments = [1e-7];
cluster_method = cce;
max_cluster_size = 2;
magnetic_field = 1.2;
remove_partial_methyls = true;

#[filter(label = tempo)]
  residues in [TEM];

#[spin_properties(label = tempo, isotope = 1H)]
  tunnel_splitting = 80e3; // Hz.
