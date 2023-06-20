
//input_structure_file = "../../assets/TEMPO_wat_gly_70A.pdb";
input_structure_file = "../../assets/TEMPO.pdb";
radius = 18; // m.
detected_spin_position = centroid_over_serials([28,29]);
//write_structure_pdb = "";
write_tensors = "out_tensors";

neighbor_cutoff_dipole_dipole = 1; // Hz.
pulse_sequence = hahn;
number_timepoints = [101];
time_increments = [1e-7];
cluster_method = cce;
max_cluster_size = 2;
magnetic_field = 1.2;
remove_partial_methyls = true;

//#[filter(label = tempo)]
//  residues in [TEM];

//#[spin_properties(label = tempo, isotope = 1H)]
//  tunnel_splitting = 80e3; // Hz.
#[filter(label = tempo_n)]
  residues in [TEM];
  elements in [N];

#[filter(label = tempo_o)]
  residues in [TEM];
  elements in [O];

#[filter(label = tempo_c)]
  residues in [TEM];
  elements in [C];

#[spin_properties(label = tempo_n, isotope = 14N)]
  hyperfine_coupling = [20,20,100]*1e6;
  //hyperfine_coupling = [-13.8725, -13.8725, 122.329]*1e6;
  hyperfine_x = diff(bonded(tempo_c),bonded(tempo_c));
  hyperfine_y = diff(particle,bonded(tempo_o));

  electric_quadrupole_coupling = [-0.28, -1.47, 1.75]*1e6;
  electric_quadrupole_x = diff(bonded(tempo_c),bonded(tempo_c));
  electric_quadrupole_y = diff(particle,bonded(tempo_o));

