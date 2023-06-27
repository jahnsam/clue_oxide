input_structure_file = "MD_300K_25ns-30ns_center_0001.pdb";
radius = 25;

detected_spin_position = centroid_over_serials([28,29]);
detected_spin_g_matrix = [2.0097, 2.0064, 2.0025];
detected_spin_g_y = diff(filter(tempo_n) , filter(tempo_o) );
detected_spin_g_x = diff(filter(tempo_c1) , filter(tempo_c19) );
number_timepoints = [101];

time_increments = [5e-8];
cluster_method = cce;
max_cluster_size = 2;
magnetic_field = 1.2;
pulse_sequence = hahn;

orientation_grid = lebedev(170);

neighbor_cutoff_3_spin_hahn_mod_depth = 3.23e-5; 
neighbor_cutoff_3_spin_hahn_taylor_4 = 1e17;


#[filter(label = tempo_n)]
  residues in [TEM];
  elements in [N];

#[filter(label = tempo_o)]
  residues in [TEM];
  elements in [O];

#[filter(label = tempo_c1)]
  serials in [1];

#[filter(label = tempo_c19)]
  serials in [19];
  
#[spin_properties(label = tempo_n, isotope = 14N)]
  hyperfine_coupling = [20,20,100]*1e6;
  hyperfine_x = diff(bonded(tempo_c1),bonded(tempo_c19));
  hyperfine_y = diff(particle,bonded(tempo_o));

  electric_quadrupole_coupling = [-0.28, -1.47, 1.75]*1e6;
  electric_quadrupole_x = diff(bonded(tempo_c1),bonded(tempo_c19));
  electric_quadrupole_y = diff(particle,bonded(tempo_o));



