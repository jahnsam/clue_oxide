rng_seed = 42;

input_structure_file = "VO(TCPP-Zn2-bpy).pdb";
radius = 40; // angstroms.

detected_spin_position = [8.322, 8.322, 4.407];
detected_spin_g_matrix = [1.9865, 1.9865, 1.963];
detected_spin_g_x = [1,0,0];
detected_spin_g_y = [0,1,0];

number_timepoints = [101];

time_increments = [5e-9]; // seconds.
cluster_method = cce;
max_cluster_size = 2;
magnetic_field = 0.35; // tesla.
pulse_sequence = hahn;

orientation_grid = lebedev(302);

neighbor_cutoff_3_spin_hahn_mod_depth = 3.23e-8; 
neighbor_cutoff_3_spin_hahn_taylor_4 = 1e14; // (rad/4)^4.

clash_distance_pbc = 0.1; // angstroms.

#[filter(label = mof_v)]
  elements in [V];
  serials in [38];
  cell_ids in [0];

#[spin_properties(label = mof_v, isotope = 51V)]
  hyperfine_coupling = [168,168,477]*1e6; // HZ
  hyperfine_x = [1,0,0];
  hyperfine_y = [0,1,0];


#[filter(label = mof_v2)]
  elements in [V];
  serials in [39];

// Remove the duplicate V.  
#[spin_properties(label = mof_v2, isotope = 51V)]
  active = false;
