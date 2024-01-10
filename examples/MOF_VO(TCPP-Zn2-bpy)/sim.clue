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

cluster_density_matrix = thermal(5);

orientation_grid = lebedev(302);

neighbor_cutoff_3_spin_hahn_mod_depth = 3.23e-8; 
neighbor_cutoff_3_spin_hahn_taylor_4 = 1e14; // (rad/4)^4.

clash_distance_pbc = 0.1; // angstroms.

#[group(mof_v)]
  primary_cell;
  elements in [V];
  serials in [38];

#[spin_properties(mof_v, 51V)]
  hyperfine_coupling = [168,168,477]*1e6; // HZ
  hyperfine_x = [1,0,0];
  hyperfine_y = [0,1,0];


#[group(mof_v2)]
  elements in [V];
  serials in [39];

// Remove the duplicate V.  
#[spin_properties(mof_v2, 51V)]
  active = false;
  
