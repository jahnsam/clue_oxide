// TODO: 
//  Each convergence function should have the option to skip convergence
//  and return the next level's result.
//
// TODO:
//  All funcs/converge_signal_for_clusters_cutoffs() should have the ability to 
//  return the cluster cutoffs so that the orientation average can pick one set
//  of cutoffs to apply to every orientation.
//
// TODO:
//  Allow saving of all, some, or only the final signal.
//  BASENAME_sructure_config_grid_ori_sstate_nCCE_cutoffs_cluster.csv
//  format!("BASENAME_str{}_cfg{}_grd{}_ori{}_ss{}_CCE{}_r{}_hf{}_dip{}_clu{}.csv",
//                "TEMPO",cfg_seed,lebedev28,13,bath_seed, 4.  12,  3e4,  1e3)
//  if name.len() > max_len{
//  name = hash_func(name);
//  save_key_value_to_file(); }
// 
// TODO: 
//  Implement a checkpoint system to restart an interupted simulation.
//  Generate systematic name, before calculations and if it exists,
//  load the data rather than recalculating.
//
//
//------------------------------------------------------------------------------
// This averages the spin decoherence signal over multiple input structures,
// multiple input PDBs for example.
fn average_structure_signal(config: &Config) -> Result<Signal,CluEError>{
} 
//------------------------------------------------------------------------------
// This function converges the spin decoherence signal over different 
// configurations of the same base structure, a Monte Carlo search over
// different isotopologues for example.
fn converge_structure_configuration_signal(config: &Config) 
  -> Result<Signal,CluEError>
{

  let mut signal = Signal::ones(config.n_timepoints);
  let mut signal0 = Signal::zeros(config.n_timepoints);
  let mut cutoffs = config.initial_cutoff.clone();

  let mut counter = 0;
  loop{
  
    let structure = build_structure(config);

    signal0 = signal;
    (signal,cutoffs) = converge_orientation_averaged_signal(structure,config);

    if !do_converge{ break; }

    let is_converged = signal.approx_eq(signal0,config.method_of_comparison);

    counter += 1;
    if is_converged || counter >= config.converge_structure.max_loops{
      break;
    }
  }

  Ok(signal);
}
//------------------------------------------------------------------------------
// This function converges the orientation grid for a particular structure
// configuration.
fn converge_orientation_averaged_signal(config: &Config) 
  -> Result<Signal,CluEError>
{
  loop{
    (Signal,cutoffs) = calculate_orientation_averaged_signal();
  }
  Ok(Signal,cutoffs)
} 
//------------------------------------------------------------------------------
// This function averages the spin decoherence signal over all orientaions
// of an input quadrature grid, such as a Lebedev grid.
fn calculate_orientation_averaged_signal(config: &Config) 
  -> Result<Signal,CluEError>
{
  (Signal,cutoffs) = converge_state_signal(config); 
  Ok(Signal,cutoffs)
} 
//------------------------------------------------------------------------------
// This function converges the spin decoherence signal of multiple initial
// states.
fn converge_state_signal(config: &Config) 
  -> Result<Signal,CluEError>
{
  mean_field_tensors = get_mean_feilds();
  (Signal,cutoffs) = converge_signal_for_clusters_cutoffs(tensors,config);
  Ok(Signal,cutoffs)
} 
//------------------------------------------------------------------------------
// This function converges the cluster cutoffs: 
// cluster size, neighbor cutoffs, and vertex cutoffs.
fn converge_signal_for_clusters_cutoffs(tensors, config: &Config) 
  -> Result<Signal,CluEError>
{

  loop{
    cutoffs.next()
    let clusters = find_clusters(tensors,cutoffs);
    signal = calculate_cluster_signal();
  }

  Ok(Signal,cutoffs)
} 
//------------------------------------------------------------------------------