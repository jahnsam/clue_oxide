use crate::config::{Config,ClusterMethod};
use crate::CluEError;
use crate::Structure;
use crate::HamiltonianTensors;
use crate::build_adjacency_list;
use crate::find_clusters;
use crate::signal::calculate_cluster_signal::calculate_cluster_signals;
use crate::signal::calculate_analytic_restricted_2cluster_signals::{
  calculate_analytic_restricted_2cluster_signals};
use crate::signal::cluster_correlation_expansion::*;
use rand_chacha::ChaCha20Rng;
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
pub fn average_structure_signal(rng: &mut ChaCha20Rng, config: &Config) 
  -> Result<(),CluEError>
{
  let structure = Structure::build_structure(rng,config)?;

  if let Some(filename) = &config.write_structure_pdb{
    structure.write_pdb(&format!("{}.pdb",filename))?;
  }

  let tensors = HamiltonianTensors::generate(&structure, config)?;

  let adjacency_list = build_adjacency_list(&tensors, config)?;

  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };

  let mut cluster_set = find_clusters(&adjacency_list, max_cluster_size)?;
  for (size_idx,clusters_of_size) in cluster_set.clusters.iter().enumerate(){
    println!("Found {} clusters of size {}.",clusters_of_size.len(),size_idx+1);
  }

  let Some(cluster_method) = &config.cluster_method else {
    return Err(CluEError::NoClusterMethod);
  };

  match cluster_method{
    ClusterMethod::AnalyticRestricted2CCE => 
      calculate_analytic_restricted_2cluster_signals(cluster_set, &tensors,
          config)?,
    ClusterMethod::CCE => calculate_cluster_signals(cluster_set, &tensors, 
          config)?,
  }

  //let signal = do_cluster_correlation_expansion_product(&cluster_set,config)?;
  Ok(())
} 
/*
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
*/
