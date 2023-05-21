use crate::config::{Config,ClusterMethod};
use crate::CluEError;
use crate::cluster::methyl_clusters::partition_cluster_set_by_exchange_groups;
use crate::cluster::find_clusters::ClusterSet;
use crate::Structure;
use crate::HamiltonianTensors;
use crate::build_adjacency_list;
use crate::find_clusters;
use crate::signal::Signal;
use crate::signal::cluster_correlation_expansion::*;
use crate::signal::calculate_analytic_restricted_2cluster_signals::{
  calculate_analytic_restricted_2cluster_signals};
use rand_chacha::ChaCha20Rng;
use crate::quantum::spin_hamiltonian::*;
use crate::math;
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
pub fn average_structure_signal(rng: &mut ChaCha20Rng, config: &Config,
    path: &str) -> Result<(),CluEError>
{
  let structure = Structure::build_structure(rng,config)?;
  
  let structure_hash = math::str_hash(&structure);

  let save_dir = format!("{}/system-{}",path,structure_hash);
  match std::fs::create_dir_all(save_dir.clone()){
    Ok(_) => (),
    Err(_) => return Err(CluEError::CannotCreateDir(save_dir)),
  }

  if let Some(filename) = &config.write_structure_pdb{
    structure.write_pdb(&format!("{}/{}.pdb",save_dir, filename))?;
  }

  let tensors = HamiltonianTensors::generate(&structure, config)?;

  let adjacency_list = build_adjacency_list(&tensors, config)?;

  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };

  let mut cluster_set = find_clusters(&adjacency_list, max_cluster_size)?;

  // TODO: add toggle in config for remove_partial_methyls
  if let Some(exchange_group_manager) = &structure.exchange_groups{
    cluster_set.remove_partial_methyls(exchange_group_manager);
  }

  for (size_idx,clusters_of_size) in cluster_set.clusters.iter().enumerate(){
    println!("Found {} clusters of size {}.",clusters_of_size.len(),size_idx+1);
  }

  let Some(cluster_method) = &config.cluster_method else {
    return Err(CluEError::NoClusterMethod);
  };

  match cluster_method{
    ClusterMethod::AnalyticRestricted2CCE => 
      calculate_analytic_restricted_2cluster_signals(&mut cluster_set, &tensors,
          config,&None)?,
    ClusterMethod::CCE => {
      let spin_multiplicity_set =
        math::unique(tensors.spin_multiplicities.clone());
      let spin_ops = ClusterSpinOperators::new(&spin_multiplicity_set,
          max_cluster_size)?;
      do_cluster_correlation_expansion(&mut cluster_set, &spin_ops, &tensors, 
          config,&None)?;
    },
  }

  // TODO: add toggle in config
  calculate_methyl_partition_cce(cluster_set, &structure, config, &save_dir)?;

  Ok(())
} 
//------------------------------------------------------------------------------
fn calculate_methyl_partition_cce(
    cluster_set: ClusterSet, structure: &Structure, config: &Config, 
    save_dir: &String)
  -> Result<(),CluEError>
{
  if let Some(exchange_group_manager) = &structure.exchange_groups{
    let cluster_partitions = partition_cluster_set_by_exchange_groups(
        cluster_set, exchange_group_manager)?;

    let part_dir = "methyl_partitions".to_string();

    let n_tot = config.number_timepoints.iter().sum::<usize>();

    for (key_name, value_cluster_set) in cluster_partitions.iter(){

      let mut part_signal = Signal::ones(n_tot);

      for size_idx in 0..(*value_cluster_set).len() {
        for cluster in value_cluster_set.clusters[size_idx].iter(){
          if let Ok(Some(aux_signal)) = &cluster.signal{
            part_signal = &part_signal * aux_signal;
          }
        }
      }

      let part_save_path = format!("{}/{}/{}.csv",
          save_dir,part_dir,key_name);
      part_signal.write_to_csv(&part_save_path)?;
    }
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
