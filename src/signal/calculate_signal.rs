use crate::config::{Config,ClusterMethod,OrientationAveraging};
use crate::CluEError;
use crate::cluster::methyl_clusters::partition_cluster_set_by_exchange_groups;
use crate::cluster::find_clusters::ClusterSet;
use crate::Structure;
use crate::HamiltonianTensors;
use crate::integration_grid::IntegrationGrid;
use crate::build_adjacency_list;
use crate::find_clusters;
use crate::signal::Signal;
use crate::signal::write_vec_signals;
use crate::signal::cluster_correlation_expansion::*;
use crate::signal::calculate_analytic_restricted_2cluster_signals::{
  calculate_analytic_restricted_2cluster_signals};
use crate::quantum::spin_hamiltonian::*;
use crate::math;
use crate::space_3d::UnitSpherePoint;

use num_complex::Complex;
use rand_chacha::ChaCha20Rng;
//------------------------------------------------------------------------------
// This averages the spin decoherence signal over multiple input structures,
// multiple input PDBs for example.
pub fn calculate_structure_signal(rng: &mut ChaCha20Rng, config: &Config,
    path_opt: &Option<String>) -> Result<Vec::<Signal>,CluEError>
{
  let structure = Structure::build_structure(rng,config)?;
  

  let save_dir_opt = match path_opt{
    Some(path) => {
      
      let structure_hash = math::str_hash(&structure);
      
      let save_dir = format!("{}/system-{}",path,structure_hash);
      match std::fs::create_dir_all(save_dir.clone()){
        Ok(_) => (),
        Err(_) => return Err(CluEError::CannotCreateDir(save_dir)),
      }
  
      if let Some(filename) = &config.write_structure_pdb{
        structure.write_pdb(&format!("{}/{}.pdb",save_dir, filename))?;
      }
      Some(save_dir)
    },
    None => None,
  };


  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };
  let mut order_n_signals = Vec::<Signal>::with_capacity(max_cluster_size);

  let n_tot = config.number_timepoints.iter().sum::<usize>();
  
  for _size_idx in 0..max_cluster_size{
    order_n_signals.push(Signal::zeros(n_tot));
  }

  let tensors = HamiltonianTensors::generate(&structure, config)?;

  let integration_grid = match config.orientation_averaging{
    Some(OrientationAveraging::Lebedev(n_ori)) 
      => IntegrationGrid::lebedev(n_ori)?,
    None => IntegrationGrid::z_3d(),
  };

  for iori in 0..integration_grid.len(){
    let rot_dir = UnitSpherePoint::from(integration_grid.xyz(iori)?);
    
    let mut ori_sigs = calculate_signal_at_orientation(
        rot_dir,tensors.clone(), &structure,config, &save_dir_opt )?;

    let weight = Complex::<f64>{ re: integration_grid.weight(iori), im: 0.0};
    for size_idx in 0..max_cluster_size{
      ori_sigs[size_idx].scale(weight);
      order_n_signals[size_idx] 
        = &order_n_signals[size_idx] + &ori_sigs[size_idx];
    }

  }

  Ok(order_n_signals)
}
//------------------------------------------------------------------------------
fn calculate_signal_at_orientation(rot_dir: UnitSpherePoint,
    mut tensors: HamiltonianTensors, structure: &Structure,
    config: &Config, path_opt: &Option<String>)
  -> Result<Vec::<Signal>,CluEError>
{
  
  let save_dir_opt = match path_opt{
    Some(path) => {

      if let Some(ori_path) = &config.write_orientation_signals{ 
      
        let save_dir = format!("{}/orientations/theta_{}_phi_{}",path,
            rot_dir.theta(), rot_dir.phi());
      
        match std::fs::create_dir_all(save_dir.clone()){
          Ok(_) => (),
          Err(_) => return Err(CluEError::CannotCreateDir(save_dir)),
        }

        Some(save_dir)
      }else{
        None
      }
    },
    None => None,
  };
  

  tensors.rotate_active(&rot_dir);

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

  let order_n_signals = match cluster_method{
    ClusterMethod::AnalyticRestricted2CCE => 
      calculate_analytic_restricted_2cluster_signals(&mut cluster_set, &tensors,
          config,&save_dir_opt)?,
    ClusterMethod::CCE => {
      let spin_multiplicity_set =
        math::unique(tensors.spin_multiplicities.clone());
      let spin_ops = ClusterSpinOperators::new(&spin_multiplicity_set,
          max_cluster_size)?;
      do_cluster_correlation_expansion(&mut cluster_set, &spin_ops, &tensors, 
          config,&save_dir_opt)?
    },
  };


  if let Some(save_dir) =  &save_dir_opt{
    let save_path = format!("{}/signal.csv", save_dir);
    write_vec_signals(&order_n_signals, &save_path)?;
  }

  // TODO: calculate nuclear contributions


  // TODO: add toggle in config
  calculate_methyl_partition_cce(cluster_set, &structure, config, 
      &save_dir_opt)?;

  Ok(order_n_signals)
} 
//------------------------------------------------------------------------------
fn calculate_methyl_partition_cce(
    cluster_set: ClusterSet, structure: &Structure, config: &Config, 
    save_dir_opt: &Option<String>)
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

      if let Some(save_dir) = save_dir_opt{
        let save_dir = format!("{}/{}",
            save_dir,part_dir);
        match std::fs::create_dir_all(save_dir.clone()){
          Ok(()) => (),
          Err(_) 
            => return Err(CluEError::CannotCreateDir(save_dir.to_string())),
        }
        let part_save_path = format!("{}/{}.csv",
            save_dir,key_name);
        part_signal.write_to_csv(&part_save_path)?;
      }
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
