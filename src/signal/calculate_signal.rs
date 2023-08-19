
use crate::build_adjacency_list;
use crate::find_clusters;
use crate::config::{Config,ClusterMethod,OrientationAveraging};
use crate::CluEError;
use crate::cluster::methyl_clusters::partition_cluster_set_by_exchange_groups;
use crate::cluster::find_clusters::ClusterSet;
use crate::HamiltonianTensors;
use crate::integration_grid::IntegrationGrid;
use crate::physical_constants::{ONE,PI};
use crate::signal::Signal;
use crate::signal::write_vec_signals;
use crate::signal::cluster_correlation_expansion::*;
use crate::signal::calculate_analytic_restricted_2cluster_signals::{
  calculate_analytic_restricted_2cluster_signals};
use crate::Structure;
use crate::quantum::spin_hamiltonian::*;
use crate::math;
use crate::space_3d::UnitSpherePoint;

use num_complex::Complex;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
use rand_distr::{Distribution,Uniform};

//------------------------------------------------------------------------------
/// This function averages of `number_system_instances` instances of building
/// a system and calculating the signal.
pub fn calculate_signals(rng: &mut ChaCha20Rng, config: &Config,
    path_opt: &Option<String>) -> Result<Vec::<Signal>,CluEError>
{
  let Some(max_cluster_size) = config.max_cluster_size else {
    return Err(CluEError::NoMaxClusterSize);
  };
  let Some(number_system_instances) = config.number_system_instances else {
    return Err(CluEError::NoNumberSystemInstances);
  }; 

  let range = Uniform::from(0..std::u64::MAX);
  let new_seeds = (0..number_system_instances).map(|_ii|
      range.sample(rng) )
    .collect::<Vec::<u64>>();

  let n_tot = config.number_timepoints.iter().sum::<usize>();

  let mut signals = (0..max_cluster_size).map(|_ii| Signal::zeros(n_tot))
    .collect::<Vec::<Signal>>();

  for ii in 0..number_system_instances{

    println!("\nSystem instance: {}/{}.", ii+1, number_system_instances);

    let mut inst_rng = ChaCha20Rng::seed_from_u64(new_seeds[ii]);
    let inst_signals = calculate_structure_signal(&mut inst_rng, 
        config, path_opt)?;

    for (order_idx, sig) in inst_signals.iter().enumerate(){
      signals[order_idx] = &signals[order_idx] + sig;
    }
  }

  for sig in signals.iter_mut(){
    sig.scale(ONE/(number_system_instances as f64));
  }
  Ok(signals)
}
//------------------------------------------------------------------------------
// This function builds the spin structure and calculates the signal over
// one or multiple orientations.
fn calculate_structure_signal(rng: &mut ChaCha20Rng, config: &Config,
    path_opt: &Option<String>) -> Result<Vec::<Signal>,CluEError>
{

  // Build input structure.
  let structure = Structure::build_structure(rng,config)?;
  

  // Determine if/where to save results.
  let save_dir_opt = match path_opt{
    Some(path) => {
      
      let structure_hash = math::str_hash(&structure);
      
      let save_dir = format!("{}/system-{}",path,structure_hash);
      match std::fs::create_dir_all(save_dir.clone()){
        Ok(_) => (),
        Err(_) => return Err(CluEError::CannotCreateDir(save_dir)),
      }
  
      if let Some(info_dir) = &config.write_info{

        let info_path = format!("{}/{}",save_dir,info_dir);

        if std::fs::create_dir_all(info_path.clone()).is_err(){
          return Err(CluEError::CannotCreateDir(info_path));
        }

        if let Some(filename) = &config.write_bath{
          structure.bath_to_csv(&format!("{}/{}.csv", info_path, filename))?;
        }

        if let Some(filename) = &config.write_structure_pdb{
          structure.write_pdb(&format!("{}/{}.pdb", info_path, filename))?;
        }

        if let Some(exchange_group_manager) = &structure.exchange_groups{
          if let Some(filename) = &config.write_exchange_groups{
            let csv_file = format!("{}/{}.csv",info_path,filename);
            exchange_group_manager.to_csv(&csv_file)?;
          }
        }
      }
      Some(save_dir)
    },
    None => None,
  };



  // Determine maximum cluster size.
  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };

  // Get number of data points per trace.
  let n_tot = config.number_timepoints.iter().sum::<usize>();
  
  // Initialize output.
  let mut order_n_signals = (0..max_cluster_size).map(|_| Signal::zeros(n_tot))
    .collect::<Vec::<Signal>>();

  
  // Generate coupling tensors.
  let tensors = HamiltonianTensors::generate(&structure, config)?;
  if let Some(save_dir) = &save_dir_opt{
    if let Some(info_dir) = &config.write_info{
      let info_path = format!("{}/{}",save_dir,info_dir);
      
      if let Some(tensor_save_name) = &config.write_tensors{
        if std::fs::create_dir_all(info_path.clone()).is_err(){
          return Err(CluEError::CannotCreateDir(info_path));
        }
        let tensor_path = format!("{}/{}.txt",info_path,tensor_save_name);

        tensors.save(&tensor_path,&structure)?;
      }
    }
  }


  // Determin orientation averaging method.
  let integration_grid = match &config.orientation_grid{
    Some(OrientationAveraging::Grid(grid)) => grid.clone(),
    Some(OrientationAveraging::Lebedev(n_ori)) 
      => IntegrationGrid::lebedev(*n_ori)?.remove_3d_hemisphere(),
    Some(OrientationAveraging::Random(n_ori)) 
      => IntegrationGrid::random_unit_sphere(*n_ori,rng),
    None => IntegrationGrid::z_3d(),
  };

  // Loop over orientation.
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

  if let Some(save_dir) =  &save_dir_opt{
    let save_path = format!("{}/signal.csv", save_dir);
    let headers = (1..=max_cluster_size)
      .map(|ii| format!("signal_{}",ii))
      .collect::<Vec::<String>>();
    write_vec_signals(&order_n_signals, headers, &save_path)?;
  }
  Ok(order_n_signals)
}
//------------------------------------------------------------------------------
// This function calculates the signal for the input structure at the input
// orientation.
fn calculate_signal_at_orientation(rot_dir: UnitSpherePoint,
    mut tensors: HamiltonianTensors, structure: &Structure,
    config: &Config, path_opt: &Option<String>)
  -> Result<Vec::<Signal>,CluEError>
{
  
  let theta_degrees = rot_dir.theta()*180.0/PI;
  let phi_degrees = rot_dir.phi()*180.0/PI;
  println!(
      "\nMagnetic field orientation: θ = {} degrees; φ = {} degrees.", 
      theta_degrees, phi_degrees);
  // Determine if/where to save results.
  let save_dir_opt = match path_opt{
    Some(path) => {

      if let Some(ori_path) = &config.write_orientation_signals{ 
      
        let save_dir = format!("{}/{}/theta_{}deg_phi_{}deg",path, ori_path,
            theta_degrees, phi_degrees);
      
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
  

  // Rotate the coupling tensor to the specified orientation.
  tensors.rotate_pasive(&rot_dir);

  // Determine adjacencies.
  let adjacency_list = build_adjacency_list(&tensors, structure, config)?;

  // Determine maximum cluster size.
  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };

  // Find clusters.
  let mut cluster_set = find_clusters(&adjacency_list, max_cluster_size)?;

  // Remove partial methyls.
  let Some(do_remove_partial_methyls) = &config.remove_partial_methyls else{
    return Err(CluEError::NoRemovePartialMethyls);
  };
  if *do_remove_partial_methyls{
    cluster_set.remove_partial_methyls(structure)?;
  }

  for (size_idx,clusters_of_size) in cluster_set.clusters.iter().enumerate(){
    println!("Found {} clusters of size {}.",clusters_of_size.len(),size_idx+1);
  }

  if let Some(path) = &save_dir_opt{
    if let Some(cluster_file) = &config.write_clusters{
      let cluster_save_path = format!("{}/{}.txt",path,cluster_file);
      cluster_set.save(&cluster_save_path,structure)?;
    }
  }


  // Calculate cluster signals.
  let order_n_signals = match &config.cluster_method{
    Some(ClusterMethod::AnalyticRestricted2CCE) => 
      calculate_analytic_restricted_2cluster_signals(&mut cluster_set, &tensors,
          config,&save_dir_opt)?,
    Some(ClusterMethod::CCE) => {
      let spin_multiplicity_set =
        math::unique(tensors.spin_multiplicities.clone());
      let spin_ops = ClusterSpinOperators::new(&spin_multiplicity_set,
          max_cluster_size)?;
      do_cluster_correlation_expansion(&mut cluster_set, &spin_ops, &tensors, 
          config,&save_dir_opt)?
    },
    None => return Err(CluEError::NoClusterMethod)
  };


  if let Some(save_dir) =  &save_dir_opt{
    let save_path = format!("{}/signal.csv", save_dir);
    let headers = (1..=max_cluster_size)
      .map(|ii| format!("signal_{}",ii))
      .collect::<Vec::<String>>();
    write_vec_signals(&order_n_signals, headers, &save_path)?;
  }

  // TODO: add toggle in config
  caculate_sans_spin_signals(&order_n_signals[max_cluster_size - 1],
    &cluster_set, structure, config, &save_dir_opt)?;


  // TODO: add toggle in config
  calculate_methyl_partition_cce(cluster_set, structure, config, 
      &save_dir_opt)?;

  Ok(order_n_signals)
} 
//------------------------------------------------------------------------------
// For each active spin in structure, this function calculates the signal,
// where the spin is dropped. 
fn caculate_sans_spin_signals(ref_signal: &Signal,
    cluster_set: &ClusterSet, structure: &Structure, config: &Config, 
    save_dir_opt: &Option<String>)
  -> Result<(),CluEError>
{
  
  let save_path = match save_dir_opt{
    None => return Ok(()),
    Some(path) 
      => format!("{}/sans_spin_signals.csv", path),
  };

  let mut spin_signals 
    = caculate_bath_spin_contributions(cluster_set, structure, config)?;

  for sig in spin_signals.iter_mut(){
    *sig = ref_signal / sig;
  }

  let headers = structure.bath_particles.iter().enumerate()
    .filter_map(|(idx,particle)| 
        if particle.active{
          Some( format!("sans_{}",idx) )
        }else{
          None
        }
    ).collect::< Vec::<String> >();

  write_vec_signals(&spin_signals,headers,&save_path)
}
//------------------------------------------------------------------------------
// For each active spin in structure, this function calculates the product of
// auxiliary signals for clusters that contain the spin.
// Note that since clusters contain multiple spins, the contributions from
// multiple spins are not independent of each other, 
// meaning that the product of contributions will not in general equal the 
// calculated signal.
fn caculate_bath_spin_contributions(
    cluster_set: &ClusterSet, structure: &Structure, config: &Config)
  -> Result<Vec::<Signal>,CluEError>
{

  // Get number of data points per trace.
  let n_tot = config.number_timepoints.iter().sum::<usize>();

  let mut spin_contributions = (0..structure.number_active())
    .map(|_| Signal::ones(n_tot)).collect::<Vec::<Signal>>();

  let mut index = 0;
  let vertex_to_index = structure.bath_particles.iter()
    .filter_map(|particle| 
        if particle.active{
          let idx = index;
          index += 1;
          Some(idx)
        }else{
          None
        }  
    ).collect::< Vec::<usize> >();

  for clusters_of_size in cluster_set.clusters.iter(){

    for cluster in clusters_of_size{

      for &vertex in cluster.vertices(){

        let idx = vertex_to_index[vertex-1];

        match &cluster.signal{
          Ok(Some(signal)) => 
            spin_contributions[idx] = &spin_contributions[idx] * signal,
          Ok(None) => (),
          Err(err) => return Err(err.clone()),
        }
      } 
    }

  }

  Ok(spin_contributions)
}
//------------------------------------------------------------------------------
fn calculate_methyl_partition_cce(
    cluster_set: ClusterSet, structure: &Structure, config: &Config, 
    save_dir_opt: &Option<String>)
  -> Result<(),CluEError>
{
  if let Some(exchange_group_manager) = &structure.exchange_groups{
    let cluster_partitions = partition_cluster_set_by_exchange_groups(
        cluster_set, exchange_group_manager, structure)?;

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
            => return Err(CluEError::CannotCreateDir(save_dir)),
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
