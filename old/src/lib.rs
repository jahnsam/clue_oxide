//! CluE Oxide (Cluster Evolution Oxide) is a spin dynamics simulation program 
//! for electron spin decoherence.
//!
//! CluE Oxide implements Yang an Liu's cluster correlation expansion for an
//! input structure, usually a PDB file.  
//! 
//! Documentation on CluE Oxide the program is provided in the 
//![manual](https://github.com/jahnsam/clue_oxide/blob/main/manual/CluE_Oxide.pdf).
pub mod config;
pub mod clue_errors;
pub mod cluster;
pub mod info;
pub mod integration_grid;
pub mod physical_constants;
pub mod quantum;
pub mod structure;
pub mod signal;
pub mod space_3d;
pub mod symmetric_list_2d;
pub mod math;
pub mod kmeans;
pub mod io;



use crate::config::Config;
use crate::clue_errors::CluEError;
use crate::structure::Structure;
use crate::quantum::tensors::HamiltonianTensors;
use crate::cluster::build_adjacency_list::build_adjacency_list;
use crate::cluster::find_clusters::find_clusters;
use crate::signal::calculate_signal;
use crate::signal::write_vec_signals;

use num_complex::Complex;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

/// This function follows the config instruction to central simulate spin 
/// decoherence.
pub fn run(config: Config) 
  -> Result<( Vec::<f64>, Vec::<Complex::<f64>> ),CluEError>
{

  let mut rng = match config.rng_seed{
    Some(seed) => ChaCha20Rng::seed_from_u64(seed),
    None => ChaCha20Rng::from_entropy(),
  };

  let root_path = if let Some(root_dir) = &config.root_dir{
    root_dir.to_string()
  }else{
    String::from("./")
  };

  let save_path_opt = match &config.save_name{
    Some(save_name) => {
      if save_name.is_empty(){
        None
      }else{
        Some(format!("{}/{}", root_path, save_name))
      }
    },
    None => {
     let config_hash = math::str_hash(&config);
    Some(format!("{}/CluE-{}", root_path, config_hash))
    },
  };

  if let Some(save_path) = &save_path_opt{
    match std::fs::create_dir_all(save_path.clone()){
      Ok(_) => (),
      Err(_) => return Err(CluEError::CannotCreateDir(save_path.to_string())),
    }

    config.write_time_axis(save_path.clone())?;
  }

  let order_n_signals 
    = calculate_signal::calculate_signals(&mut rng, &config,
      &save_path_opt )?;


  let time_axis = config.get_time_axis()?;

  let max_size = order_n_signals.len();

  if let Some(save_dir) =  &save_path_opt{
    let save_path = format!("{}/signal.csv", save_dir);
    let headers = (1..=max_size)
      .map(|ii| format!("signal_{}",ii))
      .collect::<Vec::<String>>();
    write_vec_signals(&order_n_signals, headers, &save_path)?;
  }

  Ok((time_axis.clone(), order_n_signals[max_size-1].data.clone()))
}

