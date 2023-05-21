//! CluE (Cluster Evolution) is a spin dynamics simulation program for
//! electron spin decoherence.
//!
//! CluE implements Yang an Liu's cluster correlation expansion for an input
//! structure, usually a PDB file.  Simulation details are specified in a config
//! file.
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
pub mod io;


use crate::config::Config;
use crate::clue_errors::CluEError;
use crate::structure::Structure;
use crate::quantum::tensors::HamiltonianTensors;
use crate::cluster::build_adjacency_list::build_adjacency_list;
use crate::cluster::find_clusters::find_clusters;
use crate::signal::calculate_signal;

use num_complex::Complex;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

/// This function follows the config instruction to central simulate spin 
/// decoherence.
pub fn run(config: Config) 
  -> Result<( Vec::<f64>, Vec::<Complex::<f64>> ),CluEError>
{


  let mut rng: ChaCha20Rng;
  match config.rng_seed{
    Some(seed) => rng = ChaCha20Rng::seed_from_u64(seed),
    None => rng = ChaCha20Rng::from_entropy(),
  }

  let config_hash = math::str_hash(&config);
  let mut save_path = String::new();
  if let Some(root_dir) = &config.root_dir{
    save_path = root_dir.to_string();
  }

  if let Some(save_name) = &config.save_name{
    if save_name.is_empty(){
      return Err(CluEError::SaveNameEmpty);
    }
    save_path = format!("{}/{}{}",save_path,save_name,config_hash);
  }else{
    return Err(CluEError::SaveNameNotSet);
  };

  match std::fs::create_dir_all(save_path.clone()){
    Ok(_) => (),
    Err(_) => return Err(CluEError::CannotCreateDir(save_path)),
  }

  config.write_time_axis(save_path.clone())?;

  let order_n_signals 
    = calculate_signal::calculate_structure_signal(&mut rng, &config,
      &Some(save_path) )?;


  let time_axis = config.get_time_axis()?;

  let max_size = order_n_signals.len();

  Ok((time_axis.clone(), order_n_signals[max_size-1].data.clone()))
}

