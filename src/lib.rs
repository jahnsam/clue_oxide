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
pub mod elements;
pub mod info;
pub mod isotopes;
pub mod integration_grid;
pub mod physical_constants;
pub mod quantum;
pub mod structure;
pub mod signal;
pub mod space_3d;
pub mod symmetric_list_2d;
pub mod math;
pub mod misc;
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
pub fn run(mut config: Config) 
  -> Result<( Vec::<f64>, Vec::<Complex::<f64>> ),CluEError>
{

  config.set_time_axis()?;

  let mut rng = match config.rng_seed{
    Some(seed) => ChaCha20Rng::seed_from_u64(seed),
    None => ChaCha20Rng::from_entropy(),
  };

  let save_path_opt = match &config.output_directory{
    Some(output_directory) => {
      if output_directory.is_empty(){
        None
      }else{
        Some(format!("{}", output_directory))
      }
    },
    None => {
     let config_hash = math::str_hash(&config);
    Some(format!("CluE-{}", config_hash))
    },
  };


  let mut time_axis = config.get_time_axis()?;

  let seconds_to_unit_of_time = match config.unit_of_time_to_seconds{
    Some(unit_of_time_to_seconds) => 1.0/unit_of_time_to_seconds,
    None => return Err(CluEError::NoUnitOfTime),
  };

  for t in time_axis.iter_mut(){
    *t *= seconds_to_unit_of_time; 
  }

  if let Some(save_path) = &save_path_opt{
    match std::fs::create_dir_all(save_path.clone()){
      Ok(_) => (),
      Err(_) => return Err(CluEError::CannotCreateDir(save_path.to_string())),
    }

    io::write_data(&[time_axis.clone()],
        &format!("{}/time_axis.csv",save_path), vec!["time_axis".to_string()])?;
  }

  let order_n_signals 
    = calculate_signal::calculate_signals(&mut rng, &config,
      &save_path_opt )?;

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

