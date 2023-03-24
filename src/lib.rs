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


use crate::config::Config;
use crate::clue_errors::CluEError;
use crate::structure::Structure;

use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

/// This function follows the config instruction to central simulate spin 
/// decoherence.
pub fn run(config: Config) -> Result<(),CluEError>{

  let mut rng: ChaCha20Rng;
  match config.rng_seed{
    Some(seed) => rng = ChaCha20Rng::seed_from_u64(seed),
    None => rng = ChaCha20Rng::from_entropy(),
  }

  let structure = Structure::build_structure(&mut rng,&config)?;
  
  if let Some(filename) = config.write_structure_pdb{
    structure.write_pdb(&format!("{}.pdb",filename))?;
  }

  Ok(())
}

