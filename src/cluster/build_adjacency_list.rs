use crate::config::Config;
//use crate::config::{Config,NeighborCutoff};
use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::physical_constants::PI;
use crate::quantum::tensors::HamiltonianTensors;
use crate::signal::calculate_analytic_restricted_2cluster_signals::{
  hahn_three_spin_modulation_depth,
  hahn_three_spin_modulation_frequency,
};

pub fn build_adjacency_list(tensors: &HamiltonianTensors, config: &Config) 
  -> Result<AdjacencyList, CluEError>
{

  let n_spins = tensors.len();
  let mut adjacency_list = AdjacencyList::with_capacity(n_spins);

  for idx0 in 1..n_spins{
    adjacency_list.connect(idx0,idx0);


    for idx1 in idx0+1..n_spins{

      if are_spins_neighbors(idx0,idx1, tensors,config){
        adjacency_list.connect(idx0,idx1);
      }

    }
  }
  Ok(adjacency_list)
}
//------------------------------------------------------------------------------
fn are_spins_neighbors(idx0: usize,idx1: usize,
    tensors: &HamiltonianTensors, config: &Config) -> bool
{

  if tensors.spin_multiplicities[idx0] != tensors.spin_multiplicities[idx1]{
    return false;
  }

  // TODO: add toggle in config {
  let Some(zeeman0) = tensors.spin1_tensors.get(idx0) else {
    return false;
  };

  let Some(zeeman1) = tensors.spin1_tensors.get(idx1) else {
    return false;
  };

  let delta_zeeman = zeeman1 - zeeman0;
  let sum_zeeman = zeeman1 + zeeman0;

  if 2.0*delta_zeeman.norm()/sum_zeeman.norm() > 1e-12 {
    return false;
  }
  // TODO: }

  let Some(dipdip) = tensors.spin2_tensors.get(idx0,idx1) else {
    return false;
  };

  let Some(hf0) = tensors.spin2_tensors.get(0,idx0) else {
    return false;
  };
  
  let Some(hf1) = tensors.spin2_tensors.get(0,idx1) else {
    return false;
  };

  let delta_hf = (hf0.zz() -hf1.zz()).abs();
  let b = 0.5*( dipdip.xx() + dipdip.yy() ).abs();

  if let Some(cutoff) = &config.neighbor_cutoff_delta_hyperfine{
    if delta_hf < *cutoff {return false;}
  }
  if let Some(cutoff) = &config.neighbor_cutoff_dipole_dipole{
    if b < *cutoff {return false;}
  }
  if let Some(cutoff) = &config.neighbor_cutoff_3_spin_hahn_mod_depth{
    let k = hahn_three_spin_modulation_depth(delta_hf,b);
    if k < *cutoff {return false;}
  }
  if let Some(cutoff) = &config.neighbor_cutoff_3_spin_hahn_taylor_4{
    let k = hahn_three_spin_modulation_depth(delta_hf,b);
    let omega = 2.0*PI*hahn_three_spin_modulation_frequency(delta_hf,b);
    let kom4 = k*omega.powi(4);
    if kom4 < *cutoff {return false;}
  }
  true
}
