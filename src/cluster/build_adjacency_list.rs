use crate::config::Config;
//use crate::config::{Config,NeighborCutoff};
use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::quantum::tensors::HamiltonianTensors;

pub fn build_adjacency_list(tensors: &HamiltonianTensors, config: &Config) 
  -> Result<AdjacencyList, CluEError>
{

  let n_spins = tensors.len();
  let mut adjacency_list = AdjacencyList::with_capacity(n_spins);

  // TODO: required_spins -> config option
  let required_spins = vec![0];
  for idx0 in 1..n_spins{
    if required_spins.contains(&idx0) {continue; }

    for idx1 in idx0+1..n_spins{
      if required_spins.contains(&idx1) {continue;}

      if are_spins_neighbors(idx0,idx1, tensors,config){
        adjacency_list.connect(idx0,idx1);
      }

    }
  }
  Ok(adjacency_list)
}
//------------------------------------------------------------------------------
pub fn are_spins_neighbors(idx0: usize,idx1: usize,
    tensors: &HamiltonianTensors, config: &Config) -> bool
{

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
    let k = (2.0*b*delta_hf/(delta_hf*delta_hf + b*b)).powi(2);
    if k < *cutoff {return false;}
  }
  if let Some(cutoff) = &config.neighbor_cutoff_3_spin_hahn_taylor_4{
    let k = (2.0*b*delta_hf/(delta_hf*delta_hf + b*b)).powi(2);
    let kom4 = k*(delta_hf*delta_hf + b*b).powi(2);
    if k < *cutoff {return false;}
  }
  /*
  for cutoff in &config.neighbor_cutoffs{
    match cutoff{
      NeighborCutoff::DeltaHyperfine(cutoff) => { 
        if delta_hf < *cutoff {return false;}
      },
      NeighborCutoff::DipoleDipole(cutoff) => {
        if b < *cutoff { return false;}
      },
      NeighborCutoff::HahnThreeSpinModulationDepth(cutoff) => {
        let k = (2.0*b*delta_hf/(delta_hf*delta_hf + b*b)).powi(2);
        if k < *cutoff {return false;}
      },
      NeighborCutoff::HahnThreeSpinFourthOrderTaylorCoefficient(cutoff) => {
        let k = (2.0*b*delta_hf/(delta_hf*delta_hf + b*b)).powi(2);
        let kom4 = k*(delta_hf*delta_hf + b*b).powi(2);
        if kom4 < *cutoff { return false; }
      },
    }
  }
  */
  true
}
