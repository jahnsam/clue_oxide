use crate::Config;
use crate::clue_errors::CluEError;
use crate::cluster::Cluster;
use crate::signal::Signal;
use crate::HamiltonianTensors;

use rayon::prelude::*;

use std::collections::HashMap;

pub fn calculate_cluster_signals(
    clusters: &mut Vec::<HashMap::<Vec::<usize>, Cluster>>, 
    tensors: &HamiltonianTensors, config: &Config,
    ) -> Result<(),CluEError>
{

  

  let max_size = clusters.len();

  for cluster_size in 0..max_size{
    clusters[cluster_size].par_iter_mut().for_each(
        |(_idx, cluster)| 
        (*cluster).signal = calculate_full_signal(tensors,config)
        );
  }
  Ok(())
}

fn calculate_full_signal(tensors: &HamiltonianTensors, config: &Config) 
  -> Option<Signal>
{
  let nt = tensors.len();
  Some(Signal::ones(nt))
}

