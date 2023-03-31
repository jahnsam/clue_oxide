use crate::Config;
use crate::clue_errors::CluEError;
 use crate::cluster::find_clusters::ClusterSet;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::math;

use rayon::prelude::*;

use std::collections::HashMap;

pub fn calculate_cluster_signals(
    mut cluster_set: ClusterSet, 
    tensors: &HamiltonianTensors, config: &Config,
    ) -> Result<(),CluEError>
{

  
  let clusters = &mut cluster_set.clusters;

  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };

  let max_size = clusters.len();

  for cluster_size in 0..max_size{
    let n_clusters = clusters[cluster_size].len();
    let n_batches 
      = math::ceil( (n_clusters as f64)/(batch_size as f64)) as usize;
      for ibatch in 0..n_batches{
        let idx = ibatch*batch_size;
        clusters[cluster_size].par_iter_mut().skip(idx).take(batch_size)
          .for_each(
            |cluster| 
              (*cluster).signal = calculate_full_signal(tensors,config)
        );
      }
  }
  Ok(())
}

fn calculate_full_signal(tensors: &HamiltonianTensors, config: &Config) 
  -> Option<Signal>
{
  let nt = tensors.len();
  Some(Signal::ones(nt))
}

