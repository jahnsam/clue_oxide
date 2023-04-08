use crate::Config;
use crate::clue_errors::CluEError;
 use crate::cluster::find_clusters::ClusterSet;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::math;
use crate::quantum::spin_hamiltonian::*;

use rayon::prelude::*;

use std::collections::HashMap;

pub fn calculate_cluster_signals(
    mut cluster_set: ClusterSet, spin_ops: &ClusterSpinOperators,
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
          .for_each(|cluster| 
              (*cluster).signal = calculate_full_signal(&cluster.vertices(),
                spin_ops,tensors,config)
        );

      }
  }
  Ok(())
}

fn calculate_full_signal(tensor_indices: &Vec::<usize>, 
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config) 
  -> Result<Option<Signal>,CluEError>
{

  let hamiltonian = build_hamiltonian(tensor_indices,spin_ops, tensors,config)?;
  let density_matrix = get_density_matrix(&hamiltonian, config)?;
  let signal = propagate_pulse_sequence(&density_matrix, &hamiltonian, config)?;
  Ok(Some(signal))
}

