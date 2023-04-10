use crate::Config;
use crate::clue_errors::CluEError;
 use crate::cluster::find_clusters::ClusterSet;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::math;
use crate::quantum::spin_hamiltonian::*;
use crate::cluster::get_subclusters::build_subclusters;
use crate::cluster::Cluster;

use rayon::prelude::*;

use std::collections::HashMap;

pub fn calculate_cluster_signals(
    mut cluster_set: ClusterSet, spin_ops: &ClusterSpinOperators,
    tensors: &HamiltonianTensors, config: &Config,
    ) -> Result<(),CluEError>
{

  
  let clusters = &mut cluster_set.clusters;
  let cluster_indices = &cluster_set.cluster_indices;

  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };

  let n_tot = config.number_timepoints.iter().sum::<usize>();
  let mut signal = Signal::ones(n_tot);

  let max_size = clusters.len();
  let mut order_n_signals = Vec::<Signal>::with_capacity(max_size);

  for cluster_size in 1..=max_size{
    let n_clusters = clusters[cluster_size-1].len();
    let n_batches 
      = math::ceil( (n_clusters as f64)/(batch_size as f64)) as usize;
      for ibatch in 0..n_batches{
        let idx = ibatch*batch_size;

        clusters[cluster_size-1].par_iter_mut().skip(idx).take(batch_size)
          .for_each(|cluster| 
              (*cluster).signal = calculate_full_signal(&cluster.vertices(),
                spin_ops,tensors,config)
        );


        let idx_end: usize;
        if idx+batch_size < clusters[cluster_size-1].len() {
          idx_end = idx + batch_size;
        }else{
         idx_end = clusters[cluster_size-1].len();
        }
 
        for ii in idx..idx_end{
          let cluster = &clusters[cluster_size-1][ii];
          let mut aux_signal: Signal;
          match &cluster.signal{
            Ok(Some(sig)) => aux_signal = sig.clone(),
            Ok(None) => return Err(
                CluEError::ClusterHasNoSignal(cluster.to_string())),
            Err(err) => return Err(err.clone()),
          }

          let subclusters = build_subclusters(cluster.vertices());

          for subcluster_vertices in subclusters.iter(){
            if subcluster_vertices.is_empty(){ continue; }
            let subcluster_size = subcluster_vertices.len();

            if subcluster_size >= cluster_size{
              let subcluster = Cluster::from(subcluster_vertices.clone());
              return Err(CluEError::NotAProperSubset(
                    subcluster.to_string(),cluster.to_string()));
            }

            if let Some(subcluster_idx) = cluster_indices[subcluster_size-1]
              .get(subcluster_vertices)
            {
              let subcluster = &clusters[subcluster_size-1][*subcluster_idx];
              match &subcluster.signal{
                Ok(Some(subsignal)) =>  aux_signal = &aux_signal/&subsignal,
                Ok(None) => continue,
                Err(err) => return Err(err.clone()),
              }

            }
          
          }

          signal = &signal * &aux_signal;
      }
    }

    order_n_signals.push(signal.clone());
  }
  signal.write_to_csv("out.csv");
  Ok(())
}
//------------------------------------------------------------------------------
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

