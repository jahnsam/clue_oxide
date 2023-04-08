use crate::Config;
use crate::clue_errors::CluEError;
use crate::cluster::find_clusters::ClusterSet;
use crate::signal::Signal;

/// This function calculates the product of every `Signal` in `cluster_set`.
pub fn do_cluster_correlation_expansion_product(cluster_set: &ClusterSet,
    config: &Config) -> Result<Signal, CluEError>
{

  let n_tot = config.number_timepoints.iter().sum::<usize>();

  let mut signal = Signal::ones(n_tot);

  for clusters in cluster_set.clusters.iter(){

    for cluster in clusters.iter(){
      if let Ok(Some(v)) =  &cluster.signal{
        signal = &signal * v;
      }
    }

  }

  Ok(signal)
}
