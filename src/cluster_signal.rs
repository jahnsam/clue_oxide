use crate::cluster::Cluster;
use crate::signal::Signal;
use rayon::prelude::*;

use std::collections::HashMap;

pub fn calculate_cluster_signals(
    clusters: &mut Vec::<HashMap::<Vec::<usize>, Cluster>>
    ){

  
  let tensors = get_tensors();

  let max_size = clusters.len();

  for cluster_size in 0..max_size{
    clusters[cluster_size].par_iter_mut().for_each(
        |(_idx, cluster)| (*cluster).signal = calculate_full_signal(&tensors)
        );
  }
}

fn calculate_full_signal(tensors: &Vec::<usize>) -> Option<Signal>{
  let nt = tensors.len();
  Some(Signal::ones(nt))
}

fn get_tensors() -> Vec::<usize>{
  let mut tensors = Vec::<usize>::with_capacity(100);
  for ii in 0..100{
    tensors.push(ii);
  }
  tensors
}
