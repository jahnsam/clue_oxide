use crate::Config;
use crate::cluster::find_clusters::ClusterSet;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::CluEError;
use crate::physical_constants::ONE;
use crate::math;

use rayon::prelude::*;
use std::collections::HashMap;
use num_complex::Complex;

pub fn calculate_analytic_restricted_2cluster_signals(
    cluster_set: &mut ClusterSet, 
    tensors: &HamiltonianTensors, config: &Config,
    ) ->Result<(), CluEError>
{

  
  let clusters = &mut cluster_set.clusters;

  if clusters.len() != 2 {
    return Err(CluEError::WrongClusterSizeForAnalyticCCE(clusters.len()));
  }

  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };

  let n_clusters = clusters[1].len();
  let n_batches = math::ceil( (n_clusters as f64)/(batch_size as f64)) as usize;

  for ibatch in 0..n_batches{
    let idx = ibatch*batch_size;
    clusters[1].par_iter_mut().skip(idx).take(batch_size).for_each(
        |cluster| {
          let aux_signal = analytic_restricted_2cluster_signal(
            &cluster.vertices(), tensors, config).unwrap(); // TODO: fix unwrap
          (*cluster).signal = Some(aux_signal);
        });
  }

  Ok(())
}

fn analytic_restricted_2cluster_signal(vertices: &Vec::<usize>,
    tensors: &HamiltonianTensors, config: &Config) -> Result<Signal,CluEError>
{
  if vertices.len() != 2 {
    return Err(CluEError::WrongClusterSizeForAnalyticCCE(vertices.len()));
  }
  let idx0 = vertices[0];
  let idx1 = vertices[1];

  let Some(dipdip) = &tensors.spin2_tensors.get(idx0,idx1) else {
    return Ok(Signal::new());
  };
  let Some(hf0) =  &tensors.spin2_tensors.get(0,idx0) else {
    return Ok(Signal::new());
  }; 
  let Some(hf1) =  &tensors.spin2_tensors.get(0,idx1) else {
    return Ok(Signal::new());
  }; 

  let b = dipdip.zz();
  let delta_hf = (hf0.zz() - hf1.zz()).abs();
  
  let omega = (delta_hf*delta_hf + b*b).sqrt();
  let k = (2.0*b*delta_hf/( omega*omega )).powi(2);

  let time_axis = config.get_time_axis()?;
  let nt = time_axis.len();
  let mut data = Vec::<Complex<f64>>::with_capacity(nt);

  for t in time_axis.iter(){
    let s4 = (omega*(*t)).sin().powi(4);
    data.push(ONE - k*s4 );
  }

  Ok(Signal{data})


}

