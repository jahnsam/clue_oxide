use crate::Config;
use crate::clue_errors::CluEError;
use crate::cluster::{Cluster,
  get_subclusters::build_subclusters,
  cluster_set::ClusterSet,
  unit_of_clustering::UnitOfClustering,
};
use crate::signal::{Signal, load_batch_signals, write_batch_signals};
use crate::structure::Structure;
use crate::HamiltonianTensors;
use crate::math;
use crate::quantum::spin_hamiltonian::*;

use rayon::prelude::*;
use std::path::Path;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function calculates the cluster correlation expansion (CCE)
/// approximation of the echo decay.
/// The CCE method is developed  in
/// W. Yang and R.-B. Liu, “Quantum many-body theory of qubit decoherence
/// in a finite-size spin bath baths,” Phys. Rev. B 78, 085315 (2008).
/// and 
/// W. Yang and R.-B. Liu, “Quantum many-body theory of qubit decoherence
/// in a finite-size spin bath. II. Ensemble dynamics,” Phys. Rev. B 79, 115320
/// (2009).
pub fn do_cluster_correlation_expansion(
    cluster_set: &mut ClusterSet, spin_ops: &ClusterSpinOperators,
    tensors: &HamiltonianTensors, config: &Config, 
    save_path_opt: &Option<String>,structure: &Structure,
    ) -> Result<Vec::<Signal>,CluEError>
{

  
  calculate_auxiliary_signals(cluster_set, spin_ops, tensors, config, 
      save_path_opt,structure)?; 

  let n_tot = config.number_timepoints.iter().sum::<usize>();
  let max_size = cluster_set.clusters.len();


  let mut order_n_signals = Vec::<Signal>::with_capacity(max_size);
  let mut signal = Signal::ones(n_tot);

  for cluster_of_size in cluster_set.clusters.iter(){
    for cluster in cluster_of_size.iter(){
      match &cluster.signal{
        Ok(Some(aux_sig)) => signal = &signal * aux_sig,
        Ok(None) => (),
        Err(err) => return Err(err.clone()),
      }
    }
    order_n_signals.push(signal.clone());
  }

  Ok(order_n_signals)
}
//------------------------------------------------------------------------------
// For each cluster in cluster_set this function replaces cluster.signal
// with that cluster's auxiliary signal.
fn calculate_auxiliary_signals(
    cluster_set: &mut ClusterSet, spin_ops: &ClusterSpinOperators,
    tensors: &HamiltonianTensors, config: &Config, 
    save_path_opt: &Option<String>, structure: &Structure,
    ) 
  -> Result<(),CluEError>
{

  let clusters = &mut cluster_set.clusters;
  let cluster_indices = &cluster_set.cluster_indices;

  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };

  let n_tot = config.number_timepoints.iter().sum::<usize>();
  let max_size = clusters.len();

  // Loop over cluster sizes.
  for cluster_size in 1..=max_size{

    // Split the cluster into batches.
    let n_clusters = clusters[cluster_size-1].len();
    let n_batches 
      = math::ceil( (n_clusters as f64)/(batch_size as f64)) as usize;

    // Loop over batches
    for ibatch in 0..n_batches{
      let idx = ibatch*batch_size;

      // Check for saved data.
      if let Some(path) = &save_path_opt{
        if let Some(aux_dir) = &config.write_auxiliary_signals{
          let load_dir = format!("{}/{}",path, aux_dir);
          let aux_filename = format!("{}/cluster_size_{}_batch_{}.csv",
              load_dir, cluster_size,ibatch);

          if Path::new(&aux_filename).exists(){
             load_batch_signals(&mut clusters[cluster_size-1],
                 idx,batch_size,&aux_filename)?;

             continue;
          }
        }
      }

      // Calculate cluster signals for batch.
      // if !do_gCCE{
      clusters[cluster_size-1].par_iter_mut().skip(idx).take(batch_size)
        .for_each(|cluster| 
            cluster.signal = calculate_cluster_signal(cluster.vertices(),
              spin_ops,tensors,config)
      );
      // } else{
      //clusters[cluster_size-1].par_iter_mut().skip(idx).take(batch_size)
      //  .for_each(|cluster| 
      //      cluster.signal = calculate_general_cluster_signal(cluster.vertices(),
      //        spin_ops,tensors,config)
      //);
      //}
      //


      if cluster_size <= 1 {
        continue;
      }
      // Loop over cluster in this batch and calculate the auxiliary signals
      // from the cluster signals.
      for iclu in (0..n_clusters).skip(idx).take(batch_size){
        let mut cluster = clusters[cluster_size-1][iclu].clone();
        let cluster_num_spins = cluster.vertices().len();

        // Find subclusters.
        let subclusters = build_subclusters(cluster.vertices());

        // Unpack signal.
        let aux_signal: &mut Signal = match &mut cluster.signal{
          Ok(Some(sig)) => sig,
          Ok(None) => return Err(
              CluEError::ClusterHasNoSignal(cluster.to_string())),
          Err(err) => return Err(err.clone()),
        };


        // Loop over subclusters.
        for subcluster_vertices in subclusters.iter(){
          if subcluster_vertices.is_empty(){ continue; }

          let subcluster_num_spins = subcluster_vertices.len();

          // error check
          if subcluster_num_spins >= cluster_num_spins{
            let subcluster = Cluster::from(subcluster_vertices.clone());
            return Err(CluEError::NotAProperSubset(
                  subcluster.to_string(),cluster.to_string()));
          }


          let (start_idx,end_idx) = match config.unit_of_clustering{
            Some(UnitOfClustering::Spin) 
                => (subcluster_num_spins-1,subcluster_num_spins-1),
            Some(UnitOfClustering::Set) 
              => (0,std::cmp::min(subcluster_num_spins-1,cluster_size-2)),
            None => return Err(CluEError::NoUnitOfClustering),
          };
          
          for subcluster_size_idx in start_idx..=end_idx{
            // Devide out subcluster auxiliary signals.
            if let Some(subcluster_idx) 
                = cluster_indices[subcluster_size_idx]
                .get(subcluster_vertices)
            {
              let subcluster = &clusters[subcluster_size_idx][*subcluster_idx];
              match &subcluster.signal{
                Ok(Some(subsignal)) =>  {
                  *aux_signal = &(*aux_signal)/subsignal;
                  break;
                },
                Ok(None) => continue,
                Err(err) => return Err(err.clone()),
              }
            }
          }
        
        }
        clusters[cluster_size-1][iclu] = cluster;
      }


      // Decide if the auxiliary signals should be saved.
      if let Some(path) = &save_path_opt{
        if let Some(aux_dir) = &config.write_auxiliary_signals{
          let save_dir = format!("{}/{}",path, aux_dir);
          match std::fs::create_dir_all(save_dir.clone()){
            Ok(_) => (),
            Err(_) => return Err(CluEError::CannotCreateDir(save_dir)),
          }
          let aux_filename = format!("{}/cluster_size_{}_batch_{}.csv",
              save_dir, cluster_size,ibatch);

          write_batch_signals(&clusters[cluster_size-1], n_tot, idx, batch_size,
              &aux_filename, structure)?; 
        }
      }
    }
  }

  Ok(())
}
//------------------------------------------------------------------------------
// This function calculate the cluster signal for the cluster specified by
// tensor_indices.
fn calculate_cluster_signal(tensor_indices: &[usize], 
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config) 
  -> Result<Option<Signal>,CluEError>
{

  let hamiltonian = build_hamiltonian(tensor_indices,spin_ops, tensors,config)?;
  let density_matrix = get_density_matrix(&hamiltonian, config)?;
  let signal = propagate_pulse_sequence(&density_matrix, &hamiltonian, config)?;
  Ok(Some(signal))
}
//------------------------------------------------------------------------------
// This function calculate the cluster signal for the cluster specified by
// tensor_indices.
/*
fn calculate_general_cluster_signal(tensor_indices: &Vec::<usize>, 
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config) 
  -> Result<Option<Signal>,CluEError>
{

  let hamiltonian = build_general_hamiltonian(tensor_indices,spin_ops, tensors,config)?;
  let density_matrix = get_general_density_matrix(&hamiltonian, config)?;
  let signal = propagate_pulse_sequence_gCCE(&density_matrix, &hamiltonian, config)?;
  Ok(Some(signal))
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  use crate::cluster::adjacency::AdjacencyList;
  use crate::config::PulseSequence;
  use crate::find_clusters;
  use crate::physical_constants::ONE;
  use crate::quantum::tensors::*;
  use crate::signal::calculate_analytic_restricted_2cluster_signals::*;
  use crate::structure::particle::Particle;
  use crate::space_3d::{SymmetricTensor3D,Vector3D};

  //----------------------------------------------------------------------------
  #[test]
  fn test_do_cluster_correlation_expansion(){

    let z0 = 33.0e9;
    let z1 = 80.0e6;
    let a1 = 10.0e6;
    let a2 = -10.0e6;
    let b = 10.0e3;
    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  
    
    let spin_ops = ClusterSpinOperators::new(&vec![2],3).unwrap();

    let mut config = Config::new();
    config.number_timepoints = vec![21];
    let delta_hf = a1 - a2;
    let freq = hahn_three_spin_modulation_frequency(delta_hf,b);
    config.time_increments = vec![0.05/freq];
    config.pulse_sequence = Some(PulseSequence::CarrPurcell(1));
    config.unit_of_clustering = Some(UnitOfClustering::Spin);

    config.set_defaults().unwrap();
    config.construct_time_axis().unwrap();
  

    let mut adjacency_list = AdjacencyList::with_capacity(5);
    adjacency_list.connect(1,2);
    adjacency_list.connect(1,3);
    adjacency_list.connect(1,4);
    adjacency_list.connect(2,3);
    adjacency_list.connect(2,4);
    adjacency_list.connect(3,4);

    let mut cluster_set = find_clusters(&adjacency_list, 3).unwrap();

    let bath_particles = Vec::<Particle>::with_capacity(4);
    let connections = AdjacencyList::with_capacity(4);
    let cell_offsets =  vec![Vector3D::zeros()];
    
    let structure = Structure::new( bath_particles, connections, cell_offsets);

    let order_n_signals = do_cluster_correlation_expansion(&mut cluster_set, 
        &spin_ops, &tensors, &config, &None, &structure).unwrap();

    let spin_indices = vec![1,2];


    let ref_signal_opt = analytic_restricted_2cluster_signal(
        &spin_indices,&tensors,&config).unwrap();
    let Some(ref_signal) = ref_signal_opt else{
      panic!("Could not calculate reference signal.");
    };

    for (ii,v) in order_n_signals[1].data.iter().enumerate(){
      let v0 = ref_signal.data[ii];
      assert!((v-v0*v0).norm() < 1e-12);
      assert!((order_n_signals[0].data[ii]-ONE).norm() < 1e-12);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_calculate_cluster_signal(){


    let z0 = 33.0e9;
    let z1 = 80.0e6;
    let a1 = 10.0e6;
    let a2 = -10.0e6;
    let b = 10.0e3;
    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  
    let spin_indices = vec![1,2];
    
    let spin_ops = ClusterSpinOperators::new(&vec![2],2).unwrap();

    let mut config = Config::new();
    config.number_timepoints = vec![21];
    let delta_hf = a1 - a2;
    let freq = hahn_three_spin_modulation_frequency(delta_hf,b);
    config.time_increments = vec![0.05/freq];
    config.pulse_sequence = Some(PulseSequence::CarrPurcell(1));

    config.set_defaults().unwrap();
    config.construct_time_axis().unwrap();
  
    let signal_opt  = calculate_cluster_signal(&vec![1,2], &spin_ops, &tensors, 
        &config).unwrap();

    let Some(signal) = signal_opt else{
      panic!("Could not calculate signal.");
    }; 
    assert_eq!(signal.data.len(),21);


    let ref_signal_opt = analytic_restricted_2cluster_signal(
        &spin_indices,&tensors,&config).unwrap();
    let Some(ref_signal) = ref_signal_opt else{
      panic!("Could not calculate reference signal.");
    };

    for (ii,v) in signal.data.iter().enumerate(){
      let v0 = ref_signal.data[ii];
      let err: f64;
      if (v+v0).norm() < 1e12{
        err = (v-v0).norm();
      }else{
        err = (2.0+(v-v0)/(v+v0)).norm();
      }
      assert!(err < 1e-9);
    }
  }
  //----------------------------------------------------------------------------
  fn build_restricted_three_spin_tensors(z0: f64, z1: f64, a1: f64, a2: f64, 
      b: f64) -> HamiltonianTensors{
    let spin_multiplicities = vec![2,2,2,2,2];

    let mut spin1_tensors = Spin1Tensors::new(5);
    let zeeman0 = Vector3D::from([0.0, 0.0, z0]);
    let zeeman1 = Vector3D::from([0.0, 0.0, z1]);
    spin1_tensors.set(0,zeeman0);
    spin1_tensors.set(1,zeeman1.clone());
    spin1_tensors.set(2,zeeman1.clone());
    spin1_tensors.set(3,zeeman1.clone());
    spin1_tensors.set(4,zeeman1);

    let mut spin2_tensors = Spin2Tensors::new(5);
    let hf1 = SymmetricTensor3D::from([ 0.0, 0.0, 0.0,
                                             0.0, 0.0,
                                                    a1]);
    let hf2 = SymmetricTensor3D::from([ 0.0, 0.0, 0.0,
                                             0.0, 0.0,
                                                   a2]);

    let dip = SymmetricTensor3D::from([ -b/2.0,    0.0, 0.0,
                                                -b/2.0, 0.0,
                                                         b]);

    spin2_tensors.set(0,1,hf1.clone());
    spin2_tensors.set(0,2,hf2.clone());
    spin2_tensors.set(1,2,dip.clone());

    spin2_tensors.set(0,3,hf1);
    spin2_tensors.set(0,4,hf2);
    spin2_tensors.set(3,4,dip);

    let ge = -1.7609e11;
    HamiltonianTensors{
      spin_multiplicities,
      spin1_tensors,
      spin2_tensors,
      detected_gamma_matrix: SymmetricTensor3D::from([ge, 0.0, 0.0,
                                                           ge, 0.0,
                                                                ge]),
      magnetic_field: Vector3D::from([0.0,0.0,1.2]),
      }
  }  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
