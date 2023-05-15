use crate::clue_errors::CluEError;
use crate::cluster::Cluster;
use crate::cluster::find_clusters::ClusterSet;
use crate::math;
use crate::structure::exchange_groups::ExchangeGroupManager;

use std::collections::HashMap;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn partition_cluster_set_by_exchange_groups(cluster_set: ClusterSet,
    exchange_group_manager: &ExchangeGroupManager)
  -> Result<HashMap::<String,ClusterSet>,CluEError>
{

  let max_size = cluster_set.clusters.len();


  // Count the number of partitions. 
  let n_groups = exchange_group_manager.exchange_groups.len() + 1;

  let mut n_cluster_partitions 
    = HashMap::<String,Vec::<usize>>::with_capacity(n_groups);
  for clu_size in 0..max_size{
    for (idx,cluster) in cluster_set.clusters[clu_size].iter().enumerate(){
      let key = get_cluster_partition_key(cluster,exchange_group_manager);

      let mut counts: Vec::<usize> = vec![0;max_size];
      if let Some(n) = n_cluster_partitions.get(&key){
        counts = n.clone();
      }
      counts[clu_size] += 1;

      n_cluster_partitions.insert(key,counts);
    }
  }


  // Collect the partitions.
  let n_groups = n_cluster_partitions.len();
  let mut cluster_partitions 
    = HashMap::<String,Vec::<Vec<Cluster>>>::with_capacity(n_groups);

    for clu_size in 0..cluster_set.clusters.len(){
      for (idx,cluster) in cluster_set.clusters[clu_size].iter().enumerate(){
        let key = get_cluster_partition_key(cluster,exchange_group_manager);

        let Some(counts) = n_cluster_partitions.get(&key) else
        {
          return Err(CluEError::InvalidClusterPartitionKey);
        };

        if let Some(clusters) = cluster_partitions.get_mut(&key){
          clusters[clu_size].push(cluster.clone());
        }else{
          let mut clusters = (0..max_size)
            .map(|size| Vec::<Cluster>::with_capacity(counts[size]))
            .collect::<Vec::<Vec::<Cluster>>>();

          clusters[clu_size].push(cluster.clone());
          cluster_partitions.insert(key.clone(), clusters);
        }
      }
    }


  // Transformation: Cluster -> ClusterSet.
  let mut cluster_set_partitions 
    = HashMap::<String,ClusterSet>::with_capacity(n_groups);
  
  for (key,clusters) in cluster_partitions{
    let cluster_set = ClusterSet::from(clusters);
    cluster_set_partitions.insert(key,cluster_set);
  }


  Ok(cluster_set_partitions)
}
//------------------------------------------------------------------------------
fn get_cluster_partition_key(cluster: &Cluster,
    exchange_group_manager: &ExchangeGroupManager) -> String
{
  let mut key_ids = Vec::<usize>::with_capacity(cluster.len());

  for &vertex in cluster.vertices().iter(){
    
    if let Some(id) = exchange_group_manager.exchange_group_ids[vertex]{
      key_ids.push(id);
    } 
  }

  key_ids = math::unique(key_ids);

  if key_ids.is_empty(){
    return String::from("bath");
  }

  let mut key = exchange_group_manager.exchange_groups[key_ids[0]].to_string();
  for (ii,&id) in key_ids.iter().skip(1).enumerate(){
    key = format!("{}_{}",key,
        exchange_group_manager.exchange_groups[id].to_string());
  }

  key

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// TODO: generalize to other exchange groups
/// This function removes all clusters from a cluster set that contain partial
/// methyl group, leaving only clusters that contain either no methyl hydrogens 
/// or all three hydrogens from any methyl.
impl ClusterSet{
  pub fn remove_partial_methyls(&mut self,
      exchange_group_manager: &ExchangeGroupManager){

    for clu_size in 0..self.clusters.len(){

      let n_clusters = self.clusters[clu_size].len();
      let mut to_keep = Vec::<usize>::with_capacity(n_clusters);
      
      for (idx,cluster) in self.clusters[clu_size].iter().enumerate(){
        if cluster.contains_partial_methyl(&exchange_group_manager){ continue; }
        to_keep.push(idx);
      }

      let mut kept_clusters = Vec::<Cluster>::with_capacity(to_keep.len());
      let mut kept_cluster_indices = HashMap::with_capacity(to_keep.len());
      for (new_idx, old_idx) in to_keep.iter().enumerate(){
        let cluster = self.clusters[clu_size][*old_idx].clone();
        kept_cluster_indices.insert(cluster.vertices.clone(),new_idx);
        kept_clusters.push(cluster);
      }

      self.clusters[clu_size] = kept_clusters;
      self.cluster_indices[clu_size] = kept_cluster_indices;
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// TODO: generalize to other exchange groups
impl Cluster{
pub fn contains_partial_methyl(&self,exchange_group_manager: 
    &ExchangeGroupManager) -> bool
{
  let mut count_map = HashMap::<usize,usize>::with_capacity(self.len());
  for &idx in self.vertices().iter(){
    match exchange_group_manager.exchange_group_ids[idx]{
      Some(key) => {
        let mut count: usize = 0;
        if let Some(n) = count_map.get(&key){
          count = *n;
        }
        count += 1;

        count_map.insert(key,count);
      },
      None => (),
    }
  }

  for (_key, count) in count_map.iter(){
    if count % 3 != 0 {
      return true;
    }
  }
  false
}
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::space_3d::Vector3D;
  use crate::structure::exchange_groups::{C3Rotor,ExchangeGroup};
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_partition_cluster_set_by_exchange_groups(){
    let exchange_group_manager =  generate_test_exchange_group_manager();

    let clusters = vec![
      vec![Cluster::from(vec![0]),Cluster::from(vec![4])],
      vec![Cluster::from(vec![0,1]),Cluster::from(vec![0,4]),
        Cluster::from(vec![4,8])],
      vec![Cluster::from(vec![0,1,2]),Cluster::from(vec![0,3,4]),
        Cluster::from(vec![3,6,7])],
      vec![Cluster::from(vec![0,1,2,3])],
      vec![Cluster::from(vec![0,1,2,3,4])],
      vec![Cluster::from(vec![0,1,2,3,4,5])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7,8])],
    ];

    let mut ref_partitions 
      = HashMap::<String,ClusterSet>::with_capacity(4);


    ref_partitions.insert("bath".to_string(), ClusterSet::from(vec![
      vec![Cluster::from(vec![4])],
      vec![Cluster::from(vec![4,8])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_0_1_2".to_string(), ClusterSet::from(vec![
      vec![Cluster::from(vec![0])],
      vec![Cluster::from(vec![0,1]),Cluster::from(vec![0,4])],
      vec![Cluster::from(vec![0,1,2])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_6_3_7".to_string(), ClusterSet::from(vec![
      vec![],
      vec![],
      vec![Cluster::from(vec![3,6,7])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_0_1_2_methyl_6_3_7".to_string(), 
     ClusterSet::from(vec![
      vec![],
      vec![],
      vec![Cluster::from(vec![0,3,4])],
      vec![Cluster::from(vec![0,1,2,3])],
      vec![Cluster::from(vec![0,1,2,3,4])],
      vec![Cluster::from(vec![0,1,2,3,4,5])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7,8])],
    ]));

    let ref_counts = clusters.iter().map(|v| v.len())
      .collect::<Vec::<usize>>();

    let cluster_set = ClusterSet::from(clusters);

    let cluster_partitions: HashMap::<String,ClusterSet> = 
      partition_cluster_set_by_exchange_groups(cluster_set,
          &exchange_group_manager).unwrap();


    let mut counts = ref_counts.iter().map(|x| 0)
      .collect::<Vec::<usize>>();

    for (_key, value_cluster_set) in cluster_partitions.iter(){
      for (clu_size,clusters) in value_cluster_set.clusters.iter().enumerate(){
        counts[clu_size] += clusters.len();
      }
    }
    assert_eq!(counts,ref_counts);


    for (ref_key,ref_cluster_set) in ref_partitions.iter(){
      let cluster_set = &cluster_partitions[ref_key];
      assert_eq!(*cluster_set, *ref_cluster_set);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_cluster_partition_key(){
    let exchange_group_manager =  generate_test_exchange_group_manager();

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![0,1,2]),
        &exchange_group_manager), "methyl_0_1_2".to_string());

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![0,1,2,4]),
        &exchange_group_manager), "methyl_0_1_2".to_string());

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![4]),
        &exchange_group_manager), "bath".to_string());

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![1]),
        &exchange_group_manager), "methyl_0_1_2".to_string());

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![0,1,2,3]),
        &exchange_group_manager), "methyl_0_1_2_methyl_6_3_7".to_string());

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![3,2]),
        &exchange_group_manager), "methyl_0_1_2_methyl_6_3_7".to_string());

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_remove_partial_methyls(){
    let exchange_group_manager =  generate_test_exchange_group_manager();

    let clusters = vec![
      vec![Cluster::from(vec![0]),Cluster::from(vec![4])],
      vec![Cluster::from(vec![0,1]),Cluster::from(vec![0,4]),
        Cluster::from(vec![4,8])],
      vec![Cluster::from(vec![0,1,2]),Cluster::from(vec![0,3,4]),
        Cluster::from(vec![3,6,7])],
      vec![Cluster::from(vec![0,1,2,3])],
      vec![Cluster::from(vec![0,1,2,3,4])],
      vec![Cluster::from(vec![0,1,2,3,4,5])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7,8])],
    ];

    let mut cluster_set = ClusterSet::from(clusters);

    let clusters = vec![
      vec![Cluster::from(vec![4])],
      vec![Cluster::from(vec![4,8])],
      vec![Cluster::from(vec![0,1,2]),Cluster::from(vec![3,6,7])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7])],
      vec![Cluster::from(vec![0,1,2,3,4,5,6,7,8])],
    ];


    cluster_set.remove_partial_methyls(&exchange_group_manager);
    let ref_cluster_set = ClusterSet::from(clusters);
    assert_eq!(cluster_set,ref_cluster_set);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_contains_partial_methyl(){
    let exchange_group_manager =  generate_test_exchange_group_manager();

    let cluster = Cluster::from(vec![4]);
    assert!(!cluster.contains_partial_methyl(&exchange_group_manager));

    let cluster = Cluster::from(vec![0,1,2]);
    assert!(!cluster.contains_partial_methyl(&exchange_group_manager));

    let cluster = Cluster::from(vec![0]);
    assert!(cluster.contains_partial_methyl(&exchange_group_manager));

    let cluster = Cluster::from(vec![1,2]);
    assert!(cluster.contains_partial_methyl(&exchange_group_manager));

    let cluster = Cluster::from(vec![0,1,3]);
    assert!(cluster.contains_partial_methyl(&exchange_group_manager));
   
    let cluster = Cluster::from(vec![0,1,2,3,4,5,6,7,8]);
    assert!(!cluster.contains_partial_methyl(&exchange_group_manager));
  }
  //----------------------------------------------------------------------------
  fn generate_test_exchange_group_manager() -> ExchangeGroupManager{
    
    let exchange_groups = vec![
      ExchangeGroup::Methyl(C3Rotor{
        center: Vector3D::from([0.0,0.0,0.0]),
        normal: Vector3D::from([0.0,0.0,1.0]),
        indices: [0,1,2]
      }),
      ExchangeGroup::Methyl(C3Rotor{
        center: Vector3D::from([1.0,1.0,1.0]),
        normal: Vector3D::from([1.0,0.0,0.0]),
        indices: [6,3,7]
      })
    ];

    let exchange_group_ids = vec![
      Some(0),Some(0),Some(0),
      Some(1),None, None, 
      Some(1),Some(1),None];

    let exchange_couplings = exchange_group_ids.iter()
      .map(|_id| 0.0).collect::<Vec::<f64>>();

    ExchangeGroupManager{
      exchange_groups,
      exchange_group_ids,
      exchange_couplings
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
