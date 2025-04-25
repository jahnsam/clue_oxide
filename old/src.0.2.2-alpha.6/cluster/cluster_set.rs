use crate::clue_errors::CluEError;
use crate::structure::Structure;

use crate::cluster::Cluster;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// The `ClusterSet` is a data structure that holds clusters of various sizes
/// as well as indices to find the cluster.
/// The first field is `clusters: Vec::<Vec::<Cluster>>`, where
/// `clusters[0]` is a  vector of 1-clusters, `clusters[1]` a vector of
/// 2-clusters, and so on. 
/// The other field is `cluster_indices: Vec::<HashMap::<Vec::<usize>,usize>>`,
/// where `cluster_indices[0]` is a `HashMap` that has the vertices of 
/// 1-clusters as keys and the index where the cluster with those vertices can
/// be found as the value.
/// For example if (key,value) = (`vec![a,b]`,idx) in `cluster_indices[1]`,
/// then `clusters[1][idx]` has vertices `vec![a,b]`.
#[derive(Debug,Clone,PartialEq)]
pub struct ClusterSet{
  pub clusters: Vec::<Vec::<Cluster>>,
  pub cluster_indices: Vec::<HashMap::<Vec::<usize>,usize>>, 
}
//------------------------------------------------------------------------------
impl ClusterSet{
  //----------------------------------------------------------------------------
  /// This function returns the number of clusters in the set.
  pub fn len(&self) -> usize{ self.clusters.len() }
  //----------------------------------------------------------------------------
  /// This function returns true iff there are no clusters in the set.
  pub fn is_empty(&self) -> bool{ self.clusters.is_empty() }
  //----------------------------------------------------------------------------
  /// This function converts a `Vec::<Vec::<Cluster>>` to a `ClusterSet`. 
  pub fn from(clusters:Vec::<Vec::<Cluster>>) -> Self
  {
    let mut cluster_indices 
      = Vec::<HashMap::<Vec::<usize>,usize>>::with_capacity(clusters.len());

    for clusters_of_size in clusters.iter(){

      let mut new_cluster_indices 
        = HashMap::<Vec::<usize>,usize>::with_capacity(clusters_of_size.len());

      for (idx,cluster) in clusters_of_size.iter().enumerate(){
        new_cluster_indices.insert(cluster.vertices().clone(),idx);
      }

      cluster_indices.push(new_cluster_indices);
    }

    ClusterSet{clusters,cluster_indices}
  }
  //----------------------------------------------------------------------------
  /// This function writes the `ClusterSet` to the supplied file.
  /// The `Structure` is needed to convert the internal clusters vertices
  /// to match the indices across output files.
  pub fn save(&self, filename: &str,structure: &Structure) 
    -> Result<(),CluEError>
  {
    let Ok(file) = File::create(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };

    let max_size = self.clusters.len();
    let n_clusters: usize = self.clusters.iter().map(|c| c.len()).sum();

    let chars_per_line = 4 + 2*max_size;
    let bytes_per_char = 32;
    
    let n_bytes = n_clusters*chars_per_line*bytes_per_char + 3200;

    let mut stream = BufWriter::with_capacity(n_bytes,file);

    let mut line = "#[clusters, number_clusters = [".to_string();
    for (ii,cluster_of_size) in self.clusters.iter().enumerate(){
      if ii == 0{
        line = format!("{}{}",line,cluster_of_size.len());
      }else{
        line = format!("{},{}",line,cluster_of_size.len());
      }
    }
    line = format!("{}] ]\n\n",line);
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    for (size_idx,cluster_of_size) in self.clusters.iter().enumerate(){
      let line = format!("#[cluster_size = {}]\n",size_idx +1);
      if stream.write(line.as_bytes()).is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }
      for cluster in cluster_of_size.iter(){
        let line = format!("{}\n",
            cluster.to_string_result(structure)? );
          if stream.write(line.as_bytes()).is_err(){
            return Err(CluEError::CannotWriteFile(filename.to_string()) );
          }
      }
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function removes all cluster with more than `max_spins` spins.
  /// Note that prune_large_clusters does not care about the nominal
  /// cluster size: whether a methyl group is considered a 1-cluster
  /// or a 3-cluster, prune_large_clusters recognizes it as having 3 spins.
  pub fn prune_large_clusters(&mut self, max_size: usize)
    -> Result<(),CluEError>
  {
    let mut clusters_to_keep 
        = Vec::<Vec::<usize>>::with_capacity(self.clusters.len());

    for size_idx in 0..self.clusters.len(){
      clusters_to_keep.push(
          Vec::<usize>::with_capacity(self.clusters[size_idx].len()));

      for (idx, cluster) in self.clusters[size_idx].iter().enumerate(){
        if cluster.len() <= max_size{
          clusters_to_keep[size_idx].push(idx);
        }
      }
    }

    self.prune_clusters(clusters_to_keep)
  }
  //----------------------------------------------------------------------------
  /// This function prunes clusters from a `ClusterSet` adjusting the 
  /// `cluster_indices` to keep them correct.  
  /// The `clusters_to_keep` input is a `Vec::<Vec::<usize>>`.  
  /// This corresponds to the `clusters Vec::<Vec::<Cluster>>` in the 
  /// `ClusterSet`.
  /// Specifically, if `clusters_to_keep[n][i] = a` for some `i`, then 
  /// `clusters[n][a]` is kept, and if `a not in clusters_to_keep[n]`,
  /// then `clusters[n][a]` is pruned out.
  pub fn prune_clusters(&mut self, clusters_to_keep: Vec::<Vec::<usize>>)
    -> Result<(),CluEError>
  {
    if clusters_to_keep.len() != self.len(){
      return Err(CluEError::CannotPruneClustersMisMatchedSizes);
    }
    for (size_idx, to_keep) in clusters_to_keep.iter().enumerate(){  

      let mut kept_clusters = Vec::<Cluster>::with_capacity(to_keep.len());
      let mut kept_cluster_indices = HashMap::with_capacity(to_keep.len());

      for (new_idx, old_idx) in to_keep.iter().enumerate(){
        let cluster = self.clusters[size_idx][*old_idx].clone();

        kept_cluster_indices.insert(cluster.vertices.clone(),new_idx);
        kept_clusters.push(cluster);
      }

      self.clusters[size_idx] = kept_clusters;
      self.cluster_indices[size_idx] = kept_cluster_indices;
    }

    Ok(())
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  //use super::*;
  use crate::find_clusters;
  use crate::cluster::adjacency::AdjacencyList;

  //----------------------------------------------------------------------------
  #[test]
  fn test_prune_large_clusters(){
    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    let mut cube = AdjacencyList::with_capacity(8);
    cube.connect(0,1);
    cube.connect(1,2);
    cube.connect(2,3);
    cube.connect(3,0);

    cube.connect(4,5);
    cube.connect(5,6);
    cube.connect(6,7);
    cube.connect(7,4);

    cube.connect(0,4);
    cube.connect(1,5);
    cube.connect(2,6);
    cube.connect(3,7);

    let mut cluster_set = find_clusters(&cube,4).unwrap();

    assert_eq!(cluster_set.clusters[0].len(), 8);
    assert_eq!(cluster_set.clusters[1].len(), 12);
    assert_eq!(cluster_set.clusters[2].len(), 24);
    assert_eq!(cluster_set.clusters[3].len(), 38);

    cluster_set.prune_large_clusters(4).unwrap();

    assert_eq!(cluster_set.clusters[0].len(), 8);
    assert_eq!(cluster_set.clusters[1].len(), 12);
    assert_eq!(cluster_set.clusters[2].len(), 24);
    assert_eq!(cluster_set.clusters[3].len(), 38);

    cluster_set.prune_large_clusters(3).unwrap();

    assert_eq!(cluster_set.clusters[0].len(), 8);
    assert_eq!(cluster_set.clusters[1].len(), 12);
    assert_eq!(cluster_set.clusters[2].len(), 24);
    assert_eq!(cluster_set.clusters[3].len(), 0);

    cluster_set.prune_large_clusters(2).unwrap();

    assert_eq!(cluster_set.clusters[0].len(), 8);
    assert_eq!(cluster_set.clusters[1].len(), 12);
    assert_eq!(cluster_set.clusters[2].len(), 0);
    assert_eq!(cluster_set.clusters[3].len(), 0);

    let clusters_0 = cluster_set.clusters[0].clone();
    cluster_set.clusters[0] = cluster_set.clusters[1].clone();
    cluster_set.clusters[1] = clusters_0;

    assert_eq!(cluster_set.clusters[0].len(), 12);
    assert_eq!(cluster_set.clusters[1].len(), 8);

    cluster_set.prune_large_clusters(1).unwrap();

    assert_eq!(cluster_set.clusters[0].len(), 0);
    assert_eq!(cluster_set.clusters[1].len(), 8);

  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
