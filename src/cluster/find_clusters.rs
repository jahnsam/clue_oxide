use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::math;
use crate::structure::Structure;

use crate::cluster::Cluster;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// The `ClusterSet` is a data structure that holds clusters of various sizes
/// as well as indices to find the cluster.
/// The first feild is `clusters: Vec::<Vec::<Cluster>>`, where
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

    for cluster_of_size in self.clusters.iter(){
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

}
//------------------------------------------------------------------------------
/// This function finds all connected graphs (clusters) with up to `max_size`
/// vertices.
pub fn find_clusters( adjacency_list: &AdjacencyList, max_size: usize) 
  -> Result< ClusterSet, CluEError>
{

  let mut clusters = Vec::<Vec::<Cluster>>::with_capacity(max_size);
  let mut cluster_indices 
    = Vec::<HashMap<Vec::<usize>,usize>>::with_capacity(max_size);

  let vertices = adjacency_list.get_active_vertices();
  let n_one_clusters = vertices.len();
  let mut one_clusters = Vec::<Cluster>::with_capacity(n_one_clusters); 
  let mut one_cluster_indices = HashMap::new(); 

  // Identify all 1-clusters, as those with adjacency_matrix[[ii,ii]] == true.
  for vertex in vertices.into_iter(){

    let key = vec![vertex];
    let cluster = Cluster::from(key.clone());
    one_cluster_indices.insert(key, one_clusters.len());
    one_clusters.push(cluster);
  }


  clusters.push(one_clusters);
  cluster_indices.push(one_cluster_indices);

  for clu_size in 1..max_size{

    // Build n-clusters from (n-1)-clusters.
    if let Ok(n_cluster_set)
     = build_n_clusters(&clusters[clu_size -1],adjacency_list){
      let ClusterSet{clusters: nn_clusters,cluster_indices: nn_cluster_indices} 
       = n_cluster_set; 
      
      if nn_clusters.len() != 1 {
        return Err(CluEError::ExpectedClusterSetWithNSizes(1,
              nn_clusters.len()) );
      }
      if nn_cluster_indices.len() != 1{
        return Err(CluEError::ExpectedClusterSetWithNSizes(1,
              nn_cluster_indices.len()) );
      }
      for n_clusters in nn_clusters{ clusters.push(n_clusters);}
      for n_cluster_indices in nn_cluster_indices{
        cluster_indices.push(n_cluster_indices);
      }

    }else{
        return Err(CluEError::NoClustersOfSize(clu_size+1));
    }


  }

  Ok(ClusterSet{
    clusters,
    cluster_indices,
    }) 
}

//------------------------------------------------------------------------------
// This function takes a set of valid (n-1)-clusters and an `&AdjacencyList`,
// and finds all the n-cluster that contain at least one of the (n-1)-clusters.
fn build_n_clusters(
    n_minus_1_clusters: &Vec::<Cluster>, 
    adjacency_list: &AdjacencyList)
  -> Result<ClusterSet,CluEError>
{

  let mut new_clusters = Vec::<Cluster>::new();
  let mut new_cluster_indices = HashMap::new();
  if n_minus_1_clusters.is_empty(){ 
    return Ok( ClusterSet{ 
      clusters: vec![new_clusters],
      cluster_indices: vec![new_cluster_indices],
    } );
  }

  for cluster in n_minus_1_clusters.iter() {

    for idx in cluster.vertices.iter(){
      if let Some(neighbors) = adjacency_list.get_neighbors(*idx){

        for vertex in neighbors.iter(){
      
          let mut new_indices: Vec::<usize> = cluster.vertices.clone();
          new_indices.push(*vertex);
          new_indices = math::unique(new_indices);

          if new_indices.len() == cluster.vertices.len(){ continue; }
        
          new_indices.sort();

          let  new_cluster = Cluster::from(new_indices.clone());
          if new_cluster_indices.get(&new_indices).is_some(){
            continue;
          }
          new_cluster_indices.insert(new_indices,new_clusters.len());
          new_clusters.push(new_cluster);
        } 
      }
    }
  }

  Ok( ClusterSet{ 
    clusters: vec![new_clusters],
    cluster_indices: vec![new_cluster_indices],
    } )
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  //----------------------------------------------------------------------------
  #[test]
  fn test_find_clusters(){
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

    let mut cluster_set = find_clusters(&cube,9).unwrap();
    let clusters = &mut cluster_set.clusters;
    let cluster_indices = &cluster_set.cluster_indices;

    assert!(clusters[8].is_empty());
    clusters.pop();
    assert_eq!(clusters.len(),8);

    assert_eq!(clusters[0].len(), 8);
    for ii in 0..8{
      let v = Vec::<usize>::from([ii]);
      let idx = cluster_indices[0][&v];
      assert_eq!(clusters[0][idx].vertices,*v);
    }
    
    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    let two_clusters = vec![
      vec![0,1],/*vec![0,2],*/vec![0,3],vec![0,4],//vec![0,5],vec![0,6],vec![0,7],
      vec![1,2],/*vec![1,3],vec![1,4],*/vec![1,5],//vec![1,6],vec![1,7],
      vec![2,3],/*vec![2,4],vec![2,5],*/vec![2,6],//vec![2,7],
      /*vec![3,4],vec![3,5],vec![3,6],*/vec![3,7],
      vec![4,5],/*vec![4,6],*/vec![4,7],
      vec![5,6],//vec![5,7],
      vec![6,7],
    ];

    assert_eq!(clusters[1].len(), two_clusters.len());
    for v in two_clusters.iter(){
      let idx = cluster_indices[1][v];
      assert_eq!(clusters[1][idx].vertices,*v);
    } 

    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    let three_clusters = vec![
      vec![0,1,2],vec![0,1,3],vec![0,1,4],vec![0,1,5],//vec![0,1,6],vec![0,1,7],
      vec![0,2,3],//vec![0,2,4],vec![0,2,5],vec![0,2,6],vec![0,2,7],
      vec![0,3,4],/*vec![0,3,5],vec![0,3,6],*/vec![0,3,7],
      vec![0,4,5],/*vec![0,4,6],*/vec![0,4,7],
      //vec![0,5,6],vec![0,5,7],
      //vec![0,6,7],
      vec![1,2,3],/*vec![1,2,4],*/vec![1,2,5],vec![1,2,6],//vec![1,2,7],
      //vec![1,3,4],vec![1,3,5],vec![1,3,6],vec![1,3,7],
      vec![1,4,5],//vec![1,4,6],vec![1,4,7],
      vec![1,5,6],//vec![1,5,7],
      //vec![1,6,7],
      /*vec![2,3,4],vec![2,3,5],*/vec![2,3,6],vec![2,3,7],
      //vec![2,4,5],vec![2,4,6],vec![2,4,7],
      vec![2,5,6],//vec![2,5,7],
      vec![2,6,7],
      /*vec![3,4,5],vec![3,4,6],*/vec![3,4,7],
      //vec![3,5,6],vec![3,5,7],
      vec![3,6,7],
      vec![4,5,6],vec![4,5,7],
      vec![4,6,7],
      vec![5,6,7],
    ];

    let size_idx = 2;
    assert_eq!(clusters[size_idx].len(), three_clusters.len());
    for v in three_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 

    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    let four_clusters = vec![
      vec![0,1,2,3],vec![0,1,2,4],vec![0,1,2,5],vec![0,1,2,6],//vec![0,1,2,7],
      vec![0,1,3,4],vec![0,1,3,5],/*vec![0,1,3,6],*/vec![0,1,3,7],
      vec![0,1,4,5],/*vec![0,1,4,6],*/vec![0,1,4,7],
      vec![0,1,5,6],//vec![0,1,5,7],
      //vec![0,1,6,7],
      vec![0,2,3,4],/*vec![0,2,3,5],*/vec![0,2,3,6],vec![0,2,3,7],
      //vec![0,2,4,5],vec![0,2,4,6],vec![0,2,4,7],
      //vec![0,2,5,6],vec![0,2,5,7],
      //vec![0,2,6,7],
      vec![0,3,4,5],/*vec![0,3,4,6],*/vec![0,3,4,7],
      //vec![0,3,5,6],vec![0,3,5,7],
      vec![0,3,6,7],
      vec![0,4,5,6],vec![0,4,5,7],
      vec![0,4,6,7],
      //vec![0,5,6,7],
      /*vec![1,2,3,4],*/vec![1,2,3,5],vec![1,2,3,6],vec![1,2,3,7],
      vec![1,2,4,5],//vec![1,2,4,6],vec![1,2,4,7],
      vec![1,2,5,6],//vec![1,2,5,7],
      vec![1,2,6,7],
      //vec![1,3,4,5],vec![1,3,4,6],vec![1,3,4,7],
      //vec![1,3,5,6],vec![1,3,5,7],
      //vec![1,3,6,7],
      vec![1,4,5,6],vec![1,4,5,7],
      //vec![1,4,6,7],
      vec![1,5,6,7],
      /*vec![2,3,4,5],vec![2,3,4,6],*/vec![2,3,4,7],
      vec![2,3,5,6],//vec![2,3,5,7],
      vec![2,3,6,7],
      vec![2,4,5,6],//vec![2,4,5,7],
      vec![2,4,6,7],
      vec![2,5,6,7],
      /*vec![3,4,5,6],*/vec![3,4,5,7],
      vec![3,4,6,7],
      vec![3,5,6,7],
      vec![4,5,6,7],
    ];

    let size_idx = 3;
    assert_eq!(clusters[size_idx].len(), four_clusters.len());
    for v in four_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 

    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    // 0,-,2,-,-,5,6,7
    // -,1,-,3,4,-,6,7
    // 0,-,2,-,4,5,-,7
    // 0,1,-,3,-,5,6,-
    // -,1,2,3,4,-,6,-
    // 0,-,2,3,-,5,-,7
    // 0,1,-,3,4,-,6,-
    // 0,1,2,-,-,5,-,7
    let five_clusters = vec![
      vec![0,1,2,3,4],vec![0,1,2,3,5],vec![0,1,2,3,6],vec![0,1,2,3,7],
      vec![0,1,2,4,5],vec![0,1,2,4,6],vec![0,1,2,4,7],
      vec![0,1,2,5,6],//vec![0,1,2,5,7],
      vec![0,1,2,6,7],
      vec![0,1,3,4,5],/*vec![0,1,3,4,6],*/vec![0,1,3,4,7],
      vec![0,1,3,5,6],vec![0,1,3,5,7],
      vec![0,1,3,6,7],
      vec![0,1,4,5,6],vec![0,1,4,5,7],
      vec![0,1,4,6,7],
      vec![0,1,5,6,7],
      vec![0,2,3,4,5],vec![0,2,3,4,6],vec![0,2,3,4,7],
      vec![0,2,3,5,6],//vec![0,2,3,5,7],
      vec![0,2,3,6,7],
      vec![0,2,4,5,6],//vec![0,2,4,5,7],
      vec![0,2,4,6,7],
      //vec![0,2,5,6,7],
      vec![0,3,4,5,6],vec![0,3,4,5,7],
      vec![0,3,4,6,7],
      vec![0,3,5,6,7],
      vec![0,4,5,6,7],
      vec![1,2,3,4,5],/*vec![1,2,3,4,6],*/vec![1,2,3,4,7],
      vec![1,2,3,5,6],vec![1,2,3,5,7],
      vec![1,2,3,6,7],
      vec![1,2,4,5,6],vec![1,2,4,5,7],
      vec![1,2,4,6,7],
      vec![1,2,5,6,7],
      /*vec![1,3,4,5,6],*/vec![1,3,4,5,7],
      //vec![1,3,4,6,7],
      vec![1,3,5,6,7],
      vec![1,4,5,6,7],
      vec![2,3,4,5,6],vec![2,3,4,5,7],
      vec![2,3,4,6,7],
      vec![2,3,5,6,7],
      vec![2,4,5,6,7],
      vec![3,4,5,6,7],
    ];

    let size_idx = 4;
    assert_eq!(clusters[size_idx].len(), five_clusters.len());
    for v in five_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 


    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    let six_clusters = vec![
      vec![0,1,2,3,4,5],vec![0,1,2,3,4,6],vec![0,1,2,3,4,7],
      vec![0,1,2,3,5,6],vec![0,1,2,3,5,7],
      vec![0,1,2,3,6,7],
      vec![0,1,2,4,5,6],vec![0,1,2,4,5,7],
      vec![0,1,2,4,6,7],
      vec![0,1,2,5,6,7],
      vec![0,1,3,4,5,6],
      vec![0,1,3,4,6,7],
      vec![0,1,3,4,6,7],
      vec![0,1,3,5,6,7],
      vec![0,1,4,5,6,7],
      vec![0,2,3,4,5,6], vec![0,2,3,4,5,7],
      vec![0,2,3,4,6,7],
      vec![0,2,3,5,6,7],
      vec![0,2,4,5,6,7],
      vec![0,3,4,5,6,7],
      vec![1,2,3,4,5,6],vec![1,2,3,4,5,7],
      vec![1,2,3,4,6,7],
      vec![1,2,3,5,6,7],
      vec![1,2,4,5,6,7],
      vec![1,3,4,5,6,7],
      vec![2,3,4,5,6,7],
    ];

    let size_idx = 5;
    assert_eq!(clusters[size_idx].len(), six_clusters.len());
    for v in six_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 
    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    let seven_clusters = vec![
      vec![0,1,2,3,4,5,6],
      vec![0,1,2,3,4,5,7],
      vec![0,1,2,3,4,6,7],
      vec![0,1,2,3,5,6,7],
      vec![0,1,2,4,5,6,7],
      vec![0,1,3,4,5,6,7],
      vec![0,2,3,4,5,6,7],
      vec![1,2,3,4,5,6,7],
    ];

    let size_idx = 6;
    assert_eq!(clusters[size_idx].len(), seven_clusters.len());
    for v in seven_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 

    // 0-------1
    // | \   / |
    // |  4-5  |
    // |  | |  |
    // |  7-6  |
    // | /   \ |
    // 3-------2
    let eight_clusters = vec![
      vec![0,1,2,3,4,5,6,7],
    ];

    let size_idx = 7;
    assert_eq!(clusters[size_idx].len(), eight_clusters.len());
    for v in eight_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 
  }
  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_remove_subclusters_of(){

    let mut square = AdjacencyList::with_capacity(4);
    square.connect(0,1);
    square.connect(1,2);
    square.connect(2,3);
    square.connect(3,0);
    
    let mut cluster_set = find_clusters(&square,4).unwrap();
    remove_subclusters_of(&mut cluster_set, &Cluster::from(vec![0,1]));
    let clusters = cluster_set.clusters;
    let cluster_indices = cluster_set.cluster_indices;

    let one_clusters = vec![
      //vec![0], vec![1],
      vec![2], vec![3],
    ];

    let size_idx = 0;
    assert_eq!(clusters[size_idx].len(), one_clusters.len());
    for v in one_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 

    let two_clusters = vec![
      vec![0,1], //vec![0,2], vec![0,3]
      //vec![1,2], vec![1,3],
      vec![2,3],
    ];

    let size_idx = 1;
    assert_eq!(clusters[size_idx].len(), two_clusters.len());
    for v in two_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 

    let three_clusters = vec![
      vec![0,1,2], vec![0,1,3], 
      //vec![0,2,3], 
      // vec![1,2,3],
    ];

    let size_idx = 2;
    assert_eq!(clusters[size_idx].len(), three_clusters.len());
    for v in three_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 

    let four_clusters = vec![
      vec![0,1,2,3], 
    ];

    let size_idx = 3;
    assert_eq!(clusters[size_idx].len(), four_clusters.len());
    for v in four_clusters.iter(){
      let idx = cluster_indices[size_idx][v];
      assert_eq!(clusters[size_idx][idx].vertices,*v);
    } 
  }
  */
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
