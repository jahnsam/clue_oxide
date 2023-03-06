use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::math;

use crate::cluster::Cluster;

use std::collections::HashMap;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn find_clusters(
    adjacency_list: &AdjacencyList, 
    max_size: usize) 
  -> Result< Vec::<HashMap::<Vec::<usize>,Cluster>>, CluEError>
{

  let mut clusters = Vec::<HashMap<Vec::<usize>,Cluster>>::with_capacity(max_size);

  let vertices = adjacency_list.get_active_vertices();
  let mut one_clusters = HashMap::new(); 

  // Identify all 1-clusters, as those with adjacency_matrix[[ii,ii]] == true.
  for vertex in vertices.into_iter(){

    let key = vec![vertex];
    let cluster = Cluster::from(key.clone());
    one_clusters.insert(key, cluster);
  }


  clusters.push(one_clusters);

  for clu_size in 1..max_size{

    // Build n-clusters from (n-1)-clusters.
    if let Ok(n_clusters)
     = build_n_clusters(&mut clusters[clu_size -1],&adjacency_list){
      clusters.push(n_clusters);
    }else{
        return Err(CluEError::NoClustersOfSize(clu_size+1));
    }


  }

  Ok(clusters)
}

//------------------------------------------------------------------------------
fn build_n_clusters(
    n_minus_1_clusters: &HashMap<Vec::<usize>, Cluster>, 
    adjacency_list: &AdjacencyList)
  -> Result<HashMap::<Vec::<usize>,Cluster>,CluEError>
{

  let mut new_clusters = HashMap::new();
  if n_minus_1_clusters.is_empty(){ 
    return Ok(new_clusters);
  }

  for (indices, _cluster) in n_minus_1_clusters {

    for idx in indices.iter(){
      if let Some(neighbors) = adjacency_list.get_neighbors(*idx){

        for vertex in neighbors.iter(){
      
          let mut new_indices: Vec::<usize> = (*indices).clone();
          new_indices.push(*vertex);
          new_indices = math::unique(new_indices);

          if new_indices.len() == indices.len(){ continue; }
        
          new_indices.sort();

          let  new_cluster = Cluster::from(new_indices.clone());
          new_clusters.insert(new_indices,new_cluster);
        } 
      }
    }
  }

  Ok(new_clusters)
}
//------------------------------------------------------------------------------

fn remove_subclusters_of(
    clusters: &mut Vec::<HashMap::<Vec::<usize>,Cluster>>,
    unbreakable_cluster: &Cluster){

  for clu_size in 0..clusters.len(){ 
    let mut to_remove = Vec::<Vec::<usize>>::new();
    
    for (vertices, cluster) in clusters[clu_size].iter() {
      if cluster.overlaps(unbreakable_cluster) 
       && !cluster.contains(unbreakable_cluster) {
        to_remove.push(vertices.clone());
      }
    }

    for vertices in to_remove.into_iter(){
      clusters[clu_size].remove(&vertices);
    }
  }
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

    let mut clusters = find_clusters(&cube,9).unwrap();
    assert!(clusters[8].is_empty());
    clusters.pop();
    assert_eq!(clusters.len(),8);

    assert_eq!(clusters[0].len(), 8);
    for ii in 0..8{
      let v = Vec::<usize>::from([ii]);
      assert_eq!(clusters[0][&v].vertices,v);
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
      assert_eq!(clusters[1][v].vertices,*v);
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

    let idx = 2;
    assert_eq!(clusters[idx].len(), three_clusters.len());
    for v in three_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
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

    let idx = 3;
    assert_eq!(clusters[idx].len(), four_clusters.len());
    for v in four_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
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

    let idx = 4;
    assert_eq!(clusters[idx].len(), five_clusters.len());
    for v in five_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
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

    let idx = 5;
    assert_eq!(clusters[idx].len(), six_clusters.len());
    for v in six_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
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

    let idx = 6;
    assert_eq!(clusters[idx].len(), seven_clusters.len());
    for v in seven_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
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

    let idx = 7;
    assert_eq!(clusters[idx].len(), eight_clusters.len());
    for v in eight_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
    } 
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_remove_subclusters_of(){

    let mut square = AdjacencyList::with_capacity(4);
    square.connect(0,1);
    square.connect(1,2);
    square.connect(2,3);
    square.connect(3,0);
    
    let mut clusters = find_clusters(&square,4).unwrap();
    remove_subclusters_of(&mut clusters, &Cluster::from(vec![0,1]));

    let one_clusters = vec![
      //vec![0], vec![1],
      vec![2], vec![3],
    ];

    let idx = 0;
    assert_eq!(clusters[idx].len(), one_clusters.len());
    for v in one_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
    } 

    let two_clusters = vec![
      vec![0,1], //vec![0,2], vec![0,3]
      //vec![1,2], vec![1,3],
      vec![2,3],
    ];

    let idx = 1;
    assert_eq!(clusters[idx].len(), two_clusters.len());
    for v in two_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
    } 

    let three_clusters = vec![
      vec![0,1,2], vec![0,1,3], 
      //vec![0,2,3], 
      // vec![1,2,3],
    ];

    let idx = 2;
    assert_eq!(clusters[idx].len(), three_clusters.len());
    for v in three_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
    } 

    let four_clusters = vec![
      vec![0,1,2,3], 
    ];

    let idx = 3;
    assert_eq!(clusters[idx].len(), four_clusters.len());
    for v in four_clusters.iter(){
      assert_eq!(clusters[idx][v].vertices,*v);
    } 
  }
  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
