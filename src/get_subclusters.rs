pub fn build_subclusters(cluster: &[usize]) -> Vec::<Vec::<usize>>{

  let cluster_size = cluster.len();
  
  if cluster_size <= 1 {
    // There are no proper subclusters.
    return Vec::<Vec::<usize>>::new();
  }
  
  let n_subclusters = 2_usize.pow( cluster_size as u32);
  
  let mut subclusters = Vec::<Vec::<usize>>::with_capacity(n_subclusters-2);
  for _ii in 0..n_subclusters-2{
    subclusters.push(Vec::<usize>::with_capacity(cluster.len()));
  }

  for (clu_idx, vertex) in cluster.iter().enumerate(){
    let subcluster_selection 
      = include_vertex_in_which_subcluster(
          clu_idx, cluster_size, n_subclusters);

    for idx in subcluster_selection.iter(){
      subclusters[*idx].push(*vertex);
    }
    
  }

  subclusters
}
//------------------------------------------------------------------------------
pub fn include_vertex_in_which_subcluster(
    idx: usize,cluster_size: usize, n_subclusters: usize) 
  -> Vec::<usize>
{
  let n: usize = 2_usize.pow(cluster_size as u32);

  let increment: usize = 2_usize.pow(idx as u32);

  let mut do_push = false;
  let mut indices = Vec::<usize>::with_capacity(n);
  for ii in 0..n_subclusters{
    if ii % increment == 0 { do_push = !do_push;} 
    if ii == 0 {continue;}
    if do_push {
      indices.push(ii-1);
    }
  }

  indices
}

#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_include_vertex_in_which_subcluster(){
    let expected = vec![
      vec![1,3,5],
      vec![0,3,4],
      vec![0,1,2],
    ];
    let cluster_size = 3;
    let n_subclusters = 2_usize.pow( cluster_size as u32);
    for ii in 0..cluster_size{
      let indices 
        = include_vertex_in_which_subcluster(ii,cluster_size,n_subclusters);
      assert_eq!(indices,expected[ii]);
    }
  }

  #[test]
  fn test_build_subclusters(){
    let subclusters = build_subclusters(&[1,2,3]);
    assert_eq!(subclusters,vec![
        vec![2, 3], vec![1, 3], vec![3], vec![1, 2], vec![2], vec![1]] );

    let subclusters = build_subclusters(&[1]);
    assert!(subclusters.is_empty());

    let subclusters = build_subclusters(&[]);
    assert!(subclusters.is_empty());
  }
}
