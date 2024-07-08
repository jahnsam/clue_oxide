use crate::clue_errors::CluEError;
use crate::cluster::Cluster;
use crate::cluster::find_clusters::ClusterSet;
use crate::structure::Structure;

use std::io::BufRead;
use substring::Substring;

//------------------------------------------------------------------------------
/// This function reads a cluster file and build the corresponding `ClusterSet`.
pub fn read_cluster_file(filename: &str, structure: &Structure) 
  -> Result<ClusterSet,CluEError>
{

  let mut clusters = initialize_clusters(filename)?;

  let Ok(file) = std::fs::File::open(filename) else{
    return Err(CluEError::CannotOpenFile(filename.to_string()));
  };

  let lines = std::io::BufReader::new(file).lines();

  for line_result in lines {
    let Ok(line) = line_result else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };

    if line.contains("#") || !line.contains("["){
      continue;
    }
    
    let cluster_reference_indices = parse_cluster_line(&line)?;

    let mut cluster_indices 
      = Vec::<usize>::with_capacity(cluster_reference_indices.len());

      for &ref_idx in cluster_reference_indices.iter(){
        let n = structure.get_nth_active_from_reference_index(ref_idx)?;
        cluster_indices.push(n);
      }

    let cluster = Cluster::from(cluster_indices);

    let size_idx = cluster.len() - 1;
    clusters[size_idx].push(cluster)

  }

  Ok(ClusterSet::from(clusters))
}
//------------------------------------------------------------------------------
// This function reads lines that specify single clusters.
fn parse_cluster_line(line: &str) 
  -> Result<Vec::<usize>,CluEError>
{

  let Some(open_bracket) = line.find("[") else{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };
  let Some(close_bracket) = line.find("]") else{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };

  if open_bracket + 1 >= close_bracket {
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  }

  let cluster_str = line.substring(open_bracket+1,close_bracket);
  let n: usize = cluster_str.split(",").map(|_| 1).sum();

  let mut cluster = Vec::<usize>::with_capacity(n);
  for s in cluster_str.split(","){
    let Ok(index) = s.parse::<usize>() else{
      return Err(CluEError::ClusterLineFormatError(line.to_string()));
    };
    cluster.push(index);
  }

  Ok(cluster)
}
//------------------------------------------------------------------------------
// This function uses the meta data in the cluster files to initialize 
// the clusters' `Vec::<Vec::<Cluster>` to the correct size.
fn initialize_clusters(filename: &str) 
  -> Result<Vec::<Vec::<Cluster>>,CluEError>
{
  let Ok(file) = std::fs::File::open(filename) else{
    return Err(CluEError::CannotOpenFile(filename.to_string()));
  };

  let lines = std::io::BufReader::new(file).lines();

  for line_result in lines {
    let Ok(line) = line_result else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };

    if line.contains("#[clusters"){
      let number_clusters = parse_number_clusters(&line, filename)?;
      let mut clusters = Vec::<Vec::<Cluster>>::with_capacity(
          number_clusters.len());

      for &n in number_clusters.iter(){
        clusters.push(Vec::<Cluster>::with_capacity(n));
      }

      return Ok(clusters);
    }

  }

  Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()))
}
//------------------------------------------------------------------------------
// This function reads the meta data in the cluster files to learn how many 
// clusters of each size there are.
fn parse_number_clusters(line: &str, filename: &str) 
  -> Result<Vec::<usize>,CluEError>
{
  
  let new_line = remove_char(line,' ');
  
  let Some(idx) = new_line.find("number_clusters") else{
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  };

  let new_line = new_line.substring(idx,new_line.len());

  let Some(open_bracket) = new_line.find("[") else{
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  };
  let Some(close_bracket) = new_line.find("]") else{
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  };


  if open_bracket + 1 >= close_bracket {
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  }
  let new_line = new_line.substring(open_bracket + 1,close_bracket) ;

  let max_cluster_size: usize = new_line.split(",").map(|_| 1).sum();
  let mut n_clusters = Vec::<usize>::with_capacity(max_cluster_size);

  for s in new_line.split(","){
    let Ok(n) = s.parse::<usize>() else {
      return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
    };
    n_clusters.push(n);
  }

  Ok(n_clusters)
}
//------------------------------------------------------------------------------
// This function takes a `&str` and builds a `String` identical except 
// all instances of a specified `char` removed.
fn remove_char(s: &str, c: char) -> String {
  let mut out_str = s.trim().to_owned();
  out_str.retain(|ch| {ch != c });
  out_str
}
//------------------------------------------------------------------------------


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_read_cluster_file(){
    let filename = "assets/cluster_file.txt".to_string();
    let value = read_cluster_file(&filename).unwrap();

    let clusters = vec![
      vec![
        Cluster::from(vec![1]),
        Cluster::from(vec![2]),
        Cluster::from(vec![4]),
        Cluster::from(vec![8]),
      ],
      vec![
        Cluster::from(vec![1,2]),
        Cluster::from(vec![2,4]),
        Cluster::from(vec![4,8]),
      ],
      vec![
        Cluster::from(vec![2,4,8]),
      ],
    ];
    let expected = ClusterSet::from(clusters);

    assert_eq!(value,expected);
  }
  */
  //----------------------------------------------------------------------------
  #[test]
  fn test_initialize_clusters(){
    let filename = "assets/cluster_file.txt".to_string();
    let value = initialize_clusters(&filename).unwrap();

    assert_eq!(value.len(),3);
  
    let expected_capacities = vec![4,3,1];
    
    for (ii ,&c) in expected_capacities.iter().enumerate(){
      assert!(value[ii].is_empty());
      assert_eq!(value[ii].capacity(),c);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_cluster_line(){
    let line = "[7,13,23,42]";
    let value = parse_cluster_line(&line).unwrap();
    assert_eq!(value,vec![7,13,23,42]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_number_clusters(){
    let line = "#[clusters, number_clusters = [4276,32395] ]".to_string();
    let value = parse_number_clusters(&line, &"<file>").unwrap();
    assert_eq!(value, vec![4276,32395]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_remove_char(){
    let s = "#[clusters, number_clusters = [4276,32395] ]".to_string();
    let value = remove_char(&s,' ');
    let expected = "#[clusters,number_clusters=[4276,32395]]".to_string();
    assert_eq!(value,expected);
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
