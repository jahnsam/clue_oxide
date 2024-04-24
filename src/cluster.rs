use crate::clue_errors::*;
use crate::signal::Signal;
use crate::structure::Structure;

pub mod adjacency;
pub mod build_adjacency_list;
pub mod connected_subgraphs;
pub mod find_clusters;
pub mod read_clusters;
pub mod get_subclusters;
pub mod methyl_clusters;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub struct Cluster{
  vertices: Vec::<usize>,
  pub signal: Result<Option<Signal>,CluEError>,
}

//------------------------------------------------------------------------------
impl Cluster{
  pub fn from(vertices: Vec::<usize>) -> Self{
    Cluster{
      vertices,
      signal: Ok(None),
    }
  }
}
//------------------------------------------------------------------------------
impl ToString for Cluster{
  fn to_string(&self) -> String {
    // "[1,2,3,4]"
    let cluster = &self.vertices;

    if cluster.is_empty(){
      return "[]".to_string();
    }

    let mut string = format!("[{}", cluster[0]);

    for vertex in cluster.iter().skip(1){
      string = format!("{},{}",string,vertex);
    } 
    format!("{}]",string)
  }
}
//------------------------------------------------------------------------------
impl Cluster{
  fn len(&self) -> usize{
    self.vertices.len()
  }
  //----------------------------------------------------------------------------
  fn to_string_result(&self,structure: &Structure) -> Result<String,CluEError>{
    // "[1,2,3,4]"
    let cluster = &self.vertices;

    if cluster.is_empty(){
      return Ok("[]".to_string());
    }

    let mut string = format!("[{}", 
        structure.get_reference_index_of_nth_active(cluster[0])? 
        );

    for &vertex in cluster.iter().skip(1){
      string = format!("{},{}",string,
          structure.get_reference_index_of_nth_active(vertex)?
          );
    } 
    Ok(format!("{}]",string))
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  /*
  fn contains(&self, subcluster: &Cluster) -> bool{

    for n in subcluster.vertices.iter(){
      let mut vertex = Vec::<usize>::with_capacity(1); 
      vertex.push(*n);
      let cluster = Cluster::from(vertex);
      if !self.overlaps(&cluster){
        return false;
      }
    }
    true
  }

  fn overlaps(&self, other_cluster: &Cluster) -> bool{

    for m in self.vertices.iter(){
      for n in other_cluster.vertices.iter(){
        if *m == *n{
          return true;
        }
      }
    }
    false
  }
  */
  //----------------------------------------------------------------------------
  /*
  write(filename: &str) -> Result<(),CluEError>{
    let Ok(file) = File::create(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    }

  }
  */
  //----------------------------------------------------------------------------
  pub fn vertices(&self) -> &Vec::<usize>{
    &self.vertices
  } 
  //----------------------------------------------------------------------------
  pub fn to_header(&self, structure: &Structure) -> Result<String,CluEError> {
    // "clu_1_2_3_4"
    let cluster = &self.vertices;

    if cluster.is_empty(){
      return Ok("clu".to_string());
    }

    /*
    let mut string = format!("clu_{}", cluster[0]);

    for vertex in cluster.iter().skip(1){
      string = format!("{}_{}",string,vertex);
    } 
    */
    let mut string = format!("clu_{}", 
        structure.get_reference_index_of_nth_active(cluster[0])? 
        );

    for &vertex in cluster.iter().skip(1){
      string = format!("{}_{}",string,
          structure.get_reference_index_of_nth_active(vertex)?
          );
    } 

    Ok(string)
  }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  //use super::*;

  /*
  #[allow(non_snake_case)]
  #[test]
  fn test_Cluster(){
    let cluster0 = Cluster::from(vec![0,1,2]);
    let cluster1 = Cluster::from(vec![1,2]);
    let cluster2 = Cluster::from(vec![1,2,3]);
    let cluster3 = Cluster::from(vec![3,4]);

    assert!(cluster0.overlaps(&cluster1));
    assert!(cluster0.overlaps(&cluster2));
    assert!(!cluster0.overlaps(&cluster3));

    assert!(cluster0.contains(&cluster1));
    assert!(!cluster0.contains(&cluster2));
    assert!(!cluster0.contains(&cluster3));
  }
  */
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
