use crate::signal::Signal;

pub mod find_clusters;
pub mod adjacency;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct Cluster{
  vertices: Vec::<usize>,
  signal: Option<Signal>,
  auxiliary_signal: Option<Signal>,
}
impl Cluster{
  pub fn from(vertices: Vec::<usize>) -> Self{
    Cluster{
      vertices,
      signal: None,
      auxiliary_signal: None,
    }
  }
}
impl ToString for Cluster{
  fn to_string(&self) -> String {
    // "[1,2,3,4]"
    let cluster = &self.vertices;

    if cluster.is_empty(){
      return "[]".to_string();
    }

    let mut string = format!("[{}", cluster[0] );

    for ii in 1..cluster.len(){
      string = format!("{},{}",string,cluster[ii]);
    } 
    format!("{}]",string)
  }
}
impl Cluster{
  fn len(&self) -> usize{
    self.vertices.len()
  }
  //----------------------------------------------------------------------------
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
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

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
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>