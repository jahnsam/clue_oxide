use crate::math;


#[derive(Debug,Clone)]
pub struct AdjacencyList{
  list: Vec::< Option<Vec::<usize>> >,
}
impl AdjacencyList{

  pub fn len(&self) -> usize{
    self.list.len()
  }

  pub fn with_capacity(n: usize) -> Self{
    let mut list = Vec::< Option<Vec::<usize>> >::with_capacity(n);
    for _ii in 0..n{
      list.push(None);
    }

    AdjacencyList{list}
  }
  //----------------------------------------------------------------------------

  pub fn are_connected(&self, m: usize, n: usize) -> bool{
    self.contains(m,n)
  }

  //----------------------------------------------------------------------------
  fn contains(&self,m: usize, n: usize) -> bool{

    if let Some(neighbors) = &self.list[m]{
      for ii in neighbors.iter(){
        if *ii == n{
          return true;
        }
      }
    }

    false
  }
  //----------------------------------------------------------------------------
  pub fn connect(&mut self, m: usize, n: usize){
    if m==n {
      self.activate(m);
      return
    }

    self.insert(m,n);
    self.insert(n,m);
  }

  //----------------------------------------------------------------------------
  fn activate(&mut self, n: usize){
    if self.list[n] == None{
      self.list[n] = Some(Vec::new());
    }
  }
  //----------------------------------------------------------------------------
  fn insert(&mut self, m: usize, n: usize){

    if let Some(neighbors) = &mut self.list[m]{
      neighbors.push(n);
      *neighbors = math::unique((*neighbors).clone());
      neighbors.sort();
    }else{
      self.list[m] = Some(Vec::from([n]));
    }
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbors<'a>(&'a self, n: usize) -> Option<&'a Vec::<usize>>{
    if let Some(neighbors) = &self.list[n]{
      return Some(neighbors);
    }
    None
  }
  //----------------------------------------------------------------------------
  pub fn n_active_vertices(&self) -> usize{
    let mut counter = 0;
    for ii in 0..self.len(){
      if let Some(_) = self.list[ii]{
        counter += 1;
      }
    }
    counter
  }
  //----------------------------------------------------------------------------
  pub fn get_active_vertices(&self) -> Vec::<usize>{
    let counter = self.n_active_vertices();

    let mut vertices = Vec::<usize>::with_capacity(counter);

    for ii in 0..self.len(){
      if let Some(_) = self.list[ii]{
        vertices.push(ii);
      }
    }

    vertices
  }
  //----------------------------------------------------------------------------

}

#[cfg(test)]
mod tests {
  use super::*;

  #[allow(non_snake_case)]
  #[test]
  fn test_AdjacencyList(){
    let mut adjacency_list = AdjacencyList::with_capacity(6);
    
    assert!(!adjacency_list.are_connected(0,1));
    
    let neighbors = adjacency_list.get_neighbors(0);
    assert_eq!(neighbors,None);

    adjacency_list.connect(0,1);
    assert!(adjacency_list.are_connected(0,1));
    assert!(adjacency_list.are_connected(1,0));

    assert!(!adjacency_list.are_connected(0,2));

    let neighbors = adjacency_list.get_neighbors(0).unwrap();
    assert_eq!(*neighbors,vec![1]);

    adjacency_list.connect(0,2);
    let neighbors = adjacency_list.get_neighbors(0).unwrap();
    assert_eq!(*neighbors,vec![1,2]);

    adjacency_list.connect(0,2);
    let neighbors = adjacency_list.get_neighbors(0).unwrap();
    assert_eq!(*neighbors,vec![1,2]);

    let neighbors = adjacency_list.get_neighbors(1).unwrap();
    assert_eq!(*neighbors,vec![0]);

    adjacency_list.connect(4,2);
    adjacency_list.connect(4,0);
    let neighbors = adjacency_list.get_neighbors(4).unwrap();
    assert_eq!(*neighbors,vec![0,2]);

    let vertices = adjacency_list.get_active_vertices();
    assert_eq!(vertices,vec![0,1,2,4]);
  }
}



