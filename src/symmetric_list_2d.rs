#[derive(Debug, Clone)] 
pub struct SymList2D<T>{
  dim: usize,
  elements: Vec::<Option<T>>,
}

impl<T> SymList2D<T>{

  pub fn new(dim: usize) -> Self{
    let n_elem = (dim*dim + dim)/2;
    let mut elements = Vec::<Option<T>>::with_capacity(n_elem);
   
    for _ii in 0..n_elem{
      elements.push(None);
    }

    SymList2D{dim, elements}
  }
  //----------------------------------------------------------------------------
  pub fn dim(&self) -> usize {self.dim}
  //----------------------------------------------------------------------------
  pub fn get(& self, row: usize, col: usize) -> Option<&T>{
    let idx = get_uppertriangle_index(row,col,self.dim);
    
    match &self.elements[idx] {
      Some(value) => Some(value),
      None => None,
    }
  }
  //----------------------------------------------------------------------------
  pub fn set(&mut self, row: usize, col: usize, value: T){
    let idx = get_uppertriangle_index(row,col,self.dim);
    self.elements[idx] = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn remove(&mut self, row: usize, col: usize){
    let idx = get_uppertriangle_index(row,col,self.dim);
    self.elements[idx] = None;
  }


}

fn get_uppertriangle_index(row: usize, col: usize, dim: usize) -> usize {
  
  let mut row = row as i32;
  let mut col = col as i32;
  let dim = dim as i32;
  
  if row>col {
    std::mem::swap(&mut row, &mut col)
  }

  let idx = (-row*row +(2*dim - 1)*row + 2*col)/2;

  idx as usize
}


#[cfg(test)]
mod tests{
  use super::*;

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  fn test_get_uppertriangle_index() {

    let dim = 3;
    // 0 1 2
    // - 3 4
    // - - 5
    let expected_values = vec![ 0,1,2, 1,3,4, 2,4,5 ]; 
    let mut counter = 0; 
    for row in 0..dim {
      for col in 0..dim {  
        let idx = get_uppertriangle_index(row, col, dim);
        let ev = expected_values[counter];
        assert!(ev==idx);
        counter += 1;
      }
    }  
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}
