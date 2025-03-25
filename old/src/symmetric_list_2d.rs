/// `SymList2D` is a data structure for a 2D list where `list[i][j]==list[j][i]`
/// for all `i` and `j`.
/// `SymList2D` starts out with all entries set to `None`.
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
  /// This function returns the number of unique indices for `SymList2D` as
  /// a 2D list, 
  /// equivalently the dimension of `SymList2D` viewed as a matrix.
  pub fn dim(&self) -> usize {self.dim}
  //----------------------------------------------------------------------------
  /// Viewed as a symmetric matrix, this functions gets the `(row,col)` entry.
  pub fn get(& self, row: usize, col: usize) -> Option<&T>{
    let idx = get_uppertriangle_index(row,col,self.dim);
    
    match &self.elements[idx] {
      Some(value) => Some(value),
      None => None,
    }
  }
  //----------------------------------------------------------------------------
  /// Viewed as a symmetric matrix, this functions sets the `(row,col)` entry.
  pub fn set(&mut self, row: usize, col: usize, value: T){
    let idx = get_uppertriangle_index(row,col,self.dim);
    self.elements[idx] = Some(value);
  }
  //----------------------------------------------------------------------------
  /// This functions sets the `(row,col)` entry to `None`.
  pub fn remove(&mut self, row: usize, col: usize){
    let idx = get_uppertriangle_index(row,col,self.dim);
    self.elements[idx] = None;
  }


}

// This function takes to indices refering to a row and a column and
// converts it to a 1D index. Below is an example for the 3-by-3 case.
// 1 2 3
// 2 4 5
// 3 5 6
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
