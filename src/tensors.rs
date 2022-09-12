use matrix_oxide as mox;
//use super::phys::*;

fn get_uppertriangle_index(row: usize, col: usize, dim: u32) -> usize {
  
  let mut row = row as i32;
  let mut col = col as i32;
  let dim = dim as i32;
  
  if row>col {
    let row0 = row;
    row = col;
    col = row0;
  }

  let idx = (-row*row +(2*dim - 1)*row + 2*col)/2;

  return idx as usize;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct SpinFieldTensors{
  tensors: Vec::< Option<mox::Mat> >,
}

impl<'a> SpinFieldTensors{

  pub fn new(number: u32) -> SpinFieldTensors {
    
    let mut tensors = Vec::<Option<mox::Mat>>::with_capacity(number as usize);
    for _ii in 0..(number as usize) {
      tensors.push(None);
    }

    SpinFieldTensors{
      tensors,
    }
  }

  pub fn set(&mut self, n: usize, ten: mox::Mat) -> bool {
    if n >= self.tensors.len() {
      return false;
    }

    self.tensors[n] = Some(ten);
    true
  }

  pub fn get(&'a self, n: usize) -> Option< &'a mox::Mat> {
    if n >= self.tensors.len() {
      return None;
    }

    match &self.tensors[n] {
      Some(ten) => Some(ten),
      None => None,  
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct SpinSpinTensors{
  number: u32,
  tensors: Vec::< Option<mox::Mat> >,
}

impl<'a> SpinSpinTensors{

  pub fn new(number: u32) -> SpinSpinTensors {
    let n_elems = ( number*(number + 1)/2 ) as usize;
    let mut tensors = Vec::< Option<mox::Mat> >::with_capacity(n_elems);

    for _ii in 0..n_elems {
      tensors.push(None);
    }

    SpinSpinTensors{
      number,
      tensors,
    }
  }


  pub fn get(&'a self, m: usize, n: usize) -> Option<&'a mox::Mat >{

    let idx = get_uppertriangle_index(m, n, self.number);

    if idx >= self.tensors.len() {
      return None;
    }


    match &self.tensors[idx] {
      Some(ten) => Some(ten),
      None => None,  
    }
  }

  pub fn set(&mut self, m: usize, n: usize, ten: mox::Mat) -> bool {

    let idx = get_uppertriangle_index(m, n, self.number);

    if idx >= self.tensors.len() {
      return false;
    }

    self.tensors[idx] = Some(ten);
    true
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#[cfg(test)]
mod tests{
  use super::*;
  use crate::phys::*;

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  #[allow(non_snake_case)]
  fn test_SpinFieldTensors() {
    let number: u32 = 3;
    let mut sf_tens = SpinFieldTensors::new(number);

    for ii in 0..(number as usize) {

      let ten = mox::Mat::ones(1,1);
      let val = 1.0 + (ii as f64);
      let ten = ten.scalar_multiply(val);

      assert!(sf_tens.set(ii,ten));

    }

    let expected_values = vec![ 1.0, 2.0, 3.0];

    for ii in 0..(number as usize) {

      let ten = sf_tens.get(ii).unwrap();
      let val = ten.get(0,0).unwrap();
      assert!( (val-expected_values[ii]).abs() < ERROR_THRESHOLD ); 

    }


  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


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
        let idx = get_uppertriangle_index(row, col, dim as u32);
        let ev = expected_values[counter];
        assert!(ev==idx);
        counter += 1;
      }
    }  
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  #[allow(non_snake_case)]
  fn test_SpinSpinTensors() {
    let number: u32 = 3;
    let mut ss_tens = SpinSpinTensors::new(number);

    let mut val = 0.0;
    for ii in 0..(number as usize) {
      for jj in 0..(number as usize) {

        let ten = mox::Mat::ones(1,1);
        val += 1.0;
        let ten = ten.scalar_multiply(val);

        assert!(ss_tens.set(ii,jj,ten));
      }
    }
    let expected_values = vec![ 
      1.0, 4.0, 7.0, 
      4.0, 5.0, 8.0,
      7.0, 8.0, 9.0];

    let mut idx = 0;
    for ii in 0..(number as usize) {
      for jj in 0..(number as usize) {
        let ten = ss_tens.get(ii,jj).unwrap();
        let val = ten.get(0,0).unwrap();
        assert!( (val-expected_values[idx]).abs() < ERROR_THRESHOLD ); 
        idx += 1;
      }
    }
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}
