use crate::symmetric_list_2d::SymList2D;
use crate::space_3d::SymmetricTensor3D;

//use super::phys;

//use std::fs;
//use std::fs::File;
//use std::io::{Error, Write};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct HamiltonianTensors{
  spin_field: SpinFieldTensors,
  spin_spin: SpinSpinTensors,

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct SpinFieldTensors{
  tensors: Vec::< Option<SymmetricTensor3D> >,
}

impl<'a> SpinFieldTensors{

  pub fn new(number: usize) -> SpinFieldTensors {
    
    let mut tensors 
     = Vec::<Option<SymmetricTensor3D>>::with_capacity(number);
    for _ii in 0..(number as usize) {
      tensors.push(None);
    }

    SpinFieldTensors{
      tensors,
    }
  }
  //----------------------------------------------------------------------------
  pub fn set(&mut self, n: usize, ten: SymmetricTensor3D) {
    self.tensors[n] = Some(ten);
  }
  //----------------------------------------------------------------------------
  pub fn get(&'a self, n: usize) -> Option< &'a SymmetricTensor3D> {
    match &self.tensors[n] {
      Some(ten) => Some(ten),
      None => None,  
    }
  }
  //----------------------------------------------------------------------------
  pub fn remove(&mut self, n: usize){
    self.tensors[n] = None;
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct SpinSpinTensors{
  tensors: SymList2D::<SymmetricTensor3D>,
}

impl<'a> SpinSpinTensors{

  pub fn new(number: usize) -> Self {
    SpinSpinTensors{
      tensors: SymList2D::new(number),
    }
  }
  //----------------------------------------------------------------------------
  pub fn get(&'a self, m: usize, n: usize) -> Option<&'a SymmetricTensor3D >{
    self.tensors.get(m,n)
  }
  //----------------------------------------------------------------------------
  pub fn set(&mut self, m: usize, n: usize, ten: SymmetricTensor3D){
    self.tensors.set(m,n,ten);
  }
  //----------------------------------------------------------------------------
  pub fn remove(&mut self, m: usize, n: usize){
    self.tensors.remove(m,n);
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

/*
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn construct_zeeman_tensor(
    gyromagnetic_ratio: &mox::Mat, 
    magnetic_field: &mox::Mat) -> Result<mox::Mat,String> {

  let magnetic_field = magnetic_field.trans();
  let mut zeeman = magnetic_field.multiply( &gyromagnetic_ratio )?;
  zeeman.scalar_multiply_mut( -phys::RADS_TO_HZ );
  Ok(zeeman)
}

pub fn construct_hyperfine_tensor(
    azz: f64, fc: f64, rot_hf2lab: &mox::Mat) -> Result<mox::Mat,String> {

    if rot_hf2lab.n_rows() != 3 || rot_hf2lab.n_cols() != 3 {
      let err_str = format!("hf to lab rotation matrix must be 3x3, not {}x{}",
          rot_hf2lab.n_rows(), rot_hf2lab.n_cols());
      return Err(err_str);
    }

    let ten_hf = mox::Mat::from(3,3, vec![-azz/2.0 +fc,      0.0, 0.0,
                                          0.0, -azz/2.0 +fc, 0.0,
                                          0.0,          0.0, azz +fc ]);


   rot_hf2lab.multiply( &ten_hf.multiply(&rot_hf2lab.trans())? )
}

pub fn get_perpendicular_dipole_dipole_frequency(
    gyromagnetic_ratio_1: f64,
    gyromagnetic_ratio_2: f64,
    r: f64
    ) -> f64 {
  phys::J_TO_HZ*
  phys::MU0/(4.0*phys::PI)*gyromagnetic_ratio_1*gyromagnetic_ratio_2   
  *phys::HBAR*phys::HBAR/r/r/r
}

pub fn construct_point_dipole_dipole_tensor(
    gyromagnetic_ratio_1: f64,
    gyromagnetic_ratio_2: f64,
    delta_r: &mox::Mat
    ) -> Result<mox::Mat, String> {

  if delta_r.n_rows() != 3 || delta_r.n_cols() != 1 {
    let err_str = format!("delta_r must be 3x1 not {}x{}", 
        delta_r.n_rows(),delta_r.n_cols());
    return Err(err_str);
  } 
  let r = delta_r.l2_norm()?;
  
  let n3nt = delta_r.multiply( &delta_r.trans() )?.scalar_multiply(3.0/r/r); 

  let h_perp = get_perpendicular_dipole_dipole_frequency(gyromagnetic_ratio_1,
      gyromagnetic_ratio_2,r);

  let ten = &mox::Mat::eye(3,3) - &n3nt; 

  Ok(ten)
  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

/*
#[cfg(test)]
mod tests{
  use super::*;
  //use crate::phys::*;

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
*/
