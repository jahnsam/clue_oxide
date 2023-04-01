use crate::symmetric_list_2d::SymList2D;
use crate::space_3d::{SymmetricTensor3D,Vector3D};

use crate::clue_errors::CluEError;
use crate::structure::Structure;
use crate::config::Config;

use crate::physical_constants::{HBAR, JOULES_TO_HERTZ,MU0,PI};

//use std::fs;
//use std::fs::File;
//use std::io::{Error, Write};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct HamiltonianTensors{
  pub spin1_tensors: Spin1Tensors, // O(S)
  pub spin2_tensors: Spin2Tensors, // O(S^2)

}
impl HamiltonianTensors{
  //----------------------------------------------------------------------------
  pub fn len(&self) -> usize{
    self.spin1_tensors.len()
  }
  //----------------------------------------------------------------------------
  pub fn generate(structure: &Structure, config: &Config) 
    -> Result<Self,CluEError>
  {
    let n_spins = structure.number_active() + 1;
    let mut spin1_tensors = Spin1Tensors::new(n_spins);
    let mut spin2_tensors = Spin2Tensors::new(n_spins);

    let Some(magnetic_field) = &config.magnetic_field else{
      return Err(CluEError::NoMagneticField);
    };

    let Some(detected_particle) = &structure.detected_particle else {
      return Err(CluEError::NoCentralSpin);
    };

    let eye = SymmetricTensor3D::eye();
    let gamma_e = detected_particle.isotope.gyromagnetic_ratio();


    spin1_tensors.set(0,
        construct_zeeman_tensor(&(gamma_e*&eye),magnetic_field));

    let mut idx0 = 0;
    for particle0 in structure.bath_particles.iter(){

      if !particle0.active {continue;}
      idx0 += 1;

      let gamma0 = particle0.isotope.gyromagnetic_ratio();

      // nuclear Zeeman
      spin1_tensors.set(idx0, 
          construct_zeeman_tensor(&(gamma0*&eye),magnetic_field));

      // hyperfine
      let mut hf_ten = SymmetricTensor3D::zeros();
      for ii in 0..detected_particle.weighted_coordinates.len(){
        let delta_r = &detected_particle.weighted_coordinates.xyz(ii) 
          - &particle0.coordinates;
        let w = detected_particle.weighted_coordinates.weight(ii);
        
        let ten = construct_point_dipole_dipole_tensor(gamma_e,gamma0,&delta_r);
        
        hf_ten = &hf_ten + &(w*&ten);
      }
      spin2_tensors.set(0,idx0, hf_ten);


      // nucleus-nucleus dipole-dipole
      let mut idx1 = 0;
      for particle1 in structure.bath_particles.iter(){
        if !particle1.active {continue;}
         idx1 += 1;
         if idx1 == idx0 {break;}

         let gamma1 = particle1.isotope.gyromagnetic_ratio();

         let delta_r01 = &particle0.coordinates - &particle1.coordinates;

         let dd_ten = construct_point_dipole_dipole_tensor(gamma1,gamma0,
              &delta_r01);

         spin2_tensors.set(idx1,idx0, dd_ten);
      }


    }

    Ok(HamiltonianTensors{
      spin1_tensors,
      spin2_tensors
      })

  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct Spin1Tensors{
  tensors: Vec::< Option<Vector3D> >,
}

impl<'a> Spin1Tensors{

  pub fn new(number: usize) -> Spin1Tensors {
    
    let mut tensors 
     = Vec::<Option<Vector3D>>::with_capacity(number);
    for _ii in 0..(number as usize) {
      tensors.push(None);
    }

    Spin1Tensors{
      tensors,
    }
  }
  //----------------------------------------------------------------------------
  pub fn set(&mut self, n: usize, ten: Vector3D) {
    self.tensors[n] = Some(ten);
  }
  //----------------------------------------------------------------------------
  pub fn get(&'a self, n: usize) -> Option< &'a Vector3D> {
    match &self.tensors[n] {
      Some(ten) => Some(ten),
      None => None,  
    }
  }
  //----------------------------------------------------------------------------
  pub fn remove(&mut self, n: usize){
    self.tensors[n] = None;
  }
  //----------------------------------------------------------------------------
  pub fn len(&self) -> usize{
    self.tensors.len()
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct Spin2Tensors{
  tensors: SymList2D::<SymmetricTensor3D>,
}

impl<'a> Spin2Tensors{

  pub fn new(number: usize) -> Self {
    Spin2Tensors{
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

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn construct_zeeman_tensor(
    gyromagnetic_ratio: &SymmetricTensor3D, 
    magnetic_field: &Vector3D) -> Vector3D
{

  let zeeman = magnetic_field * gyromagnetic_ratio;
  (-JOULES_TO_HERTZ) * &zeeman
}
/*
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

*/
pub fn get_perpendicular_dipole_dipole_frequency(
    gyromagnetic_ratio_1: f64,
    gyromagnetic_ratio_2: f64,
    r: f64
    ) -> f64 {
  JOULES_TO_HERTZ*
  MU0/(4.0*PI)*gyromagnetic_ratio_1*gyromagnetic_ratio_2   
  *HBAR*HBAR/r/r/r
}

pub fn construct_point_dipole_dipole_tensor(
    gyromagnetic_ratio_1: f64,
    gyromagnetic_ratio_2: f64,
    delta_r: &Vector3D
    ) -> SymmetricTensor3D {

  let r = delta_r.norm();
  
  let n = (1.0/r) * delta_r;
 
  let n3nt = 3.0 * &n.self_outer(); 

  let h_perp = get_perpendicular_dipole_dipole_frequency(gyromagnetic_ratio_1,
      gyromagnetic_ratio_2,r);

  let ten = &SymmetricTensor3D::eye() - &n3nt; 

  h_perp * &ten
  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  //use crate::phys::*;

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  #[allow(non_snake_case)]
  fn test_Spin1Tensors() {
    let number = 3;
    let mut sf_tens = Spin1Tensors::new(number);

    for ii in 0..number {

      let mut ten = Vector3D::zeros();
      let val = 1.0 + (ii as f64);
      ten.set_x(val);

      sf_tens.set(ii,ten);

    }

    let expected_values = vec![ 1.0, 2.0, 3.0];

    for ii in 0..(number as usize) {

      let ten = sf_tens.get(ii).unwrap();
      let val = ten.x();
      assert!( (val-expected_values[ii]).abs() < 1e-12 ); 

    }


  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  #[allow(non_snake_case)]
  fn test_Spin2Tensors() {
    let number = 3;
    let mut ss_tens = Spin2Tensors::new(number);

    let mut val = 0.0;
    for ii in 0..(number as usize) {
      for jj in 0..(number as usize) {

        let ten = SymmetricTensor3D::eye();
        val += 1.0;
        let ten = ten.scale(val);

        ss_tens.set(ii,jj,ten);
      }
    }
    let expected_values = vec![ 
      1.0, 4.0, 7.0, 
      4.0, 5.0, 8.0,
      7.0, 8.0, 9.0];

    let mut idx = 0;
    for ii in 0..number {
      for jj in 0..number {
        let ten = ss_tens.get(ii,jj).unwrap();
        let val = ten.xx();
        assert!( (val-expected_values[idx]).abs() < 1e-12 ); 
        idx += 1;
      }
    }
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}
