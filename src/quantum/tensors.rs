
use crate::clue_errors::CluEError;
use crate::config::Config;
use crate::config::particle_config::TensorSpecifier;
use crate::physical_constants::{HBAR, JOULES_TO_HERTZ,MU0,PI};
use crate::space_3d::{SymmetricTensor3D,Vector3D};
use crate::structure::Structure;
use crate::symmetric_list_2d::SymList2D;


//use std::fs;
//use std::fs::File;
//use std::io::{Error, Write};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct HamiltonianTensors{
  pub spin_multiplicities: Vec::<usize>,
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
    let mut spin_multiplicities = Vec::<usize>::with_capacity(n_spins);
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

    spin_multiplicities.push(detected_particle.isotope.spin_multiplicity());


    spin1_tensors.set(0,
        construct_zeeman_tensor(&(gamma_e*&eye),magnetic_field));

    let mut idx0 = 0;
    for particle0 in structure.bath_particles.iter(){

      if !particle0.active {continue;}
      idx0 += 1;

      let gamma0 = particle0.isotope.gyromagnetic_ratio();

      spin_multiplicities.push(particle0.isotope.spin_multiplicity());

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
      spin_multiplicities,
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

  (-0.5/PI)*&(magnetic_field * gyromagnetic_ratio)
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
//------------------------------------------------------------------------------
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

  let ten = h_perp*&(&SymmetricTensor3D::eye() - &n3nt); 

  ten
  
}
//------------------------------------------------------------------------------
/// This function takes three values and three vectors and constructs
/// ___T___ = ___XE___trans(___X___), where ___E___ is the digonal matrix of
/// the values and ___X___ is a 3×3 matrix with each vector as a column.
pub fn construct_symmetric_tensor_from_values_and_vectors(
    evals: &[f64;3], evecs: &[Vector3D;3]) -> SymmetricTensor3D
{

  let mut ten = SymmetricTensor3D::zeros();
  for (ii,v) in evecs.iter().enumerate(){
    ten = &ten + &(v.self_outer()).scale(evals[ii]);
  }

  ten

}
//------------------------------------------------------------------------------
fn construct_symmetric_tensor_from_tensor_specifier(
    tensor_specifier: &TensorSpecifier, particle_index: usize,
    structure: &Structure, config: &Config) 
  -> Result<SymmetricTensor3D, CluEError>
{
  let Some(values) = tensor_specifier.values else{
    return Err(CluEError::NoTensorValues);
  }; 

  let mut axes = Vec::<Vector3D>::with_capacity(2);
  let mut axis_dims = Vec::<usize>::with_capacity(2);

  const X: usize = 0;
  const Y: usize = 1;
  const Z: usize = 2;
  if let Some(axis_specifier) = &tensor_specifier.z_axis{
    let mut axis = axis_specifier.to_vector3d(particle_index,structure,
        config)?;
    axis = axis.normalize();
    axes.push(axis);
    axis_dims.push(Z);
  }

  if let Some(axis_specifier) = &tensor_specifier.x_axis{
    let mut axis = axis_specifier.to_vector3d(particle_index,structure,
        config)?;

    if axes.len() == 1{
      axis = &axis - &axes[0].scale(axes[0].dot(&axis));
    }
    axis = axis.normalize();
    axes.push(axis);
    axis_dims.push(X);
  }

  if axes.len() < 2 {
    if let Some(axis_specifier) = &tensor_specifier.y_axis{
      let mut axis = axis_specifier.to_vector3d(particle_index, structure, 
          config)?;
      if axes.len() == 1{
        axis = &axis - &axes[0].scale(axes[0].dot(&axis));
      }
      axis = axis.normalize();
      axes.push(axis);
      axis_dims.push(Y);
    }
  }


  if axis_dims.len() != 2{
    return Err(CluEError::WrongNumberOfAxes(axes.len(),2))
  }
  let axis_dims = [axis_dims[0],axis_dims[1]];

  let all_3_axes: [Vector3D; 3]; 
  match axis_dims{
    [Z,X] => {
       let mut axis = axes[0].cross(&axes[1]);
       axis = axis.scale(1.0/axis.norm());
       all_3_axes = [axes[1].clone(),axis,axes[0].clone()];
    },
    [Z,Y] => {
       let mut axis = axes[1].cross(&axes[0]);
       axis = axis.scale(1.0/axis.norm());
       all_3_axes = [axis,axes[1].clone(),axes[0].clone()];
    },
    [X,Y] => {
       let mut axis = axes[0].cross(&axes[1]);
       axis = axis.scale(1.0/axis.norm());
       all_3_axes = [axes[0].clone(),axes[1].clone(),axis];
    },
    _ => return Err(CluEError::InvalidAxes),
  }

  let ten = construct_symmetric_tensor_from_values_and_vectors(
      &values, &all_3_axes);

  Ok(ten)

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::physical_constants::*;
  use crate::space_3d::{SymmetricTensor3D, Vector3D};
  use crate::config::particle_config::ParticleConfig;
  use crate::structure::pdb;
  use crate::Config;
  use crate::config::DetectedSpinCoordinates;
  use crate::structure::particle_filter::{ParticleFilter, 
    SecondaryParticleFilter,VectorSpecifier};


  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_symmetric_tensor_from_tensor_specifier(){
    let filename = "./assets/TEMPO.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults();
    structure.build_primary_structure(&config).unwrap();

    config.particles.push( ParticleConfig::new("nitrogen".to_string()) );
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Nitrogen];
    config.particles[0].filter = Some(filter);

    let label = String::from("oxygen");
    config.particles.push( ParticleConfig::new(label.clone()) );
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Oxygen];
    config.particles[1].filter = Some(filter);

    let label = String::from("carbon");
    config.particles.push( ParticleConfig::new(label.clone()) );
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Carbon];
    config.particles[2].filter = Some(filter);

    let vector_specifier_no = VectorSpecifier::Diff(
        SecondaryParticleFilter::Particle,"nitrogen".to_string(),
        SecondaryParticleFilter::Bonded,"oxygen".to_string());

    let vector_specifier_cc = VectorSpecifier::Diff(
        SecondaryParticleFilter::Bonded,"carbon".to_string(),
        SecondaryParticleFilter::Bonded,"carbon".to_string());

    let r_n = Vector3D::from([36.440e-10, 36.900e-10,  37.100e-10]);
    let r_o = Vector3D::from([35.290e-10,  36.430e-10,  37.810e-10]);
    let r_c1 = Vector3D::from([37.700e-10, 36.150e-10, 37.340e-10]);
    let r_c19 = Vector3D::from([36.610e-10, 38.380e-10, 36.850e-10]);

    let delta_r_no = &r_o - &r_n;
    let delta_r_cc = &r_c19 - &r_c1;

    let particle_index = 27;
    assert_eq!(structure.bath_particles[particle_index].element,
        Element::Nitrogen);
    let r_no = vector_specifier_no.to_vector3d(particle_index,&structure,
        &config).unwrap();
    assert_eq!(r_no,delta_r_no);

    let r_cc = vector_specifier_cc.to_vector3d(particle_index,&structure,
        &config).unwrap();
    assert_eq!(r_cc,delta_r_cc);


    // de Oliveira, M.; Knitsch, R.; Sajid, M.; Stute, A.; Elmer, L.-M.; 
    // Kehr, G.; Erker, G.; Magon, C. J.; Jeschke, G.; Eckert, H.
    // Aminoxyl Radicals of B/P Frustrated Lewis Pairs:
    // Refinement of the Spin-Hamiltonian Parameters by Field- and
    // Temperature-Dependent Pulsed EPR Spectroscopy.
    // PLoS ONE 2016, 11 (6), e0157944.
    // https://doi.org/10.1371/journal.pone.0157944.
    let e2qQh = 3.5*1e6; // Hz
    let eta = 0.68;
    let values = [e2qQh*(eta-1.0), e2qQh*(-eta-1.0), 2.0*e2qQh];
    assert_eq!(values,[-1119999.9999999998, -5880000.000000001, 7000000.0]);

    let tensor_specifier = TensorSpecifier{
      values: Some(values.clone()),
      x_axis: Some(vector_specifier_no),
      y_axis: Some(vector_specifier_cc),
      z_axis: None,
    };

    let tensor = construct_symmetric_tensor_from_tensor_specifier(
        &tensor_specifier, particle_index,&structure, &config).unwrap();

    let x = r_no.normalize();
    let y = (&r_cc - &x.scale(x.dot(&r_cc))).normalize();

    assert_eq!(&tensor*&r_no,r_no.scale(values[0]));
    assert!((&(&tensor*&x) 
          - &x.scale(values[0])).norm()/values[0].abs() < 1e-9 );

    assert!((&(&tensor*&y) 
          - &y.scale(values[1])).norm()/values[1].abs() < 1e-9 );
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_symmetric_tensor_from_values_and_vectors(){
    let x = Vector3D::from([1.0, 0.0, 0.0]);
    let y = Vector3D::from([0.0, 1.0, 0.0]);
    let z = Vector3D::from([0.0, 0.0, 1.0]);

    let xyz = [x,y,z];
    let evals = [-1.0,2.0,10.0];

    let ten = construct_symmetric_tensor_from_values_and_vectors(
        &evals,&xyz);

    let xyz_ten = SymmetricTensor3D::from([evals[0],0.0,0.0,
                                                    evals[1], 0.0,
                                                              evals[2]]);

    assert_eq!(ten.xx(),xyz_ten.xx());
    assert_eq!(ten.xy(),xyz_ten.xy());
    assert_eq!(ten.xz(),xyz_ten.xz());
    assert_eq!(ten.yy(),xyz_ten.yy());
    assert_eq!(ten.yz(),xyz_ten.yz());
    assert_eq!(ten.zz(),xyz_ten.zz());

    let v1 = Vector3D::from([1.0,1.0,1.0]).scale(1.0/3.0f64.sqrt());
    let v2 = Vector3D::from([-1.0,-1.0,2.0]).scale(1.0/6.0f64.sqrt());
    let v3 = Vector3D::from([1.0,-1.0,0.0]).scale(1.0/2.0f64.sqrt());

    assert!( (v1.norm() - 1.0) < 1e-12);
    assert!( (v2.norm() - 1.0) < 1e-12);
    assert!( (v3.norm() - 1.0) < 1e-12);
    assert!( v1.dot(&v2) < 1e-12);
    assert!( v1.dot(&v3) < 1e-12);
    assert!( v2.dot(&v3) < 1e-12);
    
    let evecs = [v1,v2,v3];


    let ten = construct_symmetric_tensor_from_values_and_vectors(
        &evals,&evecs);

    for (ii,v) in evecs.iter().enumerate(){
      assert_eq!(&ten*v , v.scale(evals[ii]));
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_zeeman_tensor() {
    let gmr = ELECTRON_G*MUB/HBAR;
    assert!((gmr+1.76085963023e11).abs() < 0.00000000053e11);
    let gyromagnetic_ratio = SymmetricTensor3D::from([gmr,0.0, 0.0,
                                                          gmr, 0.0,
                                                               gmr]); 
    let magnetic_field = Vector3D::from([0.0, 0.0, 1.2]);

    let zeeman = construct_zeeman_tensor(&gyromagnetic_ratio,&magnetic_field);
    assert_eq!(zeeman.x(),0.0);
    assert_eq!(zeeman.y(),0.0);
    assert!( (zeeman.z() - 33629941709.0).abs() < 0.1);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_perpendicular_dipole_dipole_frequency(){
    let gmr = Isotope::Hydrogen1.gyromagnetic_ratio();
    let gmr0 = 2.6752218744e8;
    assert!( (gmr-gmr0).abs() < 0.000_000_0011e8);
    let r = 1.5e-10;
    let freq = get_perpendicular_dipole_dipole_frequency(gmr,gmr,r);
    
    let ref_freq = 35591.15890;
    assert!((freq-ref_freq).abs()/ref_freq<1e-9);
  
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_point_dipole_dipole_tensor(){
  
    let gmr = Isotope::Hydrogen1.gyromagnetic_ratio();
    let r = 1.5e-10;
    let freq = get_perpendicular_dipole_dipole_frequency(gmr,gmr,r);

    let phi_list = (0..=10).map(|x| (x as f64)/10.0*PI*2.0)
      .collect::<Vec::<f64>>();

    let theta_list = (0..=10).map(|x| (2.0*(x as f64)/10.0 - 1.0f64).acos())
      .collect::<Vec::<f64>>();

    for &theta in theta_list.iter(){
      let ct = theta.cos();
      let st = theta.sin();

      for &phi in phi_list.iter(){
        let cp = phi.cos();
        let sp = phi.sin();


        let delta_r = r*&Vector3D::from([st*cp,st*sp,ct]);

        assert!((delta_r.norm()-r).abs()/r<1e-12);

        let ten = construct_point_dipole_dipole_tensor(gmr,gmr,&delta_r);
        
        let nnt = SymmetricTensor3D::from(
            [cp*cp*st*st, cp*sp*st*st, cp*st*ct,
                          sp*sp*st*st, sp*st*ct,
                                          ct*ct]);
        let ref_ten = freq*&(&SymmetricTensor3D::eye() - &(3.0*&nnt));

        let err_ten = &ten - &ref_ten; 

        let tol = freq*1e-9;
        assert!(err_ten.xx().abs() < tol);
        assert!(err_ten.xy().abs() < tol);
        assert!(err_ten.xz().abs() < tol);
        assert!(err_ten.yy().abs() < tol);
        assert!(err_ten.yz().abs() < tol);
        assert!(err_ten.zz().abs() < tol);


      }
    }
  }

  //----------------------------------------------------------------------------
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
  //----------------------------------------------------------------------------
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
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
