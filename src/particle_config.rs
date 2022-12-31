
use crate::particle_filter::{ParticleFilter,VectorSpecifier};
//use super::particle_specifier::*;
use super::physical_constants::*;
use super::vector3::Vector3;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
pub fn find_particle_configs<'a>(
    target: &SpecifiedParticle,
    particle_configs: &'a [ParticleConfig]) 
  -> Vec::<&'a ParticleConfig>{

    let mut found_configs = Vec::<&ParticleConfig>::with_capacity(
        particle_configs.len());

    for particle_config in particle_configs.iter(){
      if particle_config.filter.covers(target){
        found_configs.push(particle_config);
      }
    }

    found_configs
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct ParticleConfig{
  pub label: String,
  pub filter: Option<ParticleFilter>,
  pub properties: Option<ParticleProperties>,
}
impl ParticleConfig{
  pub fn new(label: String) -> Self{
    ParticleConfig{
      label,
      filter: None,
      properties: None,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct ParticleProperties{
  pub isotopic_distribution:  IsotopeDistribution,
  pub exchange_coupling: f64,
  pub hyperfine: Option<TensorSpecifier>,
  pub electric_quadrupole_coupling: Option<TensorSpecifier>,
}
impl Default for ParticleProperties{
  fn default() -> Self{
    ParticleProperties{
      isotopic_distribution:  IsotopeDistribution::default(),
      exchange_coupling: 0.0,
      hyperfine: None,
      electric_quadrupole_coupling: None,

    }
  }
}
impl ParticleProperties{
  pub fn new() -> Self{
    Default::default()
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct IsotopeDistribution{
  pub isotope_abundances: Vec::<IsotopeAbundance>,
  pub force_no_pbc: bool,
  pub extracell_void_probability: f64,
}

impl Default for IsotopeDistribution{
  fn default() -> Self{
    IsotopeDistribution{
      isotope_abundances: Vec::<IsotopeAbundance>::new(),
      force_no_pbc: false,
      extracell_void_probability: 0.0,
    }
  }
}

#[derive(Debug,Clone)]
pub struct IsotopeAbundance{
  pub isotope: Isotope,
  pub abundance: f64,         
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct TensorSpecifier{
  pub values: [f64; 3],
  pub x_axis: Option<VectorSpecifier>,
  pub y_axis: Option<VectorSpecifier>,
  pub z_axis: Option<VectorSpecifier>,
}

impl Default for TensorSpecifier{
  fn default() -> Self{
  
    TensorSpecifier{
      values: [0.0; 3],
      x_axis: None,
      y_axis: None,
      z_axis: None,
   }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
// TODO: deprected
#[derive(Debug,Clone)]
pub enum AxisSpecifier{
  Connected(ParticleSelector, ParticleSelector),
  PDBIndices(Vec::<u32>),
  Random,
  SameResidueSequence( ParticleSelector,ParticleSelector ),
  Vector(Vector3),
}

#[derive(Debug,Clone)]
pub enum ParticleSelector{
  This,
  Other(Box<ParticleSpecifier>),
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb::read_pdb;

  /*
  #[test]
  fn test_find_particle_configs(){
    let filename = "./assets/TEMPO.pdb";
    let pdb = read_pdb(filename).expect("Could not read pdb file.");  

    assert_eq!(Element::Nitrogen,pdb.element(27));
    let target = SpecifiedParticle::specify(27, &pdb);


    let number = 3;
    let mut configs = Vec::<ParticleConfig>::with_capacity(number);
    for _ii in 0..number{
      configs.push(ParticleConfig::new());
    }

    configs[1].filter.elements.push(Element::Hydrogen);
    configs[2].filter.elements.push(Element::Nitrogen);

    let found_configs = find_particle_configs(&target, &configs);
    assert_eq!(2,found_configs.len());

    for ii in 0..found_configs.len() {
      assert!(found_configs[ii].filter.covers(&target));
    }
  }
  */
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

