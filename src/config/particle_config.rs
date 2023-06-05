
use crate::structure::particle_filter::{ParticleFilter,VectorSpecifier,
  SecondaryParticleFilter};
//use super::particle_specifier::*;
use crate::physical_constants::*;
use std::collections::HashMap;

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
#[derive(Debug,Clone,PartialEq)]
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
  //----------------------------------------------------------------------------
  pub fn max_possible_spin_multiplicity(&self) -> Option<usize> {
    
    let Some(properties) = &self.properties else{
      return None;
    };

    if properties.isotopic_distribution.isotope_abundances.is_empty(){
      return None;
    }

    let mut max_spin_mult = 0;
    for isotope_abundance in properties.isotopic_distribution
      .isotope_abundances.iter()
    {
      max_spin_mult = usize::max(max_spin_mult,
          isotope_abundance.isotope.spin_multiplicity());
    }

    Some(max_spin_mult)
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq,Default)]
pub struct ParticleProperties{
  pub cosubstitute: Option<SecondaryParticleFilter>,
  pub extracell_isotopic_distribution:  Option<IsotopeDistribution>,
  pub isotopic_distribution:  IsotopeDistribution,
  pub isotope_properties: HashMap::<String,IsotopeProperties>,
}

impl ParticleProperties{
  pub fn new() -> Self{
    Default::default()
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq,Default)]
pub struct IsotopeProperties{
  pub electric_quadrupole_coupling: Option<TensorSpecifier>,
  pub exchange_coupling: Option<f64>,
  pub hyperfine_coupling: Option<TensorSpecifier>,
}

impl IsotopeProperties{
  pub fn new() -> Self{ IsotopeProperties::default() }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq,Default)]
pub struct IsotopeDistribution{
  pub isotope_abundances: Vec::<IsotopeAbundance>,
  pub extracell_void_probability: Option<f64>,
}

#[derive(Debug,Clone,PartialEq)]
pub struct IsotopeAbundance{
  pub isotope: Isotope,
  pub abundance: f64,         
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq,Default)]
pub struct TensorSpecifier{
  pub values: Option<[f64; 3]>,
  pub x_axis: Option<VectorSpecifier>,
  pub y_axis: Option<VectorSpecifier>,
  pub z_axis: Option<VectorSpecifier>,
}

impl TensorSpecifier{
  pub fn new() -> Self{
    TensorSpecifier::default()
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
  /*
#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb::read_pdb;

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
}
  */
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

