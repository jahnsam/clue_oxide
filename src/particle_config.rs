use super::particle_specifier::*;
use super::physical_constants::*;
use super::vector3::Vector3;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn find_particle_configs<'a>(
    target: &SpecifiedParticle,
    particle_configs: &'a Vec::<ParticleConfig>) 
  -> Vec::<&'a ParticleConfig>{

    let mut found_configs = Vec::<&ParticleConfig>::with_capacity(
        particle_configs.len());

    for ii in 0..particle_configs.len(){
      if particle_configs[ii].filter.covers(target){
        found_configs.push(&particle_configs[ii]);
      }
    }

    found_configs
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct ParticleConfig{
  pub filter: ParticleSpecifier,
  pub config: ParticleProperties,
}
impl Default for ParticleConfig{
  fn default() -> Self{
    ParticleConfig{
      filter: ParticleSpecifier::new(),
      config: ParticleProperties::default(),
    }
  }
}
impl ParticleConfig{
  pub fn new() -> Self{
    Default::default()
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
  pub x_axis: Option<AxisSpecifier>,
  pub y_axis: Option<AxisSpecifier>,
  pub z_axis: Option<AxisSpecifier>,
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
  Other(ParticleSpecifier),
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use super::super::pdb::read_pdb;

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
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

