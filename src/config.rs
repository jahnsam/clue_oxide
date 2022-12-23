use strum::IntoEnumIterator;

use crate::clue_errors::*;
use crate::lexer::*;
use super::particle_config::*;
use super::particle_specifier::SpecifiedParticle;
use super::physical_constants::*;
use super::vector3::*;


/// Config contains all the setting for CluE.
#[derive(Debug,Clone)]
pub struct Config{
  pub radius: Option<f64>,
  pub inner_radius: Option<f64>,
  //max_number_of_cells: usize,
  //pbc_style: PBCStyle,
  //error_tolerance: f64,
  pub magnetic_field: Option<Vector3>,
  pub central_spin_coordinates: Option<CentralSpinCoordinates>,
  //use_periodic_boundary_conditions: bool,
  pub particles: Vec::<ParticleConfig>,
}

impl Default for Config{
  fn default() -> Self {
    let mut particles = Vec::<ParticleConfig>::new();
    for element in Element::iter(){

      let mut particle_config = ParticleConfig::new();
      particle_config.filter.elements.push(element);

      particle_config.config.isotopic_distribution.isotope_abundances.push(
          IsotopeAbundance{
            isotope: Isotope::most_common_for(&element),
            abundance: 1.0,
          });


      particles.push(particle_config);
    }

    Config{
      radius: None,
      inner_radius: None,
      magnetic_field: None,
      central_spin_coordinates: None,
      particles,
    }
  }
}

impl Config{
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------

  pub fn get_max_spin_multiplicity_for_any_isotope(&self, element: Element)
    -> usize{

      let mut target = SpecifiedParticle::new();
      target.element = Some(element);

      let found_configs = find_particle_configs(&target, &self.particles);

      let mut max_spin_multiplicity = 1;

      for particle_config in found_configs{
        for iso_abu in &particle_config
          .config.isotopic_distribution.isotope_abundances{

          max_spin_multiplicity = max_spin_multiplicity.max(
              iso_abu.isotope.spin_multiplicity() );
        }
      }

      max_spin_multiplicity
    }
}
#[derive(Debug,Clone)]
pub enum CentralSpinCoordinates{
  Atoms (Vec::<u32>),
  XYZ (Vector3),
}


#[derive(Debug,Clone)]
pub enum PBCSyle{
  CRYST1,
  TIGHT,
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Config{
  fn read_file(filename: &str) -> Result<Self,CluEError>{
  
    let mut config = Config::new();

    let token_stream = get_tokens_from_file(filename)?;

    let mut mode = ModeAttribute::new();


    for expression in token_stream.iter(){
      if expression.lhs.is_empty(){ continue; }

      if let Token::Mode(new_mode) = &expression.lhs[0]{
        mode = new_mode.clone();
        continue;
      }

      match mode.mode{
        ConfigMode::Clusters => (),
        ConfigMode::Config =>  config.parse_config_line(expression)?,
        ConfigMode::Filter => (),
        ConfigMode::Spins => (),
        ConfigMode::Structures => (),
        ConfigMode::Tensors => (),
      }

    }


    Ok(config)
  }
  //----------------------------------------------------------------------------
  fn parse_config_line(&mut self, expression: &TokenExpression) 
    -> Result<(),CluEError>
  {
  
    let line_number = expression.line_number;
    let already_set = ||{
      CluEError::OptionAlreadySet(line_number,expression.lhs[0].to_string()) };

    //let some_f64 = ||{line_number,expression.rhs[0].to_f64_token()}

    match expression.lhs[0]{
      Token::Radius => {
        if let Some(r) = self.radius{
          return Err(already_set())
        }
        //self.radius = some_f64()?;
      },
      _ => return Err(CluEError::InvalidToken(line_number,
            expression.lhs[0].to_string())),
    }
    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

