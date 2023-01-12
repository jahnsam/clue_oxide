

//use strum::IntoEnumIterator;

use crate::clue_errors::*;
use crate::config::lexer::*;
use crate::config::token::*;
use crate::command_line_input::CommandLineInput;
//use crate::structure::particle_filter::ParticleFilter;
use crate::config::token_algebra::*;
//use crate::config::token_stream;
use crate::config::particle_config::ParticleConfig;//, ParticleProperties,  IsotopeAbundance};
//use crate::physical_constants::*;
use crate::space_3d::Vector3D;

pub mod lexer;
pub mod token;
pub mod token_algebra;
pub mod token_stream;
pub mod particle_config;

/// Config contains all the setting for CluE.
#[derive(Debug,Clone,Default)]
pub struct Config{
  pub radius: Option<f64>,
  //pub inner_radius: Option<f64>,
  //max_number_of_cells: usize,
  //pbc_style: PBCStyle,
  //error_tolerance: f64,
  pub magnetic_field: Option<Vector3D>,
  //pub central_spin_coordinates: Option<CentralSpinCoordinates>,
  //use_periodic_boundary_conditions: bool,
  pub particles: Option<Vec::<ParticleConfig>>,
}
/*
impl Default for Config{
  fn default() -> Self {
    /*
    let mut particles = Vec::<ParticleConfig>::new();
    for element in Element::iter(){

      let mut filter = ParticleFilter::new();
      filter.elements.push(element);

      let name: String = format!("default_{}",element.to_string());
      let mut particle_config = ParticleConfig::new(name);
      particle_config.filter = Some(filter);

      let mut properties = ParticleProperties::new();
      properties.isotopic_distribution.isotope_abundances.push(
          IsotopeAbundance{
            isotope: Isotope::most_common_for(&element),
            abundance: 1.0,
          });

      particle_config.properties = Some(properties);

      particles.push(particle_config);
    }
    */

    Config{
      radius: None,
      //inner_radius: None,
      magnetic_field: None,
      central_spin_coordinates: None,
      particles: None,
    }
  }
}
*/
impl Config{
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------

  /*
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
    */
}
/*
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
*/
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Config{
  pub fn read_input(input: CommandLineInput) -> Result<Self,CluEError>{
    let filename = &input.config_file.unwrap();
  
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

    if expression.relationship != Some(Token::Equals){
      return Err(CluEError::ExpectedEquality(expression.line_number));
    }

  let already_set = ||{
    CluEError::OptionAlreadySet(
        expression.line_number, expression.lhs[0].to_string()) };

    match expression.lhs[0]{
      Token::MagneticField => {
  
        if let Some(_value) = &self.magnetic_field{
          return Err(already_set());
        }
        let mut mf_opt: Option<f64> = None;
        set_to_some_f64(&mut mf_opt, expression)?;

        if let Some(bz) = mf_opt{
          self.magnetic_field = Some( Vector3D::from([0.0,0.0,bz]) );
        }else{
          return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string()))
        }

      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Radius => set_to_some_f64(&mut self.radius,expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
    }
    Ok(())
  }
}
//------------------------------------------------------------------------------
fn check_target<T>(target: &Option<T>, expression: &TokenExpression)
  -> Result<(),CluEError>
{ 
  let already_set = ||{
    CluEError::OptionAlreadySet(
        expression.line_number, expression.lhs[0].to_string()) };

  if let Some(_value) = target{
    return Err(already_set());
  }

  Ok(())
}

//------------------------------------------------------------------------------
fn set_to_some_f64(target: &mut Option<f64>, expression: &TokenExpression) 
  -> Result<(),CluEError>
{ 

  check_target(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_f64_token(tokens, expression.line_number)?;
  if let Token::VectorF64(vec) = value_token {
    if vec.len() != 1{
      return Err(CluEError::ExpectedFloatRHS(expression.line_number));
    } 
    *target = Some(vec[0]);
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//------------------------------------------------------------------------------
fn set_to_some_vector3d(
    target: &mut Option<Vector3D>, expression: &TokenExpression) 
  -> Result<(),CluEError>
{ 

  check_target(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_f64_token(tokens, expression.line_number)?;
  if let Token::VectorF64(vec) = value_token {
    if vec.len() != 3{
      return Err(CluEError::ExpectedVecOfNFloatsRHS(expression.line_number,3));
    } 
    *target = Some(Vector3D::from([vec[0],vec[1],vec[2]]));
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

