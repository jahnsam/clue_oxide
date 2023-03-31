

//use strum::IntoEnumIterator;

use crate::clue_errors::*;
use crate::config::lexer::*;
use crate::config::token::*;
use crate::config::command_line_input::CommandLineInput;
use crate::config::token_algebra::*;
//use crate::config::token_stream;
use crate::config::token_expressions::*;
use crate::config::particle_config::ParticleConfig;//, ParticleProperties,  IsotopeAbundance};
use crate::physical_constants::Isotope;
use crate::space_3d::Vector3D;
use crate::integration_grid::IntegrationGrid;


pub mod lexer;
pub mod token;
pub mod token_algebra;
pub mod token_stream;
pub mod token_expressions;
pub mod particle_config;
pub mod command_line_input;
pub mod parse_filter;
pub mod parse_properties;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// Config contains all the setting for CluE.
#[derive(Debug,Clone,Default)]
pub struct Config{
  pub cluster_batch_size: Option<usize>,
  pub cluster_method: Option<ClusterMethod>,
  pub detected_spin_position: Option<DetectedSpinCoordinates>,
  pub detected_spin_transition: Option<[usize;2]>,
  pub detected_spin_identity: Option<Isotope>,
  pub input_structure_file: Option<String>,
  pub load_geometry: Option<LoadGeometry>,
  pub magnetic_field: Option<Vector3D>,
  pub max_cluster_size: Option<usize>,
  pub neighbor_cutoffs: Vec::<NeighborCutoff>,
  pub number_timepoints: Vec::<usize>,
  //pub inner_radius: Option<f64>,
  //pbc_style: PBCStyle,
  //error_tolerance: f64,
  //use_periodic_boundary_conditions: bool,
  pub particles: Vec::<ParticleConfig>, // TODO: why in Option?
  pub pdb_model_index: Option<usize>,
  pub pulse_sequence: Option<PulseSequence>,
  pub radius: Option<f64>,
  pub rng_seed: Option<u64>,
  time_axis: Vec::<f64>,
  pub time_increments: Vec::<f64>,
  pub write_structure_pdb: Option<String>,
}


impl Config{
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------
  /// This function sets config setting that are necessary to be `Some`,
  /// but will likely be the same for most simulations.  
  /// Any field that is already `Some` will be left alone.
  pub fn set_defaults(&mut self){
    if self.cluster_batch_size == None{
      self.cluster_batch_size = Some(10000);
    }
    if self.detected_spin_identity == None{
      self.detected_spin_identity = Some(Isotope::Electron);
    }
    if self.detected_spin_transition == None{
      self.detected_spin_transition = Some([0,1]);
    }
    if self.load_geometry == None{
      self.load_geometry = Some(LoadGeometry::Sphere);
    }
    if self.pdb_model_index == None{
      self.pdb_model_index = Some(0);
    }
  }
  //----------------------------------------------------------------------------
  pub fn max_spin_multiplicity_for_particle_config(&self, id: usize) -> usize{

   if id >= self.particles.len() {return 0;}

   self.particles[id].max_possible_spin_multiplicity()
  }
  //----------------------------------------------------------------------------
  pub fn get_time_axis(&self) -> Result<&Vec::<f64>,CluEError>
  {
    if self.time_axis.is_empty() {
      return Err(CluEError::NoTimeAxis);
    };
    Ok(&self.time_axis)
  }
  //----------------------------------------------------------------------------
  pub fn construct_time_axis(&mut self) -> Result<(),CluEError>
  {
    let dts = &self.time_increments;
    if dts.is_empty(){
      return Err(CluEError::NoTimeIncrements);
    }

    let n_dts = &self.number_timepoints;
    if n_dts.is_empty(){
      return Err(CluEError::NoTimepoints);
    }

    if dts.len() != n_dts.len(){
      return Err(CluEError::LenghMismatchTimepointsIncrements(
            n_dts.len(),dts.len()));
    }

    let Some(pulse_sequence) = &self.pulse_sequence else{
      return Err(CluEError::NoPulseSequence);
    };

    let n_tot = n_dts.iter().sum::<usize>();
    let mut time_axis = Vec::<f64>::with_capacity(n_tot);
    match pulse_sequence{
      PulseSequence::CarrPurcell(n_pi_pulses) => {
        for (idx, &n_dt) in n_dts.iter().enumerate(){
          let dt = (*n_pi_pulses as f64)*dts[idx];
          let mut t = 0.0;
          for _ii in 0..n_dt{
            time_axis.push(t);
            t += dt;
          }
        }
      }
    }

    Ok(())
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[derive(Debug,Clone,PartialEq)]
pub enum LoadGeometry{
  Cube,
  Sphere,
}
#[derive(Debug,Clone,PartialEq)]
pub enum DetectedSpinCoordinates{
  CentroidOverSerials(Vec::<u32>),
  XYZ(Vector3D),
  ProbabilityDistribution(IntegrationGrid),
}
#[derive(Debug,Clone,PartialEq)]
pub enum NeighborCutoff{
  DeltaHyperfine(f64),
  DipoleDipole(f64),
  HahnThreeSpinModulationDepth(f64),
  HahnThreeSpinFourthOrderTaylorCoefficient(f64),
}
#[derive(Debug,Clone,PartialEq)]
pub enum ClusterMethod{
  AnalyticRestricted2CCE,
  CCE,
}
#[derive(Debug,Clone,PartialEq)]
pub enum PulseSequence{
  CarrPurcell(usize),
  //RefocusedEcho,
}
/*


#[derive(Debug,Clone)]
pub enum PBCSyle{
  CRYST1,
  TIGHT,
}
*/
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Config{
  pub fn read_input(input: CommandLineInput) -> Result<Self,CluEError>{
    let Some(filename) = &input.config_file else{
      return Err(CluEError::NoInputFile);
    };
  
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
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClusterMethod => { 
        if let Some(_value) = &self.cluster_method{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }
        match rhs[0]{
          Token::CCE => self.cluster_method = Some(ClusterMethod::CCE),
          Token::R2CCE => self.cluster_method 
            = Some(ClusterMethod::AnalyticRestricted2CCE),
          _  => return Err(CluEError::InvalidArgument(expression.line_number,
                "valid cluster method are \"CCE\" and \"r2CCE\"".to_string())),
        }

      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::DetectedSpinPosition =>{
        if let Some(_value) = &self.detected_spin_position{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }

        match rhs[0]{
          Token::CentroidOverSerials =>{
            let args = get_function_arguments(rhs,expression.line_number)?;

            let mut serials = Vec::<u32>::new();
            let value_token = to_i32_token(args, expression.line_number)?;
            match value_token {
              Token::Int(a) => serials.push(a as u32),
              Token::VectorI32(v) => serials = v.iter().map(|&a| a as u32)
                .collect(),  
              _  => return Err(CluEError::InvalidArgument(expression.line_number,
                    "list of positive integers".to_string())),
            }
            self.detected_spin_position = Some(
                DetectedSpinCoordinates::CentroidOverSerials(serials));

          },
          _ => return Err(CluEError::InvalidToken(expression.line_number,
                rhs[0].to_string())),
        }
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::MaxClusterSize 
        => set_to_some_usize(&mut self.max_cluster_size, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NumberTimepoints 
        => set_to_vec_usize(&mut self.number_timepoints, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::InputStructureFile 
        => set_to_some_string(&mut self.input_structure_file, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      Token::PulseSequence => {
        if let Some(_value) = &self.pulse_sequence{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }
        match rhs[0]{
          Token::Hahn 
            => self.pulse_sequence = Some(PulseSequence::CarrPurcell(1)),
          Token::CarrPurcell => {
            if rhs.len() != 3 || rhs[1] != Token::Minus{
              return Err(CluEError::InvalidPulseSequence(
                    expression.line_number));
            }
            let mut n_pi_opt: Option<usize> = None;
            let rhs = token_stream::read_strings_as_integers(rhs.clone(), 
                expression.line_number)?;
            if let Token::Int(n_pi) = rhs[2]{
              self.pulse_sequence 
                = Some(PulseSequence::CarrPurcell(n_pi as usize));
            }else{
              return Err(CluEError::InvalidPulseSequence(
                    expression.line_number));
            }

          },
          _ => return Err(CluEError::InvalidPulseSequence(
                expression.line_number)),
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Radius => set_to_some_f64(&mut self.radius,expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::RNGSeed => {
        if let Some(_value) = &self.rng_seed{
          return Err(already_set());
        }
        let mut rng_seed: Option<i32> = None;
        set_to_some_i32(&mut rng_seed,expression)?;

        if let Some(seed) = rng_seed{
          self.rng_seed = Some( seed as u64);
        }else{
          return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string()))
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::TimeIncrements 
        => set_to_vec_f64(&mut self.time_increments, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteStructurePDB 
        => set_to_some_string(&mut self.write_structure_pdb, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
    }
    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  #[test]
  fn test_parse_config_line(){
    let expressions = get_tokens_from_line("\
        cluster_method = CCE;
        detected_spin_position = centroid_over_serials([28,29]);
        input_structure_file = \"../../assets/TEMPO_wat_gly_70A.pdb\";
        max_cluster_size = 4;
        number_timepoints = [40,60];
        pulse_sequence = CP-1;
        radius = 80e-10;
        time_increments = [1e-9, 5e-7];
        write_structure_pdb = out.pdb;
        ").unwrap();

    let mut config = Config::new();
    for expression in expressions.iter(){
      config.parse_config_line(expression).unwrap();
    }


    assert_eq!(config.cluster_method,Some(ClusterMethod::CCE));
    assert_eq!(config.detected_spin_position, 
        Some(DetectedSpinCoordinates::CentroidOverSerials(vec![28,29])));
    assert_eq!(config.input_structure_file, 
         Some("../../assets/TEMPO_wat_gly_70A.pdb".to_string()));
    
    assert_eq!(config.max_cluster_size, Some(4));
    assert_eq!(config.number_timepoints, vec![40,60]);
    assert_eq!(config.pulse_sequence, Some(PulseSequence::CarrPurcell(1)));
    assert_eq!(config.radius, Some(80.0e-10));
    assert_eq!(config.time_increments, vec![1e-9,5e-7]);
    assert_eq!(config.write_structure_pdb, Some("out.pdb".to_string()));
    
  }
}

