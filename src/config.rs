use crate::clue_errors::*;

use crate::config::command_line_input::CommandLineInput;
use crate::config::lexer::*;
use crate::config::token::*;
use crate::config::token_algebra::*;
use crate::config::token_expressions::*;
use crate::config::particle_config::{ParticleConfig,TensorSpecifier};

use crate::physical_constants::{ANGSTROM,Isotope};
use crate::space_3d::Vector3D;
use crate::integration_grid::IntegrationGrid;



use crate::io;


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
  pub clash_distance_pbc: Option<f64>,
  pub cluster_batch_size: Option<usize>,
  pub cluster_method: Option<ClusterMethod>,
  pub density_matrix: Option<DensityMatrixMethod>,
  pub detected_spin_g_matrix: Option<TensorSpecifier>,
  pub detected_spin_identity: Option<Isotope>,
  pub detected_spin_multiplicity: Option<usize>,
  pub detected_spin_position: Option<DetectedSpinCoordinates>,
  pub detected_spin_transition: Option<[usize;2]>,
  pub input_structure_file: Option<String>,
  pub load_geometry: Option<LoadGeometry>,
  pub magnetic_field: Option<Vector3D>,
  pub max_cluster_size: Option<usize>,
  pub neighbor_cutoff_delta_hyperfine: Option<f64>,
  pub neighbor_cutoff_dipole_dipole: Option<f64>,
  pub neighbor_cutoff_distance: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_mod_depth: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_taylor_4: Option<f64>,
  pub number_timepoints: Vec::<usize>,
  //pub inner_radius: Option<f64>,
  //pbc_style: PBCStyle,
  //error_tolerance: f64,
  //use_periodic_boundary_conditions: bool,
  pub orientation_grid: Option<OrientationAveraging>,
  pub particles: Vec::<ParticleConfig>,
  pub pdb_model_index: Option<usize>,
  pub pulse_sequence: Option<PulseSequence>,
  pub remove_partial_methyls: Option<bool>,
  pub root_dir: Option<String>,
  pub radius: Option<f64>,
  pub rng_seed: Option<u64>,
  pub save_name: Option<String>,
  pub temperature: Option<f64>,
  time_axis: Vec::<f64>,
  pub time_increments: Vec::<f64>,
  pub write_auxiliary_signals: Option<String>,
  pub write_bath: Option<String>,
  pub write_clusters: Option<String>,
  pub write_info: Option<String>,
  pub write_exchange_groups: Option<String>,
  pub write_orientation_signals: Option<String>,
  pub write_structure_pdb: Option<String>,
  pub write_tensors: Option<String>,
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
    if self.cluster_batch_size.is_none(){
      self.cluster_batch_size = Some(10000);
    }
    if self.detected_spin_identity.is_none(){
      self.detected_spin_identity = Some(Isotope::Electron);
    }
    if self.detected_spin_transition.is_none(){
      self.detected_spin_transition = Some([0,1]);
    }
    if self.load_geometry.is_none(){
      self.load_geometry = Some(LoadGeometry::Sphere);
    }
    if self.pdb_model_index.is_none(){
      self.pdb_model_index = Some(0);
    }
    if self.density_matrix.is_none(){
      match self.temperature{
        Some(_) => self.density_matrix = Some(DensityMatrixMethod::Thermal),
        None => self.density_matrix = Some(DensityMatrixMethod::Identity),
      }
    }

    if self.remove_partial_methyls.is_none(){
      'part_met : for particle_config in self.particles.iter(){
        if let Some(properties) = &particle_config.properties{
          for (_key, isotope_props) in properties.isotope_properties.iter(){
            if isotope_props.exchange_coupling.is_some(){
              self.remove_partial_methyls = Some(true);
              break 'part_met;
            }
          }
        }
      }
    }

    if self.root_dir.is_none(){
      self.root_dir = Some("./".to_string());
    }
    if self.save_name.is_none(){
      self.save_name = Some("CluE-".to_string());
    }


    let set_default_write_path = |write_opt: &mut Option<String>, default: &str|
    {
      match write_opt{
        None => *write_opt = Some(default.to_string()),
        Some(path) => if path.is_empty(){ *write_opt = None;},
      }

    };

    set_default_write_path(&mut self.write_auxiliary_signals,
        "auxiliary_signals");

    set_default_write_path(&mut self.write_bath,
        "bath");

    set_default_write_path(&mut self.write_clusters,
        "clusters");

    set_default_write_path(&mut self.write_info,
        "info");

    set_default_write_path(&mut self.write_exchange_groups,
        "exchange_groups");

    set_default_write_path(&mut self.write_orientation_signals,
        "orientations");

    set_default_write_path(&mut self.write_structure_pdb,
        "spin_system");

    // The tensors file is large, so do not write it by default.
    /*
    set_default_write_path(&mut self.write_tensors,
        "tensors");
    */
  }
  //----------------------------------------------------------------------------
  pub fn max_spin_multiplicity_for_particle_config(&self, id: usize) 
    -> Option<usize>
  {

   if id >= self.particles.len() {return None;}

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
    self.time_axis = Vec::<f64>::with_capacity(n_tot);
    match pulse_sequence{
      PulseSequence::CarrPurcell(n_pi_pulses) => {
        let mut t = 0.0;
        for (idx, &n_dt) in n_dts.iter().enumerate(){
          let dt = (*n_pi_pulses as f64 + 1.0)*dts[idx];
          for _ii in 0..n_dt{
            self.time_axis.push(t);
            t += dt;
          }
        }
      }
    }

    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn write_time_axis(&self,save_path: String) -> Result<(),CluEError>{
    io::write_data(& vec![self.time_axis.clone()], 
        &format!("{}/time_axis.csv",save_path), vec!["time_axis".to_string()])
  }
  //----------------------------------------------------------------------------
  pub fn find_particle_config(&self, label: &str) 
    -> Option<(usize,&ParticleConfig)>
  {
    for (idx, p_cfg) in self.particles.iter().enumerate(){
      if *label == p_cfg.label{
        return Some((idx,p_cfg));
      }
    }
    None
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[derive(Debug,Clone,PartialEq)]
pub enum DensityMatrixMethod{
  Identity,
  Thermal,
}

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
/*
#[derive(Debug,Clone,PartialEq)]
pub enum NeighborCutoff{
  DeltaHyperfine(f64),
  DipoleDipole(f64),
  HahnThreeSpinModulationDepth(f64),
  HahnThreeSpinFourthOrderTaylorCoefficient(f64),
}
*/
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

#[derive(Debug,Clone,PartialEq)]
pub enum OrientationAveraging{
  Grid(IntegrationGrid),
  Lebedev(usize),
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

    config.parse_token_stream(token_stream)?;

    Ok(config)
  }
  //----------------------------------------------------------------------------
  /// This function updates config from a Vec<TokenExpression>.
  pub fn parse_token_stream(&mut self,token_stream: Vec<TokenExpression>) 
    -> Result<(),CluEError>
  {
    let mut mode = ModeAttribute::new();

    for expression in token_stream.iter(){
      if expression.lhs.is_empty(){ continue; }

      if let Token::Mode(new_mode) = &expression.lhs[0]{
        mode = new_mode.clone();
        continue;
      }

      match mode.mode{
        ConfigMode::Clusters => (),
        ConfigMode::Config =>  self.parse_config_line(expression)?,
        ConfigMode::Filter => self.parse_filter_line(expression,&mode.label)?,
        ConfigMode::SpinProperties => {
          let Some(label) = &mode.label else{
            return Err(CluEError::SpinPropertiesNeedsALabel);
          };
          if mode.isotope.is_none(){
            return Err(CluEError::SpinPropertiesNeedsAnIsotope(label.clone()));
          }
          self.parse_properties_line(expression,&mode.label,mode.isotope)?;
        },
        ConfigMode::Spins => (),
        ConfigMode::StructureProperties => {
          let Some(label) = &mode.label else{
            return Err(CluEError::StructurePropertiesNeedsALabel);
          };
          if mode.isotope.is_some(){
            return Err(
               CluEError::StructurePropertiesDoesNotNeedAnIsotope(
                 label.clone()));
          }
          self.parse_properties_line(expression,&mode.label,mode.isotope)?;
        },
        ConfigMode::Structures => (),
        ConfigMode::Tensors => (),
      }

    }


    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function parses a line of tokens and modifies the `Config`. 
  pub fn parse_config_line(&mut self, expression: &TokenExpression) 
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
      Token::ClashDistancePBC => {
        set_to_some_f64(&mut self.clash_distance_pbc, expression)?;
        if let Some(r) = &mut self.clash_distance_pbc{
            *r *= ANGSTROM;
        }else{ return Err(CluEError::NoClashDistancePBC);}
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClusterBatchSize 
        => set_to_some_usize(&mut self.cluster_batch_size, expression)?, 
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
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
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
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ReadCSV => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let grid = IntegrationGrid::read_from_csv(&(args[0].to_string()))?;

            if grid.dim() != 3{
              return Err(CluEError::WrongProbabilityDistributionDim(
                    expression.line_number,3,grid.dim()));
            }
            self.detected_spin_position = Some(
                DetectedSpinCoordinates::ProbabilityDistribution(grid));

          }
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::SquareBracketOpen =>{
            let mut r_opt: Option<Vector3D> = None;
            set_to_some_vector3d(&mut r_opt, expression)?;
            if let Some(r) = r_opt{
              self.detected_spin_position 
                = Some(DetectedSpinCoordinates::XYZ(r));
            }else{
              return Err(CluEError::NoCentralSpinCoor);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
                rhs[0].to_string())),
        }
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::InputStructureFile 
        => set_to_some_string(&mut self.input_structure_file, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::LoadGeometry =>{
        if let Some(_value) = &self.load_geometry{
          return Err(already_set());
        }

        let Some(rhs) = &expression.rhs else { 
          return Err(CluEError::NoRHS(expression.line_number));
        };

        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }else if rhs.len() > 1{
          return Err(CluEError::TooManyRHSArguments(expression.line_number));
        }

        match rhs[0]{
          Token::Cube => self.load_geometry = Some(LoadGeometry::Cube),
          Token::Sphere => self.load_geometry = Some(LoadGeometry::Sphere),
          _ => return Err(CluEError::InvalidGeometry(expression.line_number,
                rhs[0].to_string())),
        }

      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::MaxClusterSize 
        => set_to_some_usize(&mut self.max_cluster_size, expression)?,
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
      Token::NeighborCutoffDeltaHyperfine 
        => set_to_some_f64(&mut self.neighbor_cutoff_delta_hyperfine,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffDipoleDipole 
        => set_to_some_f64(&mut self.neighbor_cutoff_dipole_dipole,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffDistance 
        => {
          set_to_some_f64(&mut self.neighbor_cutoff_distance,expression)?;
          if let Some(r) = &mut self.neighbor_cutoff_distance{
            *r *= ANGSTROM;
          }else{
            return Err(CluEError::NoNeighborCutoffDistance);
          }
        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoff3SpinHahnModDepth 
        => set_to_some_f64(&mut self.neighbor_cutoff_3_spin_hahn_mod_depth,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoff3SpinHahnTaylor4 
        => set_to_some_f64(&mut self.neighbor_cutoff_3_spin_hahn_taylor_4,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NumberTimepoints 
        => set_to_vec_usize(&mut self.number_timepoints, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::OrientationGrid => {
        if let Some(_value) = &self.orientation_grid{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }

        match rhs[0]{
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::Lebedev => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let mut n_grid_opt: Option<usize> = None;
            set_to_some_usize(&mut n_grid_opt,expression)?;
            if let Some(n_grid) = n_grid_opt{
              self.orientation_grid 
                = Some(OrientationAveraging::Lebedev(n_grid));
            }else{
              return Err(CluEError::NoOrientationGrid);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ReadCSV => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let grid = IntegrationGrid::read_from_csv(&(args[0].to_string()))?;

            if grid.dim() != 3{
              return Err(CluEError::WrongOrientationGridDim(
                    expression.line_number,3,grid.dim()));
            }
            self.orientation_grid = Some(
                OrientationAveraging::Grid(grid));

          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::SquareBracketOpen =>{
            let mut r_opt: Option<Vector3D> = None;
            set_to_some_vector3d(&mut r_opt, expression)?;
            if let Some(r) = r_opt{
              let mut grid = IntegrationGrid::new(3);
              grid.push(vec![r.x(),r.y(),r.z()],1.0)?;
              self.orientation_grid
                = Some(OrientationAveraging::Grid(grid));
            }else{
              return Err(CluEError::NoOrientationGrid);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
                rhs[0].to_string())),
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
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
      Token::RemovePartialMethyls 
        => set_to_some_bool(&mut self.remove_partial_methyls,expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Radius => {
        set_to_some_f64(&mut self.radius,expression)?;
        if let Some(r) = &mut self.radius{
            *r *= ANGSTROM;
        }else{ return Err(CluEError::NoRadius);}
      },
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
      Token::Temperature
        => set_to_some_f64(&mut self.temperature, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::TimeIncrements 
        => set_to_vec_f64(&mut self.time_increments, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteAuxiliarySignals
        => set_to_some_string(&mut self.write_auxiliary_signals, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteBath
        => set_to_some_string(&mut self.write_bath, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteClusters
        => set_to_some_string(&mut self.write_clusters, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteExchangeGroups
        => set_to_some_string(&mut self.write_exchange_groups, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteInfo
        => set_to_some_string(&mut self.write_info, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteOrientationSignals
        => set_to_some_string(&mut self.write_orientation_signals, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteStructurePDB 
        => set_to_some_string(&mut self.write_structure_pdb, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteTensors 
        => set_to_some_string(&mut self.write_tensors, expression)?,
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
        clash_distance_pbc = 0.1;
        cluster_batch_size = 20000;
        cluster_method = cce;
        detected_spin_position = centroid_over_serials([28,29]);
        input_structure_file = \"../../assets/TEMPO_wat_gly_70A.pdb\";
        load_geometry = cube;
        max_cluster_size = 4;
        number_timepoints = [40,60];
        neighbor_cutoff_delta_hyperfine = 1e4;
        neighbor_cutoff_dipole_dipole = 1e3;
        neighbor_cutoff_3_spin_hahn_mod_depth = 1e-10;
        neighbor_cutoff_3_spin_hahn_taylor_4 = 1e-9;
        pulse_sequence = cp-1;
        radius = 80;
        temperature = 20;
        time_increments = [1e-9, 5e-7];
        write_auxiliary_signals = out_aux;
        write_bath = out_bath;
        write_clusters = out_clu;
        write_exchange_groups = out_ex;
        write_info = out_info;
        write_orientation_signals = out_ori;
        write_structure_pdb = out_pdb;
        write_tensors = \"tensors\";
        ").unwrap();

    let mut config = Config::new();
    for expression in expressions.iter(){
      config.parse_config_line(expression).unwrap();
    }


    let clash_distance_pbc = config.clash_distance_pbc.unwrap();
    assert!( (clash_distance_pbc - 0.1e-10).abs()/(0.1e-10) < 1e-12 );
    assert_eq!(config.cluster_batch_size,Some(20000));
    assert_eq!(config.cluster_method,Some(ClusterMethod::CCE));
    assert_eq!(config.detected_spin_position, 
        Some(DetectedSpinCoordinates::CentroidOverSerials(vec![28,29])));
    assert_eq!(config.input_structure_file, 
         Some("../../assets/TEMPO_wat_gly_70A.pdb".to_string()));
    
    assert_eq!(config.load_geometry, Some(LoadGeometry::Cube));
    assert_eq!(config.max_cluster_size, Some(4));
    assert_eq!(config.neighbor_cutoff_delta_hyperfine, Some(1e4));
    assert_eq!(config.neighbor_cutoff_dipole_dipole, Some(1e3));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_mod_depth, Some(1e-10));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_taylor_4, Some(1e-9));
    assert_eq!(config.number_timepoints, vec![40,60]);
    assert_eq!(config.pulse_sequence, Some(PulseSequence::CarrPurcell(1)));
    assert_eq!(config.radius, Some(80.0e-10));
    assert_eq!(config.temperature, Some(20.0));
    assert_eq!(config.time_increments, vec![1e-9,5e-7]);
    assert_eq!(config.write_auxiliary_signals, Some("out_aux".to_string()));
    assert_eq!(config.write_bath, Some("out_bath".to_string()));
    assert_eq!(config.write_clusters, Some("out_clu".to_string()));
    assert_eq!(config.write_exchange_groups, Some("out_ex".to_string()));
    assert_eq!(config.write_info, Some("out_info".to_string()));
    assert_eq!(config.write_orientation_signals, Some("out_ori".to_string()));
    assert_eq!(config.write_structure_pdb, Some("out_pdb".to_string()));
    assert_eq!(config.write_tensors, Some("tensors".to_string()));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_time_axis(){
    let expressions = get_tokens_from_line("\
        pulse_sequence = hahn;
        number_timepoints = [100,91];
        time_increments = [0.5e-8, 0.5e-7];
        ").unwrap();

    let mut config = Config::new();
    for expression in expressions.iter(){
      config.parse_config_line(expression).unwrap();
    }

    config.construct_time_axis().unwrap();
    let time_axis = config.get_time_axis().unwrap();
    assert_eq!(time_axis.len(),191);
    assert_eq!(time_axis[0],0.0);
    let t = 1e-8;
    assert!((time_axis[1]-t).abs()/t<1e-12);
    let t = 1e-6;
    assert!((time_axis[100]-t).abs()/t<1e-12);
    let t = 1.1e-6;
    assert!((time_axis[101]-t).abs()/t<1e-12);
    let t = 10e-6;
    assert!((time_axis[190]-t).abs()/t<1e-12);
    
  }
}

