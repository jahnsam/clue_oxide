use crate::clue_errors::*;

use crate::config::command_line_input::CommandLineInput;
use crate::config::lexer::*;
use crate::config::token::*;
use crate::config::token_algebra::*;
use crate::config::token_expressions::*;
use crate::config::particle_config::{CellType, ParticleConfig,TensorSpecifier};

use crate::physical_constants::{ANGSTROM,ELECTRON_G, Isotope};
use crate::space_3d::Vector3D;
use crate::structure::particle_filter::VectorSpecifier;
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
  pub clash_distance: Option<f64>, 
  pub clash_distance_pbc: Option<f64>,
  pub cluster_batch_size: Option<usize>, 
  pub cluster_method: Option<ClusterMethod>,
  pub clusters_file: Option<String>,
  pub density_matrix: Option<DensityMatrixMethod>,
  pub detected_spin_g_matrix: Option<TensorSpecifier>,
  pub detected_spin_identity: Option<Isotope>,
  pub detected_spin_multiplicity: Option<usize>,
  pub detected_spin_position: Option<DetectedSpinCoordinates>,
  pub detected_spin_transition: Option<[usize;2]>,
  pub extracell_particles: Vec::<ParticleConfig>,
  pub input_structure_file: Option<String>,
  pub load_geometry: Option<LoadGeometry>,
  pub magnetic_field: Option<Vector3D>,
  pub max_cluster_size: Option<usize>,
  pub neighbor_cutoff_coupling: Option<f64>,
  pub neighbor_cutoff_delta_hyperfine: Option<f64>,
  pub neighbor_cutoff_dipole_perpendicular: Option<f64>,
  pub neighbor_cutoff_distance: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_mod_depth: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_taylor_4: Option<f64>,
  pub number_system_instances: Option<usize>, 
  pub number_timepoints: Vec::<usize>,
  pub orientation_grid: Option<OrientationAveraging>,
  pub particles: Vec::<ParticleConfig>,
  pub pdb_model_index: Option<usize>,
  pub pulse_sequence: Option<PulseSequence>,
  pub remove_partial_methyls: Option<bool>,
  pub root_dir: Option<String>, 
  pub radius: Option<f64>,
  pub rng_seed: Option<u64>,
  pub save_name: Option<String>,
  pub system_name: Option<String>,
  time_axis: Vec::<f64>,
  pub time_increments: Vec::<f64>,
  pub write_auxiliary_signals: Option<String>, 
  pub write_bath: Option<String>,
  pub write_clusters: Option<String>,
  pub write_info: Option<String>,
  pub write_exchange_groups: Option<String>,
  pub write_methyl_partitions: Option<String>,
  pub write_orientation_signals: Option<String>,
  pub write_sans_spin_signals: Option<String>,
  pub write_structure_pdb: Option<String>,
  pub write_tensors: Option<String>,
}


impl Config{
  //----------------------------------------------------------------------------
  /// This function generates a empty `Config`.
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------
  /// This function generates a `Config` from an input string.
  pub fn from(input: &str) -> Result<Self,CluEError>
  {
    let expressions
      = match get_tokens_from_line(input){
        Ok(exps) => exps,
        Err(err) => return Err(err),
      };

    let mut config = Config::new();

    match  config.parse_token_stream(expressions){
      Ok(_) => (),
      Err(err) => return Err(err),
    }

    config.set_extracell_particles();

    config.set_defaults()?;

    match config.construct_time_axis(){
      Ok(_) => (),
      Err(err) => return Err(err),
    }

    Ok(config)
  }
  //----------------------------------------------------------------------------
  /// This function sets config setting that are necessary to be `Some`,
  /// but will likely be the same for most simulations.  
  /// Any field that is already `Some` will be left alone.
  pub fn set_defaults(&mut self) -> Result<(),CluEError> {

    if self.clash_distance.is_none(){
      self.clash_distance = Some(1e-12);
    }
    //if self.clash_distance_pbc.is_none(){
    //  self.clash_distance_pbc = Some(1e-12);
    //}


    if self.cluster_batch_size.is_none(){
      self.cluster_batch_size = Some(1000);
    }
    if self.load_geometry.is_none(){
      self.load_geometry = Some(LoadGeometry::Sphere);
    }
    if self.pdb_model_index.is_none(){
      self.pdb_model_index = Some(0);
    }
    if self.density_matrix.is_none(){
      self.density_matrix = Some(DensityMatrixMethod::Identity);
    }

    if self.number_system_instances.is_none(){
      self.number_system_instances = Some(1);
    }
    if self.remove_partial_methyls.is_none(){
      self.remove_partial_methyls = Some(false);
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


    let set_default_write_path = 
      |write_opt: &mut Option<String>,default: Option<String>|
    {
      match write_opt{
        None => *write_opt = default,
        Some(path) => if path.is_empty(){ *write_opt = None;},
      }

    };

    set_default_write_path(&mut self.write_auxiliary_signals,
        None );

    set_default_write_path(&mut self.write_bath,
        Some("bath".to_string()) );

    set_default_write_path(&mut self.write_clusters,
        None );

    set_default_write_path(&mut self.write_info,
        Some("info".to_string()) );

    set_default_write_path(&mut self.write_exchange_groups,
        Some("exchange_groups".to_string()) );

    set_default_write_path(&mut self.write_methyl_partitions,
        None );

    set_default_write_path(&mut self.write_orientation_signals,
        Some("orientations".to_string()) );

    set_default_write_path(&mut self.write_sans_spin_signals,
        None );

    set_default_write_path(&mut self.write_structure_pdb,
        Some("structure_pdb".to_string()) );

    // The tensors file is large, so do not write it by default.
    
    set_default_write_path(&mut self.write_tensors,
        None);
    

    self.set_detected_spin()
  }
  //----------------------------------------------------------------------------
  // This function sets the detected spin.
  fn set_detected_spin(&mut self) -> Result<(),CluEError>{

    if self.detected_spin_identity.is_none(){
      self.detected_spin_identity = Some(Isotope::Electron);
    }

    // Set g-matrix
    if self.detected_spin_g_matrix.is_some(){
     
    }else{
      self.detected_spin_g_matrix = Some(TensorSpecifier{
        values: Some([ELECTRON_G,ELECTRON_G,ELECTRON_G]),
        x_axis: Some(VectorSpecifier::Vector(Vector3D::from([1.0, 0.0, 0.0]))),
        y_axis: Some(VectorSpecifier::Vector(Vector3D::from([0.0, 1.0, 0.0]))),
        z_axis: None,
          });
    }
    if self.detected_spin_multiplicity.is_none(){
      self.detected_spin_multiplicity = match self.detected_spin_identity{
        Some(isotope) => Some(isotope.spin_multiplicity()),
        None => return Err(CluEError::NoDetectedSpinIdentity),
      };
    }
    if self.detected_spin_transition.is_none(){
      self.detected_spin_transition = Some([0,1]);
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function the highest spin multiplicity for a "#[group(...)]". 
  pub fn max_spin_multiplicity_for_particle_config(&self, id: usize) 
    -> Option<usize>
  {


   let mult0 = if id < self.particles.len() {
     self.particles[id].max_possible_spin_multiplicity()
   } else{
     None
   };

   let mult1 = if id < self.extracell_particles.len() {
     self.extracell_particles[id].max_possible_spin_multiplicity()
   } else{
     None
   };

   match (mult0,mult1){
     (Some(m0),Some(m1)) => { if m0 >= m1 {Some(m0)}else{Some(m1)}  },
     (Some(m0),None) => Some(m0),
     (None,Some(m1)) => Some(m1),
     (None,None) => None,
   }

  }
  //----------------------------------------------------------------------------
  /// This function gets the simulated experiment's time axis.
  pub fn get_time_axis(&self) -> Result<&Vec::<f64>,CluEError>
  {
    if self.time_axis.is_empty() {
      return Err(CluEError::NoTimeAxis);
    };
    Ok(&self.time_axis)
  }
  //----------------------------------------------------------------------------
  /// This function constructs the simulated experiment's time axis.
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
  /// This function writes the simulated experiment's time axis to a file.
  pub fn write_time_axis(&self,save_path: String) -> Result<(),CluEError>{
    io::write_data(& vec![self.time_axis.clone()], 
        &format!("{}/time_axis.csv",save_path), vec!["time_axis".to_string()])
  }
  //----------------------------------------------------------------------------
  /// This function finds the `ParticleConfig` that has `label` if it exists.
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
/// `DensityMatrixMethod` specifies different methods for determining the
// /density ,atrix.
#[derive(Debug,Clone,PartialEq)]
pub enum DensityMatrixMethod{
  ApproxThermal(f64),
  Identity,
  Thermal(f64),
}

/// `LoadGeometry` specifies how the system is trimmed/
#[derive(Debug,Clone,PartialEq)]
pub enum LoadGeometry{
  Cube,
  Sphere,
}

/// `DetectedSpinCoordinates` lists options for indicating the coordinates
/// of the detected spin.
#[derive(Debug,Clone,PartialEq)]
pub enum DetectedSpinCoordinates{
  CentroidOverSerials(Vec::<u32>),
  XYZ(Vector3D),
  ProbabilityDistribution(IntegrationGrid),
}

/// `ClusterMethod` list different cluster simulation methods.
#[derive(Debug,Clone,PartialEq)]
pub enum ClusterMethod{
  AnalyticRestricted2CCE,
  CCE,
  GCCE,
}

/// `PulseSequence` lists the options for pulse sequences to simulate.
#[derive(Debug,Clone,PartialEq)]
pub enum PulseSequence{
  CarrPurcell(usize),
  //RefocusedEcho,
}

/// `OrientationAveraging` lists the different orientation averaging options.
#[derive(Debug,Clone,PartialEq)]
pub enum OrientationAveraging{
  Grid(IntegrationGrid),
  Lebedev(usize),
  Random(usize),
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Config{

  /// This function takes the parsed command line input and develops the main
  /// configuration data structure.
  pub fn read_input(input: CommandLineInput) -> Result<Self,CluEError>{
    let Some(filename) = &input.config_file else{
      return Err(CluEError::NoInputFile);
    };
  
    let mut config = Config::new();

    let mut token_stream = get_tokens_from_file(filename)?;

    if let Some(options) = &input.config_options{
      let mut expressions = get_tokens_from_line(options)?;
      token_stream.append(&mut expressions);
    }

    config.parse_token_stream(token_stream)?;

    config.set_extracell_particles();

    Ok(config)
  }
  //----------------------------------------------------------------------------
  // This functions sets properties for particles outside the primary cell.
  fn set_extracell_particles(&mut self){
    
    let mut n_extra = 0; 
    let mut n_primary = 0; 
    for particle_config in self.particles.iter(){
      match particle_config.cell_type{
      CellType::AllCells => {
        n_extra += 1;
        n_primary += 1;
      },
      CellType::PrimaryCell => n_primary += 1,
      CellType::Extracells => n_extra += 1,
      }
    }

    let mut particles = Vec::<ParticleConfig>::with_capacity(n_primary);
    self.extracell_particles = Vec::<ParticleConfig>::with_capacity(n_extra);

    for particle_config in self.particles.iter(){
      match particle_config.cell_type{
      CellType::AllCells => {
        self.extracell_particles.push(particle_config.clone());
        particles.push(particle_config.clone());
      },
      CellType::PrimaryCell => particles.push(particle_config.clone()),
      CellType::Extracells 
        => self.extracell_particles.push(particle_config.clone()),
      }
    }

    self.particles = particles;

  }
  //----------------------------------------------------------------------------
  /// This function updates config from a `Vec<TokenExpression>`.
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
        ConfigMode::Filter => {
          let Some(_label) = &mode.label else{
            return Err(CluEError::FilterNeedsALabel);
          };
          self.parse_filter_line(expression,&mode.label)?;
        },
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
      Token::ClusterDensityMatrix => {
        if let Some(_value) = &self.density_matrix{
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
          Token::ApproxThermal | Token::Thermal => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let new_expression = TokenExpression{
              lhs: expression.lhs.clone(),
              rhs: Some(args),
              relationship: Some(Token::Equals),
              line_number: expression.line_number
            };

            let mut temperature_opt: Option<f64> = None;
            set_to_some_f64(&mut temperature_opt,&new_expression)?;
            if let Some(temperature) = temperature_opt{
              match rhs[0]{
                Token::ApproxThermal => self.density_matrix 
                  = Some(DensityMatrixMethod::ApproxThermal(temperature)),
                Token::Thermal => self.density_matrix 
                  = Some(DensityMatrixMethod::Thermal(temperature)),
                _ => return Err(CluEError::InvalidToken(expression.line_number,
                  expression.lhs[0].to_string())),
              }
            }else{
              return Err(CluEError::NoDensityMatrixMethod);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::Identity => self.density_matrix 
            = Some(DensityMatrixMethod::Identity),
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
                rhs[0].to_string())),
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        }
      },
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
          Token::GCCE => self.cluster_method = Some(ClusterMethod::GCCE),
          Token::R2CCE => self.cluster_method 
            = Some(ClusterMethod::AnalyticRestricted2CCE),
          _  => return Err(CluEError::InvalidArgument(expression.line_number,
                "valid cluster method are \"CCE\" and \"r2CCE\"".to_string())),
        }

      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClustersFile 
        => set_to_some_string(&mut self.clusters_file, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::DetectedSpinGMatrix 
        | Token::DetectedSpinGX| Token::DetectedSpinGY | Token::DetectedSpinGZ
        =>{
          if self.detected_spin_g_matrix.is_none() {
            self.detected_spin_g_matrix = Some(TensorSpecifier::new())
          }
          let Some(g_matrix) = &mut self.detected_spin_g_matrix else{
            return Err(CluEError::NoGMatrixSpecifier );
          };
          let Some(rhs) = &expression.rhs else{
            return Err(CluEError::NoRHS(expression.line_number));
          }; 

          match expression.lhs[0]{
            Token::DetectedSpinGMatrix => {
              let mut g_values = Vec::<f64>::new();
        
              set_to_vec_f64(&mut g_values,expression)?;
            
              if g_matrix.values.is_some(){
                return Err(already_set());
              }

              if g_values.len() == 1{
                g_matrix.values = Some([g_values[0],g_values[0],g_values[0]]);
              }else if g_values.len() == 3{
                g_matrix.values = Some([g_values[0],g_values[1],g_values[2]]);
              }else{
                return Err(CluEError::CannotInferEigenvalues(
                    expression.line_number));
              }
            },
            Token::DetectedSpinGX 
              => set_to_some_vector_specifier(&mut g_matrix.x_axis, expression,
                  "detected_spin")?,
            Token::DetectedSpinGY 
              => set_to_some_vector_specifier(&mut g_matrix.y_axis, expression,
                  "detected_spin")?,
            Token::DetectedSpinGZ
              => set_to_some_vector_specifier(&mut g_matrix.z_axis, expression,
                  "detected_spin")?,
            _ => return Err(CluEError::InvalidToken(expression.line_number,
                  rhs[0].to_string())),
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

            let grid = IntegrationGrid::read_from_csv(&(args[0].to_string()))?
              .scale(ANGSTROM);

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
            if let Some(r_angstrom) = r_opt{
              let r = r_angstrom.scale(ANGSTROM); 
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
      Token::NeighborCutoffCoupling 
        => set_to_some_f64(&mut self.neighbor_cutoff_coupling,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffDeltaHyperfine 
        => set_to_some_f64(&mut self.neighbor_cutoff_delta_hyperfine,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffDipolePerpendicular
        => set_to_some_f64(&mut self.neighbor_cutoff_dipole_perpendicular,
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
      Token::NumberSystemInstances
        => set_to_some_usize(&mut self.number_system_instances, expression)?,
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
          Token::Lebedev | Token::Random => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let new_expression = TokenExpression{
              lhs: expression.lhs.clone(),
              rhs: Some(args),
              relationship: Some(Token::Equals),
              line_number: expression.line_number
            };

            let mut n_grid_opt: Option<usize> = None;
            set_to_some_usize(&mut n_grid_opt,&new_expression)?;
            if let Some(n_grid) = n_grid_opt{
              match rhs[0]{
                Token::Lebedev => self.orientation_grid 
                  = Some(OrientationAveraging::Lebedev(n_grid)),
                Token::Random => self.orientation_grid 
                  = Some(OrientationAveraging::Random(n_grid)),
                _ => return Err(CluEError::InvalidToken(expression.line_number,
                  expression.lhs[0].to_string())),
              }
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
      Token::SaveDir 
        => set_to_some_string(&mut self.save_name, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::SystemName 
        => set_to_some_string(&mut self.system_name, expression)?,
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
        => {
          //set_to_some_f64(&mut self.temperature, expression)?
          return Err(CluEError::DeprecatedKeywordReplaced(
                expression.line_number,
                "temperature = T;".to_string(),
                "cluster_density_matrix = thermal(T);".to_string(),
                ))
        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::TimeIncrements 
        => set_to_vec_f64(&mut self.time_increments, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteAuxiliarySignals
        //=> set_to_some_string(&mut self.write_auxiliary_signals, expression)?,
        =>{
          if let Some(_value) = &self.write_auxiliary_signals{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_auxiliary_signals = match do_write_opt{
            Some(true) => Some("auxiliary_signals".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteBath
        //=> set_to_some_string(&mut self.write_bath, expression)?,
        =>{
          if let Some(_value) = &self.write_bath{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_bath = match do_write_opt{
            Some(true) => Some("bath".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteClusters
        //=> set_to_some_string(&mut self.write_clusters, expression)?,
        =>{
          if let Some(_value) = &self.write_clusters{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_clusters = match do_write_opt{
            Some(true) => Some("clusters".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteExchangeGroups
        //=> set_to_some_string(&mut self.write_exchange_groups, expression)?,
        =>{
          if let Some(_value) = &self.write_exchange_groups{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_exchange_groups = match do_write_opt{
            Some(true) => Some("exchange_groups".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteInfo
        //=> set_to_some_string(&mut self.write_info, expression)?,
        =>{
          if let Some(_value) = &self.write_info{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_info = match do_write_opt{
            Some(true) => Some("info".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteMethylPartitions => { 
          if let Some(_value) = &self.write_methyl_partitions{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_methyl_partitions = match do_write_opt{
            Some(true) => Some("methyl_partitions".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteOrientationSignals
        //=> set_to_some_string(&mut self.write_orientation_signals, expression)?,
        =>{
          if let Some(_value) = &self.write_orientation_signals{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_orientation_signals = match do_write_opt{
            Some(true) => Some("orientations".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteSansSpinSignals => { 
          if let Some(_value) = &self.write_sans_spin_signals{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_sans_spin_signals = match do_write_opt{
            Some(true) => Some("sans_spin_signals".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteStructurePDB 
        //=> set_to_some_string(&mut self.write_structure_pdb, expression)?,
        =>{
          if let Some(_value) = &self.write_structure_pdb{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_structure_pdb = match do_write_opt{
            Some(true) => Some("structure_pdb".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteTensors 
        //=> set_to_some_string(&mut self.write_tensors, expression)?,
        =>{
          if let Some(_value) = &self.write_tensors{
            return Err(already_set());
          }

          let mut do_write_opt: Option<bool> = None;
          set_to_some_bool(&mut do_write_opt ,expression)?;
          self.write_tensors = match do_write_opt{
            Some(true) => Some("tensors".to_string()),
            Some(false) => Some("".to_string()),
            None => return Err(
                CluEError::ExpectedBoolRHS(expression.line_number)),
          };

        },
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
  use crate::structure::particle_filter::{
    SecondaryParticleFilter,VectorSpecifier};

  #[test]
  fn test_parse_config_line(){
    let expressions = get_tokens_from_line("\
        #[config]\n
        clash_distance_pbc = 0.1;
        cluster_batch_size = 20000;
        cluster_density_matrix = thermal(20);
        cluster_method = cce;
        detected_spin_g_matrix = [2.0097, 2.0064, 2.0025];
        detected_spin_g_y = [-1.1500, -0.4700, 0.7100];
        detected_spin_g_x = diff(group(tempo_c1) , filter(tempo_c19) );
        detected_spin_position = centroid_over_serials([28,29]);
        input_clusters_file = \"clusters_file.txt\";
        input_structure_file = \"../../assets/TEMPO_wat_gly_70A.pdb\";
        load_geometry = cube;
        max_cluster_size = 4;
        number_timepoints = [40,60];
        neighbor_cutoff_delta_hyperfine = 1e4;
        neighbor_cutoff_coupling = 1e3;
        neighbor_cutoff_dipole_perpendicular = 100;
        neighbor_cutoff_3_spin_hahn_mod_depth = 1e-10;
        neighbor_cutoff_3_spin_hahn_taylor_4 = 1e-9;
        orientation_grid = lebedev(170);
        pulse_sequence = cp-1;
        radius = 80;
        save_dir = \"save_directory\";
        time_increments = [1e-9, 5e-7];
        write_auxiliary_signals = true;
        write_bath = true;
        write_clusters = true;
        write_exchange_groups = true;
        write_info = true;
        write_methyl_partitions = true;
        write_orientation_signals = true;
        write_sans_spin_signals = true;
        write_structure_pdb = true;
        write_tensors = false;
        ").unwrap();

    let mut config = Config::new();
    config.parse_token_stream(expressions).unwrap();
    config.set_defaults().unwrap();

    let clash_distance_pbc = config.clash_distance_pbc.unwrap();
    assert!( (clash_distance_pbc - 0.1e-10).abs()/(0.1e-10) < 1e-12 );
    assert_eq!(config.cluster_batch_size,Some(20000));
    assert_eq!(config.density_matrix, Some(DensityMatrixMethod::Thermal(20.0)));
    assert_eq!(config.cluster_method,Some(ClusterMethod::CCE));

    let g_matrix = config.detected_spin_g_matrix.unwrap();
    assert_eq!(g_matrix.values, 
        Some([2.0097, 2.0064, 2.0025] ) );
    assert_eq!(g_matrix.z_axis, None);
    assert_eq!(g_matrix.y_axis, 
        Some(VectorSpecifier::Vector(
            Vector3D::from([-1.1500, -0.4700, 0.7100]) ) ) );
    assert_eq!(g_matrix.x_axis, 
        Some(VectorSpecifier::Diff(
            SecondaryParticleFilter::Filter, "tempo_c1".to_string(),
            SecondaryParticleFilter::Filter, "tempo_c19".to_string(),
            )));

    assert_eq!(config.detected_spin_identity, Some(Isotope::Electron));
    assert_eq!(config.detected_spin_multiplicity, Some(2));
    assert_eq!(config.detected_spin_position, 
        Some(DetectedSpinCoordinates::CentroidOverSerials(vec![28,29])));
    assert_eq!(config.input_structure_file, 
         Some("../../assets/TEMPO_wat_gly_70A.pdb".to_string()));
    assert_eq!(config.clusters_file, 
         Some("clusters_file.txt".to_string()));
    
    assert_eq!(config.load_geometry, Some(LoadGeometry::Cube));
    assert_eq!(config.max_cluster_size, Some(4));
    assert_eq!(config.neighbor_cutoff_delta_hyperfine, Some(1e4));
    assert_eq!(config.neighbor_cutoff_coupling, Some(1e3));
    assert_eq!(config.neighbor_cutoff_dipole_perpendicular, Some(1e2));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_mod_depth, Some(1e-10));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_taylor_4, Some(1e-9));
    assert_eq!(config.number_timepoints, vec![40,60]);
    assert_eq!(config.pulse_sequence, Some(PulseSequence::CarrPurcell(1)));
    assert_eq!(config.radius, Some(80.0e-10));
    assert_eq!(config.save_name, Some(String::from("save_directory")));
    assert_eq!(config.time_increments, vec![1e-9,5e-7]);
    assert_eq!(config.write_auxiliary_signals, 
        Some("auxiliary_signals".to_string()));
    assert_eq!(config.write_bath, Some("bath".to_string()));
    assert_eq!(config.write_clusters, Some("clusters".to_string()));
    assert_eq!(config.write_exchange_groups, 
        Some("exchange_groups".to_string()));
    assert_eq!(config.write_info, Some("info".to_string()));
    assert_eq!(config.write_methyl_partitions, 
        Some("methyl_partitions".to_string()));
    assert_eq!(config.write_orientation_signals, 
        Some("orientations".to_string()));
    assert_eq!(config.write_structure_pdb, Some("structure_pdb".to_string()));
    assert_eq!(config.write_sans_spin_signals, 
        Some("sans_spin_signals".to_string()));
    assert_eq!(config.write_tensors, None);
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

