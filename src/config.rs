use crate::clue_errors::*;

use crate::cluster::{
  partition::PartitioningMethod,
  unit_of_clustering::UnitOfClustering,
};
use crate::config::command_line_input::CommandLineInput;
use crate::config::lexer::*;
use crate::config::token::*;
use crate::config::token_algebra::*;
use crate::config::token_expressions::*;
use crate::config::particle_config::{CellType, ParticleConfig,
  EigSpecifier,TensorSpecifier};
use crate::physical_constants::{ELECTRON_G, Isotope};
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
pub mod parse_config;
pub mod parse_filter;
pub mod parse_properties;


pub const SAVE_DIR_AUXILIARY_SIGNALS: &str = "auxiliary_signals";
pub const SAVE_FILE_BATH: &str = "bath";
pub const SAVE_FILE_CLUSTERS: &str = "clusters";
pub const SAVE_DIR_INFO: &str = "info";
pub const SAVE_FILE_EXCHANGE_GROUPS: &str = "exchange_groups";
pub const SAVE_FILE_METHYL_PARTITIONS: &str = "methyl_partitions";
pub const SAVE_DIR_ORIENTATION_SIGNALS: &str = "orientations";
pub const SAVE_FILE_SANS_SPIN_SIGNALS: &str = "sans_spin_signals";
pub const SAVE_FILE_STRUCTURE_PDB: &str = "structure_pdb";
pub const SAVE_FILE_TENSORS: &str = "tensors";


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// Config contains all the setting for CluE.
#[derive(Debug,Clone,Default)]
pub struct Config{
  pub apply_pbc: Option<bool>,
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
  pub max_cell_size: Option<usize>,
  pub max_cluster_size: Option<usize>,
  pub max_spin_order: Option<usize>,
  pub min_cell_size: Option<usize>,
  pub neighbor_cutoff_coupling: Option<f64>,
  pub neighbor_cutoff_delta_hyperfine: Option<f64>,
  pub neighbor_cutoff_dipole_perpendicular: Option<f64>,
  pub neighbor_cutoff_distance: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_mod_depth: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_taylor_4: Option<f64>,
  pub number_system_instances: Option<usize>, 
  pub number_timepoints: Vec::<usize>,
  pub orientation_grid: Option<OrientationAveraging>,
  pub partitioning_method: Option<PartitioningMethod>,
  pub particles: Vec::<ParticleConfig>,
  pub pdb_model_index: Option<usize>,
  pub pulse_sequence: Option<PulseSequence>,
  pub root_dir: Option<String>, 
  pub radius: Option<f64>,
  pub rng_seed: Option<u64>,
  pub save_name: Option<String>,
  pub system_name: Option<String>,
  time_axis: Vec::<f64>,
  pub time_increments: Vec::<f64>,
  pub unit_of_clustering: Option<UnitOfClustering>,
  pub write_auxiliary_signals: Option<bool>, 
  pub write_bath: Option<bool>,
  pub write_clusters: Option<bool>,
  pub write_info: Option<bool>,
  pub write_exchange_groups: Option<bool>,
  pub write_methyl_partitions: Option<bool>,
  pub write_orientation_signals: Option<bool>,
  pub write_sans_spin_signals: Option<bool>,
  pub write_structure_pdb: Option<bool>,
  pub write_tensors: Option<bool>,
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

    if self.apply_pbc.is_none(){
      self.apply_pbc = Some(true);
    }
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
      self.load_geometry = Some(LoadGeometry::Ball);
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

    if self.root_dir.is_none(){
      self.root_dir = Some("./".to_string());
    }

    if self.unit_of_clustering.is_none(){
      self.unit_of_clustering = Some(UnitOfClustering::Set);
    }
    if self.write_bath.is_none(){
      self.write_bath = Some(true);
    }

    if self.write_info.is_none(){
      self.write_info = Some(true);
    }

    if self.write_exchange_groups.is_none(){
      self.write_exchange_groups = Some(true);
    }

    if self.write_orientation_signals.is_none(){
      self.write_orientation_signals = Some(true);
    }


    if self.write_structure_pdb.is_none(){
      self.write_structure_pdb = Some(true);
    }
    

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
      self.detected_spin_g_matrix = Some(TensorSpecifier::Eig(EigSpecifier{
        values: Some([ELECTRON_G,ELECTRON_G,ELECTRON_G]),
        x_axis: Some(VectorSpecifier::Vector(Vector3D::from([1.0, 0.0, 0.0]))),
        y_axis: Some(VectorSpecifier::Vector(Vector3D::from([0.0, 1.0, 0.0]))),
        z_axis: None,
          }));
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
    io::write_data(&[self.time_axis.clone()], 
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
  Ball,
}

/// `DetectedSpinCoordinates` lists options for indicating the coordinates
/// of the detected spin.
#[derive(Debug,Clone,PartialEq)]
pub enum DetectedSpinCoordinates{
  CentroidOverGroup(Vec::<String>),
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
//impl ClusterMethod{
  //pub fn from(cluster_method_str: String) => Self{
  //}
//}

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
      }

    }


    Ok(())
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;

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

