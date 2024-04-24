//pub mod cif;
pub mod pdb;
pub mod extended_structure;
pub mod primary_structure;
pub mod exchange_groups;
pub mod particle;
pub mod particle_filter;

//use crate::config::particle_config;
use crate::config::{
  Config, 
  particle_config::{ParticleConfig,ParticleProperties,IsotopeProperties,
    TensorSpecifier},
};
use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::integration_grid::IntegrationGrid;
use crate::math;
use crate::space_3d::Vector3D;
use crate::structure::{
  exchange_groups::ExchangeGroupManager,
  particle::Particle,
  particle_filter::*
};
use crate::physical_constants::{ANGSTROM,Isotope};
use crate::space_3d::SymmetricTensor3D;

use rand_chacha::ChaCha20Rng;
use std::fs::File;
use std::io::{BufWriter, Write};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// TODO: remove in favor of config properties.
#[derive(Debug,Clone)]
pub struct DetectedSpin{
  pub gamma_matrix: SymmetricTensor3D,
  pub isotope: Isotope,
  pub weighted_coordinates: IntegrationGrid,
  pub spin_multiplicity: usize,
  pub transition: [usize;2],
}
//------------------------------------------------------------------------------
impl DetectedSpin{
  pub fn gyromagnetic_ratio_matrix(&self) 
    -> &SymmetricTensor3D
  {

    &self.gamma_matrix
  }
  //----------------------------------------------------------------------------
  pub fn spin_multiplicity(&self) -> usize{
    self.spin_multiplicity
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct Structure{
  // the particle giving rise to the detected signal
  pub detected_particle: Option<DetectedSpin>,
  
  // all particles that are not directly detected   
  pub bath_particles: Vec::<Particle>,

  // a list of which bath_particles have spin
  bath_spins_indices: Vec::<usize>,

  // set of lists of all periodic copies of each particle from the primary cell 
  pub cell_indices: Vec::<Vec::<Option<usize>>>,

  // set of displacements for periodic boundary conditions
  pub cell_offsets: Vec::<Vector3D>,

  // list of chemical bonds
  pub connections: AdjacencyList,

  // set of molecule
  molecules: Vec::<Vec::<usize>>,

  // a list specifing which molecule each particle belongs to
  molecule_ids: Vec::<usize>,

  // collections of atoms that should be isotope substituted together
  cosubstitution_groups: Vec::< Vec::<usize> >,

  // information on methyls
  pub exchange_groups: Option<ExchangeGroupManager>,

  // indices for which config.particle corresponds to each particle
  particle_config_ids: Vec::<Option<usize>>,

  // indices for which config.extracell_particle corresponds to each particle
  extracell_particle_config_ids: Vec::<Option<usize>>,

  // origin of the pdb frame
  pub pdb_origin: Vector3D,

  // list indices indicating the particle that each particle is a periodic
  // boundary condition copy of  
  primary_cell_indices: Vec::<usize>,

  // indices for display and output file
  nth_active_to_reference_index: Vec::<usize>,
}

impl Structure{
  pub fn new(
      bath_particles: Vec::<Particle>,
      connections: AdjacencyList,
      cell_offsets: Vec::<Vector3D>
      ) -> Self
  {
    Structure{
      detected_particle: None,
      bath_particles,
      bath_spins_indices: Vec::<usize>::new(),
      cell_indices: Vec::<Vec::<Option<usize>>>::new(),
      cell_offsets,
      connections,
      extracell_particle_config_ids: Vec::<Option<usize>>::new(),
      molecules: Vec::<Vec::<usize>>::new(),
      molecule_ids: Vec::<usize>::new(),
      cosubstitution_groups: Vec::<Vec::<usize>>::new(),
      exchange_groups: None,
      particle_config_ids: Vec::<Option<usize>>::new(),
      pdb_origin: Vector3D::zeros(),
      primary_cell_indices: Vec::<usize>::new(),
      nth_active_to_reference_index: Vec::<usize>::new(),
    }

  }
  //----------------------------------------------------------------------------
  pub fn build_structure(rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<Self,CluEError>
  {

    let Some(filename) = &config.input_structure_file else{
      return Err(CluEError::NoStructureFile);
    };

    let Some(model_idx) = config.pdb_model_index else{
      return Err(CluEError::NoModelIndex);
    };

    let mut structure = pdb::parse_pdb(filename,model_idx)?;

    structure.build_primary_structure(config)?;

    structure.build_extended_structure(rng, config)?;

    structure.map_nth_active_to_reference_indices();

    Ok(structure)
  }
  //----------------------------------------------------------------------------
  pub fn number(&self) -> usize{
    self.bath_particles.len()
  }
  //----------------------------------------------------------------------------
  pub fn number_active(&self) -> usize{
    let mut n_active = 0;
    for particle in self.bath_particles.iter(){
      if particle.active{ n_active += 1;}
    }
    n_active
  }
  //----------------------------------------------------------------------------
  pub fn find<'a>(&'a self, particle_filter: &ParticleFilter)
    -> Vec::<&'a Particle>
  {
    let indices = particle_filter.filter(self);
    let mut out = Vec::<&Particle>::with_capacity(indices.len());

    for idx in indices{
      out.push(&self.bath_particles[idx]);
    }

    out
  }
  //----------------------------------------------------------------------------
  /// This function identifies which pbc copy of the primary cell
  /// `particle_particle` with index `idx` resides in.
  pub fn cell_id(&self, idx: usize) -> Result<usize,CluEError>
  {
    if idx >= self.primary_cell_indices.len(){
      return Err(CluEError::CannotFindCellID(idx));
    }

    let idx0 = self.primary_cell_indices[idx];
    let cell_indices = &self.cell_indices[idx0];
    let Some(cell_id) = cell_indices.iter().position(|&x| {
        match x{ Some(id) => id == idx, None => false,} }) else{
      return Err(CluEError::CannotFindCellID(idx));
    };

    Ok(cell_id)
  }
  //----------------------------------------------------------------------------
  pub fn molecule_id(&self,idx: usize) -> usize{

    let idx0 = self.primary_cell_indices[idx];
   
    self.molecule_ids[idx0]

  }
  //----------------------------------------------------------------------------
  /// This function retrieves the reference index of the bath particle with
  /// index `bath_index`.
  /// The reference index of the the detected spin in zero, and the first
  /// bath particle has index 1.
  pub fn get_reference_index_from_bath_index(&self, bath_index: usize) 
    -> Result<usize,CluEError>
  {
    if bath_index >= self.bath_particles.len(){
      return Err(CluEError::CannotFindRefIndexFromBathIndex(bath_index));
    }
    Ok(bath_index + 1)
  }
  //----------------------------------------------------------------------------
  /// This function retrieves the bath index of `n`th active particle.
  pub fn get_bath_index_of_nth_active(&self, n: usize) 
    -> Result<usize,CluEError>
  {
    if n >= self.nth_active_to_reference_index.len(){
      return Err(CluEError::CannotFindRefIndexFromNthActive(n));
    }
    Ok(self.nth_active_to_reference_index[n] - 1)
  }
  //----------------------------------------------------------------------------
  /// This function retrieves the reference index of `n`th active particle.
  /// The reference index of the the detected spin in zero, and the first
  /// bath particle has index 1.
  /// The detected spin is always active and corresponds to `n=0`.
  pub fn get_reference_index_of_nth_active(&self, n: usize)
    -> Result<usize,CluEError>
  {
    if n >= self.nth_active_to_reference_index.len(){
      return Err(CluEError::CannotFindRefIndexFromNthActive(n));
    }
    Ok(self.nth_active_to_reference_index[n])

  }
  //----------------------------------------------------------------------------
  /// This function builds `nth_active_to_reference_index`, which maps
  /// the set of active particles to the set of all system particles.
  pub fn map_nth_active_to_reference_indices(&mut self){

    let mut act_to_ref = Vec::<usize>::with_capacity(self.number_active() + 1);

    act_to_ref.push(0);

    for (bath_index,particle) in self.bath_particles.iter().enumerate(){
      if particle.active{ 
        act_to_ref.push(bath_index + 1);
      }
    }

    self.nth_active_to_reference_index = act_to_ref;
  }
  //----------------------------------------------------------------------------
  /// This function takes a reference index and returns either `Ok(n)`,
  /// where when counting active particles in bath_particles, particle
  /// reference_index is the nth active particle.
  pub fn get_nth_active_from_reference_index(&self, reference_index: usize) 
    -> Result<usize,CluEError>
  { 
    // Check that the particle is a bath particle.
    if reference_index == 0{
      return Err(CluEError::DetectedSpinDoesNotHaveAnActiveIndex);
    }
    let bath_index = reference_index - 1;

    if bath_index >= self.bath_particles.len(){
      return Err(CluEError::CannotFindParticleForRefIndex(reference_index));
    }

    if !self.bath_particles[bath_index].active{
      return Err(CluEError::ParticleIsNotActive(reference_index));
    }

    // Count active particle upto and including the bath_index particle.
    let mut nth_active = 0;

    for (ii,particle) in self.bath_particles.iter().enumerate(){
      if particle.active{
        nth_active += 1;
      }
      if ii == bath_index{
        break;
      }
    }
    /*
    let nth_active: usize = self.bath_particle.iter().take(bath_index+1)
      .map(|particle| 
          if particle.active{
            1
          }else{
            0
          }
       ).sum();
    */
    let ref_idx = self.get_reference_index_of_nth_active(nth_active)?;
    if ref_idx != reference_index {
      return Err(CluEError::CannotFindParticleForRefIndex(reference_index));
    }

    Ok(nth_active)
  }
  //----------------------------------------------------------------------------
  /// This function takes a list of PDB serial numbers an returns the 
  /// centroid position.
  pub fn centroid_over_serials(&self, serials: Vec::<u32>) 
    -> Result<Vector3D,CluEError>
  {
    let mut filter = ParticleFilter::new();
    filter.serials = serials;
    filter.indices = math::unique(self.primary_cell_indices.clone());

    let indices = filter.filter(self);
    let mut r_ave = Vector3D::zeros();
    for &idx in indices.iter(){
      let r = &self.bath_particles[idx].coordinates;
      r_ave = &r_ave + r;
    }
    r_ave = r_ave.scale(1.0/(indices.len() as f64));

    Ok(r_ave)
  }
  //----------------------------------------------------------------------------
  /*
  fn get_extended_index(&self, 
      primary_cell_index: usize, unit_cell_id: usize, guess_idx: mut usize) 
    -> Option<usize>
  {
    if unit_cell_id == 0{
      return Some(primary_cell_index);
    }

    let decend: bool
    let idx = self.primary_cell_indices[guess_idx];
    if idx0 == primary_cell_index && 

  }
  */
  //----------------------------------------------------------------------------
  /*
  pub fn molecule(&self,mol_id: usize) -> Vec::<usize>{
    let mol_id0 = mol_id%self.molecules.len();
    let cell_id = (mol_id-mol_id0)/self.molecules.len();

    let mut molecule = self.molecules[mol_id0].clone();

    for idx in molecule.iter_mut(){
      idx = 
    }

  }
  */
  //----------------------------------------------------------------------------
  pub fn extract_exchange_coupling(&self, particle_index: usize,
      config: &Config)
    -> f64
  {
    let Some(isotope_properties) = self.extract_isotope_properties(
        particle_index,config) 
    else{
      return 0.0;
    };

    if let Some(ex_coup) = isotope_properties.exchange_coupling {
      ex_coup
    }else{
      0.0
    }
  }
  //----------------------------------------------------------------------------
  pub fn extract_electric_quadrupole_specifier<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a TensorSpecifier>
  {
    let Some(isotope_properties) = self.extract_isotope_properties(
        particle_index,config) 
    else{
      return None;
    };

    isotope_properties.electric_quadrupole_coupling.as_ref()
  }
  //----------------------------------------------------------------------------
  pub fn extract_g_matrix_specifier<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a TensorSpecifier>
  {
    let Some(isotope_properties) = self.extract_isotope_properties(
        particle_index,config) 
    else{
      return None;
    };

    isotope_properties.g_matrix.as_ref()
  }
  //----------------------------------------------------------------------------
  pub fn extract_hyperfine_specifier<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a TensorSpecifier>
  {
    let Some(isotope_properties) = self.extract_isotope_properties(
        particle_index,config) 
    else{
      return None;
    };

    isotope_properties.hyperfine_coupling.as_ref()
  }
  //----------------------------------------------------------------------------
  pub fn extract_isotope_properties<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a IsotopeProperties>
  {
    let particle_index_0 = self.primary_cell_indices[particle_index];

    let Some(p_cfg_id) = self.particle_config_ids[particle_index_0] else{
      return None;
    };

    let Some(properties) = &config.particles[p_cfg_id].properties else{
      return None;
    };

    let key = self.bath_particles[particle_index].isotope.to_string();
    
    properties.isotope_properties.get(&key)

  }
  //----------------------------------------------------------------------------
  // This function searches the user specified particle configuration (if any),
  // and determines which one each particle goes with.
  // This function returns a error if multiple filter cover by the same 
  // particle. 
  fn pair_particle_configs(&self, particle_configs: &[ParticleConfig]) 
    -> Result< Vec::<Option<usize>>,CluEError>
  {
  
    // Initialize pairings.
    let n_particles = self.bath_particles.len();

    let mut particle_config_ids = (0..n_particles).map(|_| None).
      collect::<Vec::<Option<usize>>>();

    // Pair particles with configs.
    for (id,particle_config) in particle_configs.iter().enumerate(){
      if let Some(filter) = &particle_config.filter{
        let indices = filter.filter(self);

        for idx in indices{
          if let Some(prev_id) = particle_config_ids[idx]{
            return Err(CluEError::FiltersOverlap(
              particle_configs[prev_id].label.to_string(),
              particle_config.label.to_string()
            ));
          }
          particle_config_ids[idx] = Some(id);
        }
      
      } 
    }
    
    Ok(particle_config_ids)
  }
  //----------------------------------------------------------------------------
  pub fn bath_to_csv(&self, filename: &str, config: &Config) 
    -> Result<(),CluEError>
  {
    let Ok(file) = File::create(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };

    let n_char_f64 = 16;
    let n_chars_per_line = 4*n_char_f64;
    let n_spins = self.bath_particles.len();

    let bytes_per_char = 32;

    let n_bytes = bytes_per_char*(n_spins+1)*n_chars_per_line;

    let mut stream = BufWriter::with_capacity(n_bytes,file);
    
    let line = "index,particle,x,y,z,active,cell_id,group\n";
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }    

    for (ii, particle) in self.bath_particles.iter().enumerate(){

      let r = &particle.coordinates - &self.pdb_origin; 

      let idx = self.get_reference_index_from_bath_index(ii)?;

      let cell_id = self.cell_id(ii)?;

      let particle_index_0 = self.primary_cell_indices[ii];

      let group_label = if cell_id == 0{
        match self.particle_config_ids[particle_index_0]{
          Some(cfg_idx) => &config.particles[cfg_idx].label,
          None => "none",
        }
      }else{
        match self.extracell_particle_config_ids[particle_index_0]{
          Some(cfg_idx) => &config.extracell_particles[cfg_idx].label,
          None => "none",
        }
      };

      let line = format!("{},{},{},{},{},{},{},{}\n",
          idx, particle.isotope.to_string(),
          r.x()/ANGSTROM, r.y()/ANGSTROM, r.z()/ANGSTROM,
          particle.active,cell_id,group_label);
    
      if stream.write(line.as_bytes()).is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }    
    }

    Ok(())
  }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  use crate::config::lexer::get_tokens_from_line;
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_structure(){
    let token_stream = get_tokens_from_line("
        input_structure_file = \"assets/TEMPO.pdb\";
        radius = 20; // angstroms.
        detected_spin_position = centroid_over_serials([28,29]);
        number_timepoints = [101];
        time_increments = [1e-7];
        cluster_method = cce;
        max_cluster_size = 2;
        magnetic_field = 1.2;

        #[filter(label = tempo)]
          elements in [H];

        #[spin_properties(label = tempo, isotope = 1H)]
          tunnel_splitting = 80e3; // Hz.
        ").unwrap();

    let mut config = Config::new();

    config.parse_token_stream(token_stream).unwrap();

    config.set_defaults().unwrap();

    let mut rng = ChaCha20Rng::from_entropy();

    let structure = Structure::build_structure(&mut rng,&config).unwrap();

    assert_eq!(structure.number_active(),19);
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
