pub mod pdb;
pub mod extended_structure;
pub mod primary_structure;
pub mod exchange_groups;
pub mod particle;
pub mod particle_filter;

//use crate::config::particle_config;
use crate::config::{Config, particle_config::ParticleProperties};
use crate::config::particle_config::IsotopeProperties;
use crate::config::particle_config::TensorSpecifier;
use crate::clue_errors::CluEError;
use crate::integration_grid::IntegrationGrid;
use crate::structure::particle::Particle;
use crate::structure::particle_filter::*;
use crate::space_3d::Vector3D;
use crate::cluster::adjacency::AdjacencyList;
use crate::structure::exchange_groups::ExchangeGroupManager;
use crate::physical_constants::Isotope;


use rand_chacha::ChaCha20Rng;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct DetectedSpin{
  pub isotope: Isotope,
  pub weighted_coordinates: IntegrationGrid,
  pub transition: [usize;2],
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

  // origin of the pdb frame
  pdb_origin: Vector3D,

  // list indices indicating the particle that each particle is a periodic
  // boundary condition copy of  
  primary_cell_indices: Vec::<usize>,
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
      molecules: Vec::<Vec::<usize>>::new(),
      molecule_ids: Vec::<usize>::new(),
      cosubstitution_groups: Vec::<Vec::<usize>>::new(),
      exchange_groups: None,
      particle_config_ids: Vec::<Option<usize>>::new(),
      pdb_origin: Vector3D::zeros(),
      primary_cell_indices: Vec::<usize>::new(),
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
    let indices = particle_filter.filter(&self);
    let mut out = Vec::<&Particle>::with_capacity(indices.len());

    for idx in indices{
      out.push(&self.bath_particles[idx]);
    }

    out
  }
  //----------------------------------------------------------------------------
  pub fn cell_id(&self, idx: usize) -> Result<usize,CluEError>
  {
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
   
    let mol_id = self.molecule_ids[idx0];// + cell_id*self.molecules.len();

    mol_id
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
    let particle_index = self.primary_cell_indices[particle_index];

    let Some(p_cfg_id) = self.particle_config_ids[particle_index] else{
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
  // Since multiple filters may cover the same particle, 
  // the last filter defined that covers the particle is the one paired with
  // said particle.
  fn pair_particle_configs(&mut self, config: &Config){
  
    // Initialize pairings.
    let n_particles = self.bath_particles.len();

    self.particle_config_ids = Vec::<Option<usize>>::with_capacity(n_particles);
    for _ii in 0..n_particles{
      self.particle_config_ids.push(None);
    }

    // Get particle configs.
    let particle_configs = &config.particles;

    // Pair particles with configs.
    for (id,particle_config) in particle_configs.iter().enumerate(){
      if let Some(filter) = &particle_config.filter{
        let indices = filter.filter(self);

        for idx in indices{
          self.particle_config_ids[idx] = Some(id);
        }
      
      } 
    }
    
  }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



