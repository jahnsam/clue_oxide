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
/*
struct ExchangeGroupManager{
  exchange_groups: Vec::<ExchangeGroup>,
  exchange_group_ids: Vec::<Option<usize>>,
  exchange_coupling: Vec::<f64>, // one entry per exchange_group
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct Structure{
  pub detected_particle: Option<DetectedSpin>,
  pub bath_particles: Vec::<Particle>,
  bath_spins_indices: Vec::<usize>,
  cell_indices: Vec::<Vec::<Option<usize>>>,
  pub cell_offsets: Vec::<Vector3D>,
  pub connections: AdjacencyList,
  molecules: Vec::<Vec::<usize>>,
  molecule_ids: Vec::<usize>,
  //cosubstitute: Option<AdjacencyList>,
  cosubstitution_groups: Vec::< Vec::<usize> >,
  exchange_groups: Option<ExchangeGroupManager>,
  particle_config_ids: Vec::<Option<usize>>,
  pdb_origin: Vector3D,
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
      //cosubstitute: None,
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
  pub fn extract_isotope_properties<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a IsotopeProperties>
  {
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



