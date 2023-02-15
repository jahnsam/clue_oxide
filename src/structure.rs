pub mod pdb;
pub mod extended_structure;
pub mod primary_structure;
pub mod exchange_groups;
pub mod particle;
pub mod particle_filter;

//use crate::config::particle_config;
use crate::config::{Config, particle_config::ParticleConfig};
use crate::integration_grid::IntegrationGrid;
use crate::structure::particle::Particle;
use crate::structure::particle_filter::ParticleFilter;
use crate::space_3d::Vector3D;
use crate::cluster::adjacency::AdjacencyList;
use crate::structure::exchange_groups::ExchangeGroupManager;
use crate::physical_constants::Isotope;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct DetectedSpin{
  isotope: Isotope,
  weighted_coordinates: IntegrationGrid,
  transition: [usize;2],
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
  detected_particle: Option<DetectedSpin>,
  pub bath_particles: Vec::<Particle>,
  bath_spins_indices: Vec::<usize>,
  pub connections: AdjacencyList,
  pub cell_offsets: Vec::<Vector3D>,
  molecules: Option<AdjacencyList>,
  cosubstitute: Option<AdjacencyList>,
  exchange_groups: Option<ExchangeGroupManager>,
  particle_config_ids: Vec::<Option<usize>>,
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
      connections,
      cell_offsets,
      molecules: None,
      cosubstitute: None,
      exchange_groups: None,
      particle_config_ids: Vec::<Option<usize>>::new(),
    }

  }
  //----------------------------------------------------------------------------
  pub fn number(&self) -> usize{
    self.bath_particles.len()
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
  //----------------------------------------------------------------------------
  //TODO: finish apply_secondary_filters()
  // This function goes through each particle, applies secondary filters,
  // and records the results.
  fn apply_secondary_filters(&mut self, config: &Config){
  
    for (idx,particle) in self.bath_particles.iter().enumerate(){

      let Some(id) = self.particle_config_ids[idx] else{continue;};

      let Some(properties) = &config.particles[id].properties else{continue;}; 
    
      // cosubstitutions: update AdjacencyList

      // hyperfine tensors: Vec::<Option<Tensor3D>> or spare format?

      // electric quadrupole tensors: Vec::<Option<Tenso3D>> or spare format?

    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



