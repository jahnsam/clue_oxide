pub mod pdb;
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
use crate::structure::exchange_groups::ExchangeGroup;
use crate::physical_constants::Isotope;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct DetectedSpin{
  isotope: Isotope,
  weighted_coordinates: IntegrationGrid,
  transition: [usize;2],
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
struct ExchangeGroupManager{
  exchange_groups: Vec::<ExchangeGroup>,
  exchange_group_ids: Vec::<usize>,
  exchange_coupling: Vec::<f64>,
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
  fn pair_particle_configs(&mut self, config: &Config){
  
    // Initialize pairings.
    let n_particles = self.bath_particles.len();

    self.particle_config_ids = Vec::<Option<usize>>::with_capacity(n_particles);
    for _ii in 0..n_particles{
      self.particle_config_ids.push(None);
    }

    // Get particle configs.
    let particle_configs: &Vec::<ParticleConfig>;
    if let Some(pc) = &config.particles{
      particle_configs = pc;
    }else{
      return;
    }

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



