pub mod pdb;
pub mod extended_structure;
pub mod primary_structure;
pub mod exchange_groups;
pub mod particle;
pub mod particle_filter;

//use crate::config::particle_config;
use crate::config::{Config, 
  particle_config::{ParticleConfig,ParticleProperties}};
use crate::clue_errors::CluEError;
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
  //cosubstitute: Option<AdjacencyList>,
  cosubstitution_groups: Vec::< Vec::<usize> >,
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
      //cosubstitute: None,
      cosubstitution_groups: Vec::<Vec::<usize>>::new(),
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
  // TODO: Slow for large systems: find_cosubstitution_groups().
  // TODO: test find_cosubstitution_groups()
  // This function goes through each particle, applies secondary filters,
  // and records the cosubstitution_groups.
  fn find_cosubstitution_groups(&mut self, config: &Config) 
    -> Result<(),CluEError>
  {
  
    let mut cosubstitution_group_ids 
      = Vec::<Option<usize>>::with_capacity(self.bath_particles.len());
    for ii in 0..self.bath_particles.len(){
      cosubstitution_group_ids.push(None);
    }
    
    let mut current_cosub_id = 0;


    // Loop through bath particles.
    for (idx,particle) in self.bath_particles.iter().enumerate(){

      // Check if there are custom properties for this particle.
      let Some(id) = self.particle_config_ids[idx] else{continue;};
      let Some(properties) = &config.particles[id].properties else{continue;}; 
    
      // Check if this particle has a filter.
      let mut filter: ParticleFilter;
      if let Some(fltr) = &config.particles[id].filter{
        filter= fltr.clone();
      }else{
        filter = ParticleFilter::new();
      }

      // Set cosubstitution group for this particle.
      update_cosubstitution_ids(&mut cosubstitution_group_ids, 
        idx, current_cosub_id, filter.clone(),properties,self)?;
      current_cosub_id += 1;

    }

    // TODO: n_none should be findable in one line.
    let mut n_none: usize  = 0;
    for opt in cosubstitution_group_ids.iter(){
      if *opt == None{ n_none += 1;}
    }

    // cosubstitutions
    let n_groups = n_none + current_cosub_id;
    let mut cosubstitution_groups 
      = Vec::<Vec::<usize>>::with_capacity(n_groups);
    for ii in 0..current_cosub_id{
      cosubstitution_groups.push(Vec::<usize>::new());
    }

    /*
    // TODO: This loop is slow for large systems.
    for ii in 0..current_cosub_id{
      let vec_indices = cosubstitution_group_ids.iter().enumerate()
        .filter(|(_, &r)| r == Some(ii))
        .map(|(index, _)| index).collect();

      cosubstitution_groups.push(vec_indices);

    }
    */
    // TODO: this loop appears to be a faster version.  
    // Delete above loop after testing.
    for (ii,id_option) in cosubstitution_group_ids.iter().enumerate(){
      if let Some(id) = id_option{
        cosubstitution_groups[*id].push(ii);
      }else{
        cosubstitution_groups.push(vec![ii]);
      }
    }

    
    for (idx,opt) in cosubstitution_group_ids.iter().enumerate(){
      if *opt != None{ continue;}
      let vec_indices = vec![idx];
      cosubstitution_groups.push(vec_indices);
    }

    self.cosubstitution_groups = cosubstitution_groups;

    Ok(())
  }
}
//------------------------------------------------------------------------------
// This function finds all the particles that cosubstitute with particle idx,
// and records them in cosubstitution_group_ids.
fn update_cosubstitution_ids(
    cosubstitution_group_ids: &mut Vec::<Option<usize>>,
    idx: usize,
    current_cosub_id: usize,
    mut filter: ParticleFilter,
    properties: &ParticleProperties,
    structure: &Structure)
    -> Result<(),CluEError>
  {
      // Check if properties defines a cosubstitution set.
      if let Some(cosubstitute) = &properties.cosubstitute{

        // Augment filter with criteria specific to particle idx.
        filter.augment_filter(idx,&cosubstitute,structure)?;

        // Find all particle indices that fit the criteria.
        let indices = filter.filter(structure);

        // Loop through all particles found.
        for index in indices.iter(){
          // Ensure the particle is not part of a different group.
          if cosubstitution_group_ids[*index]==None{
            return Err(CluEError::MultipleCosubstitutionGroups(*index));
          }
          // Assign the particle to this cosubstitution group.
          cosubstitution_group_ids[*index] = Some(current_cosub_id)
        }
      }else{
        // Put the particle in a cosubstitution group by itself.
        cosubstitution_group_ids[idx] = Some(current_cosub_id)
      }
    Ok(())  
  }
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



