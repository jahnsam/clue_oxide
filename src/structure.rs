pub mod pdb;
//pub mod primary_structure;
pub mod exchange_groups;
use crate::particle::Particle;
use crate::particle_filter::ParticleFilter;
use crate::space_3d::Vector3D;
use crate::cluster::adjacency::AdjacencyList;
use crate::structure::exchange_groups::ExchangeGroup;
use crate::physical_constants::{Element,Isotope};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct DetectedSpin{
  isotope: Isotope,
  //weighted_coordinates: IntegrationGrid,
  weighted_coordinates: Vec::<(Vector3D,f64)>,
  transition: [usize;2],
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct Structure{
  pub detected_particle: DetectedSpin,
  pub bath_particles: Vec::<Particle>,
  pub connections: AdjacencyList,
  pub coexchange: AdjacencyList,
  pub cell_offsets: Vec::<Vector3D>,
}

impl Structure{
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
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

