use crate::space_3d::Vector3D;
use crate::cluster::adjacency::AdjacencyList;
use crate::structure::exchange_groups::ExchangeGroup;
use crate::physical_constants::{Element,Isotope};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct Particle{
  element: Element,
  coordinates: Vector3D,

  serial: Option<u32>,
  residue: Option<String>,
  residue_sequence_number: Option<u32>,
  exchange_group: Option<ExchangeGroup>,
  exchange_group_id: Option<usize>,  

  isotope: Option<Isotope>,

  cell_id: usize,
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct DetectedSpin{
  isotope: Isotope,
  coordinates: Vec::<(Vector3D,f64)>,
  transition: [usize;2],
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct Structure{
  particles: Vec::<Particle>,
  connections: AdjacencyList,
  cell_offsets: Vec::<Vector3D>,
}

impl Structure{
  pub fn find<'a>(&'a self, particle_filter: &ParticleFilter) 
    -> Vec::<&'a Particle>
  {
    let indices = particle_filter.filter(&self);
    let mut out = Vec::<&Particle>::with_capacity(indices.len());

    for idx in indices{
      out.push(&self.particles[idx]);
    }

    out
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct ParticleFilter{
  elements: Vec::<Element>, 
  not_elements: Vec::<Element>, 

  serials: Vec::<u32>, 
  not_serials: Vec::<u32>, 

  residues: Vec::<String>, 
  not_residues: Vec::<String>, 

  residue_sequence_numbers: Vec::<u32>,
  not_residue_sequence_numbers: Vec::<u32>,

  //exchange_groups: Vec::<ExchangeGroup>,
  //not_exchange_groups: Vec::<ExchangeGroup>,

  isotopes: Vec::<Isotope>,
  not_isotopes: Vec::<Isotope>,

  bonded_to: Vec::<usize>,
  not_bonded_to: Vec::<usize>,

  //within_distance_of: Vec::<usize>,
  //not_within_distance_of: Vec::<usize>,

  //distance_away_from: Vec::<usize>,
  //not_distance_away_from: Vec::<usize>,

}

impl ParticleFilter{
  pub fn filter(&self, structure: &Structure) -> Vec::<usize>{
    
    let particles = &structure.particles;

    particles.iter().enumerate().filter_map(|(idx,particle)| {

      // Element
      if (self.elements.is_empty() 
          && !self.elements.contains(&particle.element))
       || self.not_elements.contains(&particle.element){
       return None;
      }

      // Serial
      if let Some(serial) = particle.serial{ 
        if (self.serials.is_empty() && !self.serials.contains(&serial))
         || self.not_serials.contains(&serial){
          return None;
      }}

      // Residue
      if let Some(res) = &particle.residue{ 
        if (self.residues.is_empty() && !self.residues.contains(&res))
         || self.not_residues.contains(&res){
          return None;
      }}

      // Residue Sequence ID
      if let Some(res_seq) = particle.residue_sequence_number{ 
        if (self.residue_sequence_numbers.is_empty() 
            && !self.residue_sequence_numbers.contains(&res_seq))
         || self.not_residue_sequence_numbers.contains(&res_seq){
          return None;
      }}

      // Isotope
      if let Some(isotope) = &particle.isotope{ 
        if (self.isotopes.is_empty() && !self.isotopes.contains(&isotope))
         || self.not_isotopes.contains(&isotope){
          return None;
      }}


      // Bonded To
      let mut pass = self.bonded_to.is_empty();
      for neighbor in self.bonded_to.iter(){
        if *neighbor == idx { continue; }
        if structure.connections.are_connected(idx,*neighbor){
          pass = true;
          break;
        }
      }
      if !pass{  return None;}

      // Not Bonded To
      for neighbor in self.not_bonded_to.iter(){
        if *neighbor == idx { continue; }
        if structure.connections.are_connected(idx,*neighbor){
          return None
        }
      }

      

      Some(idx)
    }).collect()
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
