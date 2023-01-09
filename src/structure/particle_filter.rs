use crate::structure::Structure;
use crate::space_3d::Vector3D;
//use crate::structure::exchange_groups::ExchangeGroup;
use crate::physical_constants::{Element,Isotope};
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Filters are used to specify a set of bath particles.
// TODO: ensure filters can specify:
//  methyl groups to change their tunnel splitting,
//  hyperfine/quadrupole axes and values
//  particle isotope abundances
//  particles to co-exchange for deuteration
//  particle pbc behavior
#[derive(Debug,Clone,Default)]
pub struct ParticleFilter{
  pub elements: Vec::<Element>, 
  pub not_elements: Vec::<Element>, 

  pub serials: Vec::<u32>, 
  pub not_serials: Vec::<u32>, 

  pub residues: Vec::<String>, 
  pub not_residues: Vec::<String>, 

  pub residue_sequence_numbers: Vec::<u32>,
  pub not_residue_sequence_numbers: Vec::<u32>,

  //exchange_groups: Vec::<ExchangeGroup>,
  //not_exchange_groups: Vec::<ExchangeGroup>,

  pub isotopes: Vec::<Isotope>,
  pub not_isotopes: Vec::<Isotope>,

  pub bonded_to: Vec::<usize>,
  pub not_bonded_to: Vec::<usize>,

  //within_distance_of: Vec::<usize>,
  //not_within_distance_of: Vec::<usize>,

  //distance_away_from: Vec::<usize>,
  //not_distance_away_from: Vec::<usize>,

  pub bonded_elements: Vec::<Element>, 
  pub not_bonded_elements: Vec::<Element>, 

  pub bonded_serials: Vec::<u32>, 
  pub not_bonded_serials: Vec::<u32>, 

  pub bonded_residues: Vec::<String>, 
  pub not_bonded_residues: Vec::<String>, 

  pub bonded_residue_sequence_numbers: Vec::<u32>,
  pub not_bonded_residue_sequence_numbers: Vec::<u32>,


}

impl ParticleFilter{
  pub fn new() -> Self{
    Default::default()
  }

  pub fn filter(&self, structure: &Structure) -> Vec::<usize>{
    
    let particles = &structure.bath_particles;

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
      let isotope = &particle.isotope; 
      if (self.isotopes.is_empty() && !self.isotopes.contains(&isotope))
       || self.not_isotopes.contains(&isotope){
        return None;
      }


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

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// The String that most SecondaryParticleFilters hold is the label to a
// ParticleFilter.  Config has a HashMap<String,ParticleFilter> to access
// the filter.
//
// Each enum entry corresponds to a function that takes the intersection of
// particles found by the supplied filter and a secondary filter.
#[derive(Debug,Clone)]
pub enum SecondaryParticleFilter{
  Bonded(String), // atoms bonded to particle
  Particle, // the particle itself
  SameMolecule(String), // atoms on the same molecule as particle 
  SameResidueSequenceNumber(String), // atoms with the same ResSeq as particle.
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub enum VectorSpecifier{
  Diff(SecondaryParticleFilter,SecondaryParticleFilter),
  Vector(Vector3D)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
