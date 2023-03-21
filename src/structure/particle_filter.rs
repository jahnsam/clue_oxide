use crate::structure::Structure;
use crate::space_3d::Vector3D;
use crate::structure::Particle;
use crate::physical_constants::{Element,Isotope};
use crate::clue_errors::CluEError;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Filters are used to specify a set of bath particles.
// TODO: ensure filters can specify:
//  methyl groups to change their tunnel splitting,
//  hyperfine/quadrupole axes and values
//  particle isotope abundances
//  particles to co-exchange for deuteration
//  particle pbc behavior
#[derive(Debug,Clone,Default,PartialEq)]
pub struct ParticleFilter{

  pub indices: Vec::<usize>,
  pub not_indices: Vec::<usize>,
  
  //pub cell_ids: Vec::<usize>,
  //pub not_cell_ids: Vec::<usize>,

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

  pub bonded_indices: Vec::<usize>,
  pub not_bonded_indices: Vec::<usize>,

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
  //----------------------------------------------------------------------------
  pub fn filter_indices(&self, structure: &Structure,indices: &Vec::<usize>) 
    -> Vec::<usize>
  {
    indices.iter().filter_map( |idx| {
      let particle = &structure.bath_particles[*idx];
      self.does_pass_filter(structure, *idx,particle)
    }).collect()
  }
  //----------------------------------------------------------------------------
  pub fn filter(&self, structure: &Structure) -> Vec::<usize>{
    
    let particles = &structure.bath_particles;

    particles.iter().enumerate()
      .filter_map(|(idx,particle)| 
          self.does_pass_filter(structure, idx,particle)
    ).collect()
  }
  //----------------------------------------------------------------------------
  fn does_pass_filter(&self, structure: &Structure,
      idx: usize, particle: &Particle) -> Option<usize>
  {
    let particles = &structure.bath_particles;
      // Index
      if (!self.indices.is_empty() && !self.indices.contains(&idx))
       || self.not_indices.contains(&idx){
        return None;
      }

      /*
      // Cell ID
      if !self.cell_ids.is_empty(){
        if let Ok(cell_id) = structure.cell_id(idx){ 
          if !self.cell_ids.contains(&idx) || self.not_cell_ids.contains(&idx){
            return None;
        }}
      }
      */
      // Element
      if (!self.elements.is_empty() 
          && !self.elements.contains(&particle.element))
       || self.not_elements.contains(&particle.element){
       return None;
      }

      // Serial
      if let Some(serial) = particle.serial{ 
        if (!self.serials.is_empty() && !self.serials.contains(&serial))
         || self.not_serials.contains(&serial){
          return None;
      }}

      // Residue
      if let Some(res) = &particle.residue{ 
        if (!self.residues.is_empty() && !self.residues.contains(&res))
         || self.not_residues.contains(&res){
          return None;
      }}

      // Residue Sequence ID
      if let Some(res_seq) = particle.residue_sequence_number{ 
        if (!self.residue_sequence_numbers.is_empty() 
            && !self.residue_sequence_numbers.contains(&res_seq))
         || self.not_residue_sequence_numbers.contains(&res_seq){
          return None;
      }}

      // Isotope
      let isotope = &particle.isotope; 
      if (!self.isotopes.is_empty() && !self.isotopes.contains(&isotope))
       || self.not_isotopes.contains(&isotope){
        return None;
      }


      // Bonded To
      let mut pass = self.bonded_indices.is_empty();
      for neighbor in self.bonded_indices.iter(){
        if *neighbor == idx { continue; }
        if structure.connections.are_connected(idx,*neighbor){
          pass = true;
          break;
        }
      }
      if !pass{  return None;}

      // Not Bonded To
      for neighbor in self.not_bonded_indices.iter(){
        if *neighbor == idx { continue; }
        if structure.connections.are_connected(idx,*neighbor){
          return None
        }
      }

      let idx0 = structure.primary_cell_indices[idx];

      // Bonded Elements
      let mut pass = self.bonded_elements.is_empty();
      for bonded_element in self.bonded_elements.iter(){ 
        let Some(neighbors) = structure.connections.get_neighbors(idx0) else{
          break;
        };
        for neighbor_idx in neighbors{
          if particles[*neighbor_idx].element == *bonded_element{
            pass = true;
            break;
          }
        } 
      }
      if !pass{  return None;}

      
      // Not Bonded Elements
      for not_bonded_element in self.not_bonded_elements.iter(){ 
        let Some(neighbors) = structure.connections.get_neighbors(idx0) else{
          break;
        };
        for neighbor_idx in neighbors{
          if particles[*neighbor_idx].element == *not_bonded_element{
            return None;
          }
        } 
      }

      

      Some(idx)
    }
//------------------------------------------------------------------------------
  /// This function augments a particle filter, making it specific to the
  /// specified particle.  
  /// The secondary filter specifies how the filter should be restricted.
  pub fn augment_filter(
    &mut self,
    index_ref_particle: usize,
    secondary_filter: &SecondaryParticleFilter,
    structure: &Structure) -> Result<(),CluEError>
  {


    match secondary_filter{
      SecondaryParticleFilter::Bonded 
        => self.bonded_indices.push(index_ref_particle),
      
      SecondaryParticleFilter::SameMolecule 
        => {
        let id = structure.molecule_ids[index_ref_particle];
        for idx in structure.molecules[id].iter(){
          self.indices.push(*idx);
        }
      },
      
      /*
      SecondaryParticleFilter::SameResidueSequenceNumber 
        => {
        if let Some(res_seq_id) = structure.bath_particles[index_ref_particle]
          .residue_sequence_number{
        self.residue_sequence_numbers.push(res_seq_id);
          }else{
            return Err(CluEError::CannontAugmentFilter(index_ref_particle,
                secondary_filter.to_string()));
          }
      },
      */
    }


    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// Secondary particle filters allow particles filters to be modified so as to
/// allow the filter to filter based off the relationship of particles to a
/// specified particle.
#[derive(Debug,Clone,PartialEq)]
pub enum SecondaryParticleFilter{
  Bonded, // atoms bonded to particle
  //Particle, // the particle itself
  SameMolecule, // atoms on the same molecule as particle 
  //SameResidueSequenceNumber, // atoms with the same ResSeq as particle.
}
impl ToString for SecondaryParticleFilter{
  fn to_string(&self) -> String{
    match self{
      SecondaryParticleFilter::Bonded => String::from("Bonded"),
      SecondaryParticleFilter::SameMolecule
        => String::from("SameMolecule"),
    }
  } 
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum VectorSpecifier{
  Diff(SecondaryParticleFilter,SecondaryParticleFilter),
  Vector(Vector3D)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb;
  use crate::Config;
  use crate::config::DetectedSpinCoordinates;

  //----------------------------------------------------------------------------
  #[test]
  fn test_augment_filter(){
    let filename = "./assets/a_TEMPO_a_water_a_glycerol.pdb";
    let mut structures = pdb::parse_pdb(&filename).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::MeanOverSerials(vec![28,29]) );
    let structure = &mut structures[0];
    config.set_defaults();
    structure.build_primary_structure(&config).unwrap();


    let mut filter = ParticleFilter::new();
    filter.augment_filter(43,
        &SecondaryParticleFilter::Bonded,structure)
      .unwrap();
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![44,45]);

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    filter.augment_filter(28,
        &SecondaryParticleFilter::SameMolecule,structure)
      .unwrap();
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![2,3,4,6,7,8,10,11,13,14,16,17,20,21,22,24,25,26]);

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    filter.augment_filter(28,
        &SecondaryParticleFilter::SameMolecule,structure)
      .unwrap();
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![2,3,4,6,7,8,10,11,13,14,16,17,20,21,22,24,25,26]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_filter(){
    let filename = "./assets/a_TEMPO_a_water_a_glycerol.pdb";
    let mut structures = pdb::parse_pdb(&filename).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::MeanOverSerials(vec![28,29]) );
    let structure = &mut structures[0];
    config.set_defaults();
    structure.build_primary_structure(&config).unwrap();

  
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Oxygen];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![28,32,36,41,43]);

    let mut filter = ParticleFilter::new();
    filter.not_elements = vec![Element::Hydrogen,Element::Carbon];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![27,28,32,36,41,43]);

    let mut filter = ParticleFilter::new();
    filter.serials = vec![28,29];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![27,28]);

    let mut filter = ParticleFilter::new();
    filter.bonded_elements = vec![Element::Oxygen];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![27,29,33,34,37,38,42,44,45]);

    let mut filter = ParticleFilter::new();
    filter.not_bonded_elements = vec![Element::Carbon];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![28,33,37,42,43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.residues = vec!["SOL".to_string()];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.not_residues = vec!["TEM".to_string(),"MGL".to_string()];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.residue_sequence_numbers = vec![1502];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.not_residue_sequence_numbers = vec![1,2];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    filter.bonded_elements = vec![Element::Oxygen];
    let indices = filter.filter(structure);
    assert_eq!(indices,vec![33,37,42,44,45]);

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    filter.not_bonded_elements = vec![Element::Oxygen];
    let indices = filter.filter(structure);
    assert_eq!(indices.len(),23);
    assert_eq!(indices,vec![2,3,4,6,7,8,10,11,13,14,16,17,20,21,22,24,25,26,
    30,31,35,39,40]);
  }
}
