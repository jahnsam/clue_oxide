use crate::physical_constants::*;
use crate::structure::pdb::PDB;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,Default)]
pub struct SpecifiedParticle{
  pub serial: Option<u32>,
  pub residue: Option<String> ,
  pub element: Option<Element> ,
  pub residue_sequence_number: Option<u32>,
  pub distance: Option<f64>,
}
//------------------------------------------------------------------------------
impl SpecifiedParticle{
  pub fn new() -> Self {Default::default()}
  //----------------------------------------------------------------------------
  pub fn specify(pdb_idx: usize, pdb: &PDB) -> Self{
      SpecifiedParticle{
        serial: Some(pdb.serial(pdb_idx)),
        residue: Some(pdb.residue(pdb_idx)) ,
        element: Some(pdb.element(pdb_idx)),
        residue_sequence_number: Some(pdb.residue_sequence_number(pdb_idx)) ,
        distance: Some(pdb.coordinates(pdb_idx).magnitude()),
      }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct ParticleSpecifier{
  pub serials: Vec::<u32>,
  pub not_serials: Vec::<u32>,
  
  pub residues: Vec::<String> ,
  pub not_residues: Vec::<String> ,
  
  pub elements: Vec::<Element> ,
  pub not_elements: Vec::<Element> ,
  
  pub residue_sequence_numbers: Vec::<u32> ,
  pub not_residue_sequence_numbers: Vec::<u32>,
  
  pub max_distance: f64,
  pub min_distance: f64,
}
//------------------------------------------------------------------------------

impl PartialEq<ParticleSpecifier> for ParticleSpecifier {
    fn eq(&self, other: &ParticleSpecifier) -> bool {

      // Serials
      if !are_vecs_equal(&self.serials, &other.serials){
        return false;
      }

      // Residues
      if !are_vecs_equal(&self.residues, &other.residues){
        return false;
      }


      // Elements
      if !are_vecs_equal(&self.elements, &other.elements){
        return false;
      }


      // Residue Sequence Numbers
      if !are_vecs_equal(&self.residue_sequence_numbers, 
          &other.residue_sequence_numbers){
        return false;
      }


      // Distances
      if self.max_distance != other.max_distance{ return false;}
      if self.min_distance != other.min_distance{ return false;}

      true
    }
}

//------------------------------------------------------------------------------

impl Default for ParticleSpecifier{

    fn default() -> Self{
      ParticleSpecifier{
        serials: Vec::<u32>::new(),
        not_serials: Vec::<u32>::new(),
  
        residues: Vec::<String>::new() ,
        not_residues: Vec::<String>::new() ,
  
        elements: Vec::<Element>::new() ,
        not_elements: Vec::<Element>::new() ,
  
        residue_sequence_numbers: Vec::<u32>::new() ,
        not_residue_sequence_numbers: Vec::<u32>::new(),
  
        max_distance: f64::INFINITY,
        min_distance: 0.0,
       
      }
    }
}

//------------------------------------------------------------------------------
impl ParticleSpecifier{

  pub fn new() -> Self {Default::default()}
  //----------------------------------------------------------------------------
  pub fn covers(&self, target: &SpecifiedParticle) -> bool {
    
    // Serials
    if let Some(target_value) = target.serial{

      if !does_cover(target_value,&self.serials)
        || is_excluded(target_value, &self.not_serials){
          return false;
        }
    }
    
    // Residues
    if let Some(target_value) = target.residue.clone(){
      if !does_cover(target_value.clone(),&self.residues)
        || is_excluded(target_value, &self.not_residues){
          return false;
        }
    }
    
    // Elements
    if let Some(target_value) = target.element{
      if !does_cover(target_value,&self.elements)
        || is_excluded(target_value, &self.not_elements){
          return false;
        }
    }
    
        
    // Residue Sequence Numbers
    if let Some(target_value) = target.residue_sequence_number{
      if !does_cover(target_value,&self.residue_sequence_numbers)
        || is_excluded(target_value, &self.not_residue_sequence_numbers){
          return false;
        }
    }
    

      // Distances
      if let Some(r) = target.distance{
        if r > self.max_distance { return false;}
        if r < self.min_distance { return false;}
      }

    true
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn are_vecs_equal<T: std::cmp::PartialEq>
(vec0: &[T], vec1: &[T])-> bool{
  if vec0.len() != vec1.len() {return false;}

  for ii in 0..vec0.len(){
    if vec0[ii] != vec1[ii]{ return false;}
  }
  true
}
//------------------------------------------------------------------------------
fn does_cover<T: std::cmp::PartialEq>
(target_value: T, filter_values: &[T] ) -> bool {
  let mut doescover = filter_values.is_empty();
  for filter_value in filter_values {
    doescover |= *filter_value == target_value;
     if doescover { break; }
  }
  doescover
}
//------------------------------------------------------------------------------
fn is_excluded<T: std::cmp::PartialEq>
(target_value: T, excluded_values: &[T]) -> bool
{
  for exclude_value in excluded_values{
    if target_value == *exclude_value { return true;}
  }

  false
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_covers(){
    let mut specifier = ParticleSpecifier::new();
    let mut target = SpecifiedParticle::new();

    // Elements
    target.element = Some(Element::Hydrogen);
    assert!(specifier.covers(&target));

    specifier.elements.push(Element::Nitrogen);
    assert!(!specifier.covers(&target));

    specifier.elements.push(Element::Hydrogen);
    assert!(specifier.covers(&target));

    // Distances
    target.distance = Some(1.0);
    assert!(specifier.covers(&target));

    specifier.min_distance = 2.0;
    assert!(!specifier.covers(&target));

    target.distance = Some(3.0);
    assert!(specifier.covers(&target));

    specifier.max_distance = 2.5;
    assert!(!specifier.covers(&target));

    target.distance = Some(2.4);
    assert!(specifier.covers(&target));

    // Serials
    specifier.serials = vec![1,2,4,8];
    assert!(specifier.covers(&target));

    target.serial = Some(16);
    assert!(!specifier.covers(&target));

    specifier.serials.push(16);
    assert!(specifier.covers(&target));

    // Residues
    target.residue = Some(String::from("SOL"));
    assert!(specifier.covers(&target));
    
    specifier.residues.push(String::from("WAT"));
    assert!(!specifier.covers(&target));

    specifier.residues.push(String::from("SOL"));
    assert!(specifier.covers(&target));

    // Residue Sequence Numbers
    specifier.residue_sequence_numbers = vec![1,3,5];
    assert!(specifier.covers(&target));

    target.residue_sequence_number = Some(2);
    assert!(!specifier.covers(&target));

    specifier.residue_sequence_numbers.push(7);
    assert!(!specifier.covers(&target));

    specifier.residue_sequence_numbers.push(5);
    assert!(!specifier.covers(&target));

    specifier.residue_sequence_numbers.push(2);
    assert!(specifier.covers(&target));

  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
