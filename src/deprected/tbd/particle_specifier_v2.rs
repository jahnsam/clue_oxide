use super::physical_constants::*;

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

pub struct SpecifiedParticle{
  pub serial: Option<u32>,
  pub residues: Option<String> ,
  pub elements: Option<Element> ,
  pub residue_sequence_number: Option<u32>,
  pub distance: Option<f64>,
}
impl SpecifiedParticle{
  pub fn new() -> Self{
      SpecifiedParticle{
        serial: None,
        residue: None,
        element: None,
        residue_sequence_number: None,
        distance: None,
      }
  }
  //----------------------------------------------------------------------------
  pub fn specify(pdb_idx: usize, pdb: &PDB) -> Self{
      SpecifiedParticle{
        serial: Some(pdb.serial(pdb_idx)),
        residue: Some(pdb.residue(pdb_idx)) ,
        element: Some(pdb.element(pdb_idx)),
        residue_sequence_number: Some(pdb.residue_sequence_number(pdb_idx)) ,
        distance: Some(pdb.pos(pdb_idx).magnitude()),
      }
  }
}

impl PartialEq<ParticleSpecifier> for ParticleSpecifier {
    pub fn new() -> Self{
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

    //--------------------------------------------------------------------------
    fn eq(&self, other: &ParticleSpecifier) -> bool {

      if self.serials.len() != other.serials.len()
        || self.not_serials.len() != other.not_serials.len()
        || self.residues.len() != other.residues.len()
        || self.not_residues.len() != other.not_residues.len()
        || self.elementss.len() != other.elementss.len()
        || self.residue_sequence_numbers.len() 
          != other.residue_sequence_numbers.len()
        || self.not_residue_sequence_numbers.len() 
          != other.not_residue_sequence_numbers.len()
      {
        return false;
      }

      // Serials
      if !are_vecs_eqial(&self.serials, &other.serials){
        return false;
      }
      /*
      for ii in 0..self.serials.len() {
        if self.serials[ii] != other.serials[ii]{
          return false;
        }
      }
      for ii in 0..self.not_serials.len() {
        if self.not_serials[ii] != other.not_serials[ii]{
          return false;
        }
      }
      */

      // Residues
      if !are_vecs_eqial(&self.residues, &other.residues){
        return false;
      }

      /*
      for ii in 0..self.residues.len() {
        if self.residues[ii] != other.residues[ii]{
          return false;
        }
      }
      for ii in 0..self.not_residues.len() {
        if self.not_residues[ii] != other.not_residues[ii]{
          return false;
        }
      }
      */
      

      // Elements
      if !are_vecs_eqial(&self.elements, &other.elements){
        return false;
      }

      /*
      for ii in 0..self.elements.len() {
        if self.elements[ii] != other.elements[ii]{
          return false;
        }
      }
      for ii in 0..self.elements.len() {
        if self.not_elements[ii] != other.not_elements[ii]{
          return false;
        }
      }
      */

      // Residue Sequence Numbers
      if !are_vecs_eqial(&self.residue_sequence_numbers, 
          &other.residue_sequence_numbers){
        return false;
      }

      /*
      for ii in 0..self.residue_sequence_numbers.len() {
        if self.residue_sequence_numbers[ii] 
          != other.residue_sequence_numbers[ii]{
          return false;
        }
      }
      for ii in 0..self.not_residue_sequence_numbers.len() {
        if self.not_residue_sequence_numbers[ii] 
          != other.not_residue_sequence_numbers[ii]{
          return false;
        }
      }
      */

      // Distances
      if self.max_distance != other.max_distance{ return false;}
      if self.min_distance != other.min_distance{ return false;}

      true
    }
}

impl ParticleSpecifier{

  pub fn covers(&self, target: &SpecifiedParticle) -> bool {
    
    let mut does_serials_cover = false;

    // Serials
    if let Some(target_value) = target.serial{

      if !does_cover(target_value,&self.serials)
        || is_excluded(target_value, &self.not_serials){
          return false;
        }
      /*
      let mut does_cover = self.serials==0;
      for filter_value in &self.serials {
          does_cover |= *filter_value == target_value;
          if does_cover { break; }
      }
      if !does_cover{ return false; }

      for exclude_value in &self.not_serials{
        if target_value == *exclude_value { return false;}
      }
      */
    }
    
    // Residues
    if let Some(target_value) = target.residue{
      let mut does_cover = self.residues==0;
      for filter_value in &self.residues {
          does_cover |= *filter_value == target_value;
          if does_cover { break; }
      }
      if !does_cover{ return false; }

      for exclude_value in &self.not_residues{
        if target_value == *exclude_value { return false;}
      }
    }
    
    // Elements
    if let Some(target_value) = target.element{
      let mut does_cover = self.elements==0;
      for filter_value in &self.elements {
          does_cover |= *filter_value == target_value;
          if does_cover { break; }
      }
      if !does_cover{ return false; }

      for exclude_value in &self.not_elements{
        if target_value == *exclude_value { return false;}
      }
    }
    
        
    // Residue Sequence Numbers
    if let Some(target_value) = target.residue_sequence_number{
      let mut does_cover = self.residue_sequence_numbers==0;
      for filter_value in &self.residue_sequence_numbers {
          does_cover |= *filter_value == target_value;
          if does_cover { break; }
      }
      if !does_cover{ return false; }

      for exclude_value in &self.not_residue_sequence_numbers{
        if target_value == *exclude_value { return false;}
      }
    }
    

      // Distances
      if let Some(r) = target.distance{
        if r < self.max_distance { return false;}
        if r > self.min_distance { return false;}
      }

    true
  }
}

fn are_vecs_equal<T>(vec0: &Vec::<T>, vec1: &Vec::<T>)-> bool{
  if vec0.len() != vec1.len() {return false;}

  for ii in 0..vec0.len(){
    if vec0[ii] != vec1[ii]{ return false;}
  }
}

fn does_cover<T>(target_value: T, filter_values: Vec::<T> ) -> bool {
  let mut doescover = filter_values.len()==0;
  for filter_value in &filter_values {
    doescover |= *filter_value == target_value;
     if doescover { break; }
  }
  doescover
}

fn is_excluded<T>(target_value: T, exclude_values: &Vec::<T>) -> bool
{
  for exclude_value in &excluded_values{
    if target_value == *exclude_value { return true;}
  }

  false
}

