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
  pub fn new(pdb_idx: usize, pdb: &PDB) -> Self{
      SpecifiedParticle{
        serial: Some(pdb.serial(pdb_idx)),
        residue: Some(pdb.residue(pdb_idx)) ,
        element: Some(pdb.element(pdb_idx)),
        residue_sequence_number: Some(pdb.residue_sequence_number(pdb_idx)) ,
        distance: pdb.pos(pdb_idx).magnitude(),
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

    pub fn specify(pdb_idx: usize, pdb: &PDB) -> Self{

      let r = pdb.pos(pdb_idx).magnitude();

      ParticleSpecifier{
        serials: Vec::<u32>::from([pdb.serial(pdb_idx)]),
        not_serials: Vec::<u32>::from([pdb.serial(pdb_idx)]),
  
        residues: Vec::<String>::from(pdb.residue(pdb_idx)); ,
        not_residues: Vec::<String>::from(pdb.residue(pdb_idx)); ,
  
        elements: Vec::<Element>::from([pdb.element(pdb_idx)]) ,
        not_elements: Vec::<Element>::from([pdb.element(pdb_idx)]) ,
  
        residue_sequence_numbers: Vec::<u32>::from(
            [pdb.residue_sequence_number(pdb_idx)]) ,
        not_residue_sequence_numbers: Vec::<u32>::from([
            [pdb.residue_sequence_number(pdb_idx)]) ,
  
        max_distance: r,
        min_distance: r,
       
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
      
      // Residues
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
      

      // Elements
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

      // Residue Sequence Numbers
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

      // Distances
      if self.max_distance != other.max_distance{ return false;}
      if self.min_distance != other.min_distance{ return false;}

      true
    }
}

impl ParticleSpecifier{

  pub fn covers(&self, other: &ParticleSpecifier) -> bool {
    
    let mut does_cover = true;
    let mut does_serials_cover = false;

    // Serials
    does_cover &= other.serials.len() == 0;
    'outer: for value0 in &self.serials {
      for value1 in &other.serials {
        does_cover |= *value0 == *value;
        if does_cover { break 'outer; }
      }
    } 
    if !does_cover{ return false; }

    for value0 in &self.not_serials {
      for value1 in &other.not_serials {
         if *value0 == *value{ return false; }
      }
    }
    
    // Residues
    does_residues_cover &= other.residues.len() == 0;
    'outer: for value0 in &self.residues {
      for value1 in &other.residues {
        does_cover |= *value0 == *value;
        if does_cover { break 'outer; }
      }
    }
    if !does_cover{ return false; }

    for value0 in &self.not_residues {
      for value1 in &other.not_residues {
         if *value0 == *value{ return false; }
      }
    }
    
    // Elements
    does_cover &= other.elements.len() == 0;
    'outer: for value0 in &self.elements {
      for value1 in &other.elements {
        does_cover |= *value0 == *value;
        if does_cover { break 'outer; }
      }
    }
    if !does_cover{ return false; }

    for value0 in &self.not_elements {
      for value1 in &other.not_elements {
         if *value0 == *value{ return false; }
      }
    }
        
    // Residue Sequence Numbers
    does_cover &= other.residue_sequence_numbers.len() == 0;
    'outer: for value0 in &self.residue_sequence_numbers {
      for value1 in &other.residue_sequence_numbers {
        does_cover |= *value0 == *value;
        if does_cover { break 'outer; }
      }
    }
    if !does_cover{ return false; }

    for value0 in &self.not_residue_sequence_numbers {
      for value1 in &other.not_residue_sequence_numbers {
         if *value0 == *value{ return false; }
      }
    }

      // Distances
      if self.max_distance > other.max_distance{ return false;}
      if self.min_distance < other.min_distance{ return false;}


    does_cover
  }
}
