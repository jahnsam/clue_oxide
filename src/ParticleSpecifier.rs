use super::Element::*;

pub struct ParticleSpecifier{
  pdb_ids: Vec::<u32>,
  residue: String,
  element: Option<Element>,
}

impl PartialEq<ParticleSpecifier> for ParticleSpecifier {
    fn eq(&self, other: &ParticleSpecifier) -> bool {

      if self.pdb_ids.len() != other.pdb_ids.len(){
        return false;
      }
      for ii in 0..self.pdb_ids.len() {
        if self.pdb_ids[ii] != other.pdb_ids[ii]{
          return false;
        }
      }
      
      if self.residue != other.residue {return false;}
      
      let same_elements = match (&self.element, &other.element) {
        (Some(a), Some(b)) => *a==*b,
        (None, None) => true,
        _ => false};

      if !same_elements{ return false };

      true
    }
}

impl ParticleSpecifier{

  pub fn covers(&self, other: &ParticleSpecifier) -> bool {
    
    let mut does_cover = false;

    'outer: for id0 in &self.pdb_ids {
      for id1 in &other.pdb_ids {

        does_cover |= id0==id1;

        if does_cover {
          break 'outer; 
        }
      }
    }

    if self.residue != "all" {
      does_cover &= self.residue == other.residue;
    }

   
    does_cover &= match (&self.element, &other.element) {
      (Some(a), Some(b)) => *a==*b,
        (None, _) => true,
      _ => false,};

    does_cover
  }
}
