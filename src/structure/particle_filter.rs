
use crate::clue_errors::CluEError;
use crate::Config;
use crate::physical_constants::{Element,Isotope};
use crate::space_3d::Vector3D;
use crate::math;
use crate::structure::{Particle,Structure};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ParticleFilter` specify information to filter out a select subset of
/// particles from a `Structure`.
#[derive(Debug,Clone,Default,PartialEq)]
pub struct ParticleFilter{

  pub indices: Vec::<usize>,
  pub not_indices: Vec::<usize>,
  
  pub cell_ids: Vec::<usize>,
  pub not_cell_ids: Vec::<usize>,

  pub elements: Vec::<Element>, 
  pub not_elements: Vec::<Element>, 

  pub serials: Vec::<u32>, 
  pub not_serials: Vec::<u32>, 

  pub names: Vec::<String>, 
  pub not_names: Vec::<String>, 

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

  pub within_distance: Option<f64>,
  pub not_within_distance: Option<f64>,

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
  /// This function implements `new()` for `ParticleFilter`.
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------
  /// This function find the indices of the particles that pass the filter
  /// that also have indices within `indices`.
  pub fn filter_indices(&self, structure: &Structure,indices: &[usize]) 
    -> Vec::<usize>
  {
    indices.iter().filter_map( |idx| {
      let particle = &structure.bath_particles[*idx];
      self.does_pass_filter(structure, *idx,particle)
    }).collect()
  }
  //----------------------------------------------------------------------------
  /// This function find the indices of the particles that pass the filter.
  pub fn filter(&self, structure: &Structure) -> Vec::<usize>{
    
    let particles = &structure.bath_particles;

    particles.iter().enumerate()
      .filter_map(|(idx,particle)| 
          self.does_pass_filter(structure, idx,particle)
    ).collect()
  }
  //----------------------------------------------------------------------------
  // This function check if the indicated particle passes the filter.
  fn does_pass_filter(&self, structure: &Structure,
      idx: usize, particle: &Particle) -> Option<usize>
  {
    let particles = &structure.bath_particles;
      // Index
      if (!self.indices.is_empty() && !self.indices.contains(&idx))
       || self.not_indices.contains(&idx){
        return None;
      }

      
      // Cell ID
      if !self.cell_ids.is_empty() || !self.not_cell_ids.is_empty(){
        if let Ok(cell_id) = structure.cell_id(idx){ 
          if !self.cell_ids.contains(&cell_id) 
            || self.not_cell_ids.contains(&cell_id){
            return None;
        }}
      }
      
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

      // Atom Name
      if let Some(name) = &particle.name{ 
        if (!self.names.is_empty() && !self.names.contains(name))
         || self.not_names.contains(name){
          return None;
      }}

      // Residue
      if let Some(res) = &particle.residue{ 
        if (!self.residues.is_empty() && !self.residues.contains(res))
         || self.not_residues.contains(res){
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
      if (!self.isotopes.is_empty() && !self.isotopes.contains(isotope))
       || self.not_isotopes.contains(isotope){
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

      // Within Distance
      if let Some(r_max) = self.within_distance{
        let r = particle.coordinates.norm();
        if r > r_max{
          return None;
        }
      }

      // Not Within Distance
      if let Some(r_min) = self.not_within_distance{
        let r = particle.coordinates.norm();
        if r < r_min{
          return None;
        }
      }
      

      Some(idx)
    }
//------------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// Secondary particle filters allow particles filters to be modified so as to
/// allow the filter to filter based off the relationship of particles to a
/// specified particle.
#[derive(Debug,Clone,PartialEq)]
pub enum SecondaryParticleFilter{
  Bonded, // atoms bonded to particle
  Filter,
  Particle, // the particle itself
  SameMolecule, // atoms on the same molecule as particle 
  //SameResidueSequenceNumber, // atoms with the same ResSeq as particle.
}
impl ToString for SecondaryParticleFilter{
  fn to_string(&self) -> String{
    match self{
      SecondaryParticleFilter::Bonded => String::from("bonded"),
      SecondaryParticleFilter::Filter => String::from("filter"),
      SecondaryParticleFilter::Particle => String::from("particle"),
      SecondaryParticleFilter::SameMolecule
        => String::from("same_molecule"),
    }
  } 
}
impl SecondaryParticleFilter{
  pub fn from(secondary_filter: &str) 
    -> Result<SecondaryParticleFilter,CluEError>
  {
    match secondary_filter{
      "bonded" => Ok(SecondaryParticleFilter::Bonded),
      "filter" | "group" => Ok(SecondaryParticleFilter::Filter),
      "particle" => Ok(SecondaryParticleFilter::Particle),
      "same_molecule" => Ok(SecondaryParticleFilter::SameMolecule),
      _ => Err(CluEError::CannotParseSecondaryParticleFilter(
            secondary_filter.to_string())),
    }
  }
  //----------------------------------------------------------------------------
  pub fn filter(&self, particle_index_opt: Option<usize>, label: &str, 
      structure: &Structure, config: &Config ) -> Result<Vec::<usize>,CluEError>
  {

    let Some((_id, p_cfg)) = config.find_particle_config(label) else{
      return Err(CluEError::MissingFilter(label.to_string()));
    };

    let Some(filter) = &p_cfg.filter else{
      return Err(CluEError::MissingFilter(label.to_string()));
    };

    let get_particle_index = |sec_filt: &SecondaryParticleFilter|{
      let Some(particle_index) = particle_index_opt else{
        return
          Err(CluEError::SecondaryFilterRequiresAnIndex(sec_filt.to_string()));
      };
      Ok(particle_index)
    };
    //let mut indices = Vec::<usize>::new();
  
    let indices = match self{
      //
      SecondaryParticleFilter::Bonded => {
        let particle_index =get_particle_index(self)?;
        if let Some(bonded) = structure.connections
          .get_neighbors(particle_index){
          filter.filter_indices(structure,bonded)
        }else{
          return Err(CluEError::BondsAreNotDefined);
        }
      },
      SecondaryParticleFilter::Filter => filter.filter(structure),
      //
      SecondaryParticleFilter::Particle => {
        let particle_index =get_particle_index(self)?;
        vec![particle_index]
      },
      //
      SecondaryParticleFilter::SameMolecule => {
        let particle_index =get_particle_index(self)?;
        let mol_id = structure.molecule_ids[particle_index];
        filter.filter_indices(structure,&structure.molecules[mol_id])
      },
      //
    };
    Ok(indices)
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum VectorSpecifier{
  //CentroidOverSerials(Vec::<u32>),
  Diff(SecondaryParticleFilter,String,SecondaryParticleFilter,String),
  Vector(Vector3D)
}
//------------------------------------------------------------------------------
impl VectorSpecifier{
  pub fn to_vector3d(&self, particle_index_opt: Option<usize>, 
      structure: &Structure, config: &Config) -> Result<Vector3D,CluEError>
  {
    match self{ 
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      //VectorSpecifier::CentroidOverSerials(serials) 
      //  => structure.centroid_over_serials(serials.clone()),
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      VectorSpecifier::Diff(sec_fltr_0,label_0,sec_fltr_1,label_1) => {

        let mut indices_0 = sec_fltr_0.filter(particle_index_opt, label_0, 
            structure, config)?;

        let mut indices_1 = sec_fltr_1.filter(particle_index_opt, label_1, 
            structure, config)?;


        if indices_0.len() == 2 && indices_1.len() == 2{
          indices_0.append(&mut indices_1);
          indices_0 = math::unique(indices_0);
          indices_0.sort();
          
          if indices_0.len() != 2 {
            return Err(CluEError::VectorSpecifierDoesNotSpecifyUniqueVector(
                format!("...diff({}..., {}...",
                  sec_fltr_0.to_string(), sec_fltr_1.to_string(), )));
          }

          indices_1 = vec![indices_0[1]];
          indices_0 = vec![indices_0[0]];

        }

        if indices_0.len() != 1 {
          return Err(CluEError::VectorSpecifierDoesNotSpecifyUniqueVector(
                format!("...diff({}...",sec_fltr_0.to_string() )));
        }
        let idx0 = indices_0[0]; 

        if indices_1.len() != 1 {
          return Err(CluEError::VectorSpecifierDoesNotSpecifyUniqueVector(
                format!("...diff(..., {}...",sec_fltr_1.to_string() )));
        }
        let idx1 = indices_1[0]; 

        if idx0 == idx1{
          return Err(CluEError::VectorSpecifierDoesNotSpecifyUniqueVector(
                format!("...diff({}..., {}...",
                  sec_fltr_0.to_string(), sec_fltr_1.to_string(), )));
        }

        let vector3d = &structure.bath_particles[idx1].coordinates 
          - &structure.bath_particles[idx0].coordinates;
        
        Ok(vector3d)
      },
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      VectorSpecifier::Vector(vector3d) => Ok(vector3d.clone()),
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  use crate::config::particle_config::ParticleConfig;
  use crate::structure::pdb;
  use crate::Config;
  use crate::config::DetectedSpinCoordinates;

  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_VectorSpecifier_to_vector3d(){
    let filename = "./assets/TEMPO.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&config).unwrap();

    config.particles.push( ParticleConfig::new("nitrogen".to_string()) );
    let label = String::from("oxygen");
    config.particles.push( ParticleConfig::new(label.clone()) );

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Nitrogen];
    config.particles[0].filter = Some(filter);

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Oxygen];
    config.particles[1].filter = Some(filter);

    let vector_specifier = VectorSpecifier::Diff(
        SecondaryParticleFilter::Particle,"nitrogen".to_string(),
        SecondaryParticleFilter::Bonded,"oxygen".to_string());

    let r_n = Vector3D::from([36.440e-10, 36.900e-10,  37.100e-10]);
    let r_o = Vector3D::from([35.290e-10,  36.430e-10,  37.810e-10]);

    let delta_r = &r_o - &r_n;

    let particle_index = 27;
    assert_eq!(structure.bath_particles[particle_index].element, 
        Element::Nitrogen);
    let r = vector_specifier.to_vector3d(Some(particle_index),
        &structure,&config).unwrap();

    assert_eq!(r,delta_r);

  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_SecondaryParticleFilter_filter(){
  
    let filename = "./assets/a_TEMPO_a_water_a_glycerol.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&config).unwrap();

    let label = String::from("test_label");
    config.particles.push( ParticleConfig::new(label.clone()) );
    let mut filter = ParticleFilter::new();

    config.particles[0].filter = Some(filter.clone());

    let particle_index = 43;
    let secondary_filter = SecondaryParticleFilter::Bonded;
    let indices = secondary_filter.filter(Some(particle_index), &label, 
        &structure, &config ).unwrap(); 
    assert_eq!(indices,vec![44,45]);

    let secondary_filter = SecondaryParticleFilter::Particle;
    let indices = secondary_filter.filter(Some(particle_index), &label, 
        &structure, &config ).unwrap(); 
    assert_eq!(indices,vec![43]);


    let particle_index = 28;
    let secondary_filter = SecondaryParticleFilter::SameMolecule;
    filter.elements = vec![Element::Hydrogen];
    config.particles[0].filter = Some(filter);
    let indices = secondary_filter.filter(Some(particle_index), &label, 
        &structure, &config ).unwrap(); 
    assert_eq!(indices,vec![2,3,4,6,7,8,10,11,13,14,16,17,20,21,22,24,25,26]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_filter(){
    let filename = "./assets/a_TEMPO_a_water_a_glycerol.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&config).unwrap();

  
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Oxygen];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![28,32,36,41,43]);

    let mut filter = ParticleFilter::new();
    filter.not_elements = vec![Element::Hydrogen,Element::Carbon];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![27,28,32,36,41,43]);

    let mut filter = ParticleFilter::new();
    filter.serials = vec![28,29];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![27,28]);

    let mut filter = ParticleFilter::new();
    filter.bonded_elements = vec![Element::Oxygen];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![27,29,33,34,37,38,42,44,45]);

    let mut filter = ParticleFilter::new();
    filter.not_bonded_elements = vec![Element::Carbon];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![28,33,37,42,43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.residues = vec!["SOL".to_string()];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.names = vec!["H15".to_string()];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![26]);

    let mut filter = ParticleFilter::new();
    filter.not_names = vec!["OW".to_string()];
    filter.residues = vec!["SOL".to_string()];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![44,45]);

    let mut filter = ParticleFilter::new();
    filter.not_residues = vec!["TEM".to_string(),"MGL".to_string()];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.residue_sequence_numbers = vec![1502];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.not_residue_sequence_numbers = vec![1,2];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![43,44,45]);

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    filter.bonded_elements = vec![Element::Oxygen];
    let indices = filter.filter(&structure);
    assert_eq!(indices,vec![33,37,42,44,45]);

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    filter.not_bonded_elements = vec![Element::Oxygen];
    let indices = filter.filter(&structure);
    assert_eq!(indices.len(),23);
    assert_eq!(indices,vec![2,3,4,6,7,8,10,11,13,14,16,17,20,21,22,24,25,26,
    30,31,35,39,40]);


    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen,Element::Nitrogen];
    filter.within_distance = Some(4e-10); 
    filter.not_within_distance = Some(1e-10); 
    let indices = filter.filter(&structure);
    assert_eq!(indices.len(),16);
    assert_eq!(indices,vec![2,3,4,6,7,8,10,11,14,16,20,21,22,24,25,26]);
  }
}
