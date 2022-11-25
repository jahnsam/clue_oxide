use super::config::*;
use super::pdb::PDB;
use super::vector3::Vector3;


pub struct PrimaryStructure{
  exchange_groups: Vec<ExchangeGroup>,
  exchange_group_ids: Vec<Option<usize>>,
  pdb_indices: Vec<usize>,
  pdb: PDB,
}

pub struct ExchangeGroup{
  indices: Vec<usize>
  exchange_coupling: f64
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl PrimaryStructure {
// This function takes in a pdb, and
pub fn from(mut pdb: PDB, config: &Config) -> Self {

  let origin = get_electron_coordinates(&pdb, &config).unwrap();
  pdb.set_origin(&origin);

  let number = determine_max_number_of_spinful_particles(&pdb,&config);

  let mut pdb_indices = Vec::<usize>::with_capacity(number);


  for ipdb in 0..pdb.number(){
    
    
    let element = pdb.element(ipdb);
    if config.get_max_spin_multiplicity_for_any_isotope(element) == 1 {
      continue;
    }
    

    pdb_indices.push(ipdb);

  }

  let exchange_groups = ExchangeGroup::from( pdb.find_methyls(), &config );
  
  let exchange_group_ids = find_exchange_group_ids(
      &pdb_indices, &exchange_groups);

  PrimaryStructure{
      exchange_groups,
      exchange_group_ids,
      pdb_indices,
      pdb,
  }
}
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl ExchangeGroup{
  pub fn from(methyls; Methyls, config: &Config) -> Vec<Self>{

      for pdb_idx in exchange_groups[ii].indices{
        let spin_idx = pdb_index_to_spin_index(pdb_idx);
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn find_exchange_group_ids(
    pdb_indices: &Vec<usize>, exchange_groups: Vec<ExchangeGroup>) 
  -> Vec<Option<usize>>{

    let number = pdb_indices.len()
    let mut exchange_group_ids = Vec::<Option<usize>>::with_capacity();
    for _ii in 0..number{
      exchange_group_ids.push(None);
    }

    for ii in exchange_groups.len() {

      for spin_idx in exchange_groups[ii].indices{
        exchange_group_ids[spin_idx] = Some(ii);  

      }
    
    } 

}

//------------------------------------------------------------------------------

fn determine_max_number_of_spinful_particles(pdb: &PDB,config: &Config) 
  -> usize{
  let mut counter = 0;  
  for ipdb in 0..pdb.number(){
    let element = pdb.element(ipdb);
    if config.get_max_spin_multiplicity_for_any_isotope(element) > 1 {
      counter += 1;
    }
  }
  counter
}
//------------------------------------------------------------------------------
fn get_electron_coordinates(pdb: &PDB, config: &Config) -> Option<Vector3>{

  let coor = config.central_spin_coordinates.clone()?;

  match coor{
    CentralSpinCoordinates::Atoms(atoms) => Some(pdb.get_centroid(&atoms)),
    CentralSpinCoordinates::XYZ(r) => Some(r), 
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>







