use super::config::*;
use super::clue_errors::*;
use super::pdb::PDB;
use super::vector3::Vector3;

#[derive(Debug,Clone)]
pub struct PrimaryStructure{
  //exchange_groups: Vec<ExchangeGroup>,
  //exchange_group_ids: Vec<Option<usize>>,
  pdb_indices: Vec<usize>,
  pdb: PDB,
}

/*
#[derive(Debug,Clone)]
pub struct ExchangeGroup{
  indices: Vec<usize>,
  exchange_coupling: f64,
}
*/
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl PrimaryStructure {
// This function takes in a pdb, and
pub fn from(mut pdb: PDB, config: &Config) -> Result<Self,CluEError> {

  let origin = get_electron_coordinates(&pdb, &config)?;
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

  //let exchange_groups = ExchangeGroup::from( pdb.find_methyls(), &config );
  
  //let exchange_group_ids = find_exchange_group_ids(
  //    &pdb_indices, &exchange_groups);

  Ok(PrimaryStructure{
      //exchange_groups,
      //exchange_group_ids,
      pdb_indices,
      pdb,
  })
}
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
impl ExchangeGroup{
  pub fn from(methyls: Methyls, config: &Config) -> Vec<Self>{

      for pdb_idx in exchange_groups[ii].indices{
        let spin_idx = pdb_index_to_spin_index(pdb_idx);
      }
  }
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
fn find_exchange_group_ids(
    pdb_indices: &Vec<usize>, exchange_groups: Vec<ExchangeGroup>) 
  -> Vec<Option<usize>>{

    let number = pdb_indices.len();
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
*/
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
fn get_electron_coordinates(pdb: &PDB, config: &Config) 
  -> Result<Vector3, CluEError>{

  let coor = config.central_spin_coordinates.clone();

  match coor{
    Some(CentralSpinCoordinates::Atoms(atoms)) => Ok(pdb.get_centroid(&atoms)),
    Some(CentralSpinCoordinates::XYZ(r)) => Ok(r), 
    _ => Err(CluEError::NoCentralSpinCoor)
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use super::super::pdb::*;
  use super::super::physical_constants::*;

  #[test]
  #[allow(non_snake_case)]
  fn test_PrimaryStructure() {
    let filename = "./assets/TEMPO.pdb";
    let mut pdb = read_pdb(filename).expect("Could not read pdb file.");
    let r0 = pdb.get_centroid(&vec![28,29]);
    pdb.set_origin(&r0);

    let mut config = Config::new();


    config.central_spin_coordinates = Some(
        CentralSpinCoordinates::Atoms(vec![28,29]) );

    let structure = PrimaryStructure::from(pdb, &config)
      .unwrap_or_else(|e| { panic!("error = {:?}",e);});

    assert_eq!(structure.pdb_indices.len(),19);
    assert_eq!(structure.pdb.element(structure.pdb_indices[0]),
        Element::Hydrogen);
    assert_eq!(structure.pdb.element(structure.pdb_indices[18]),
        Element::Nitrogen);

  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


