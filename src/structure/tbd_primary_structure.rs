use crate::config::*;
use crate::clue_errors::*;
use crate::structure::pdb::PDB;
use crate::vector3::Vector3;

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

  let origin = get_electron_coordinates(&pdb, config)?;
  pdb.set_origin(&origin);

  let number = determine_max_number_of_spinful_particles(&pdb,config);

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
//------------------------------------------------------------------------------
fn reconnect_bonds(&mut self){

    const MAXBOND: f64 = 3.0;
    for connections in self.connections.iter() {

      let mut it = connections.indices.iter();
      let idx0: &usize = it.next().unwrap();

      let r0 = self.coordinates[*idx0].clone();

      for idx in it{

        let r = self.coordinates[*idx].clone();

        if (&r-&r0).magnitude() < MAXBOND { continue; }

        self.coordinates[*idx].x
         = minimize_absolute_difference_for_step(r.x, r0.x, self.crystal.a);

        self.coordinates[*idx].y
         = minimize_absolute_difference_for_step(r.y, r0.y, self.crystal.b);

        self.coordinates[*idx].z
         = minimize_absolute_difference_for_step(r.z, r0.z, self.crystal.c);
      }





    }
 }
//------------------------------------------------------------------------------
  fn get_methyl_hydrogen_indices(&self, index: usize) ->Option<[usize;3] >{

    if self.elements[index] != Element::Carbon
     && self.elements[index] != Element::Nitrogen {
      return None;
    }

    let connection_idx: usize;
    if let Some(idx)
      = self.get_connections_index_from_serial(self.serials[index]){
        connection_idx = idx;
      }else{
        return None;
    }

    if self.connections[connection_idx].serials.len() != 5 {
      return None;
    }

    let mut hydrons = Vec::<usize>::with_capacity(3);

    for ii in 1..self.connections[connection_idx].indices.len() {
      let idx = self.connections[connection_idx].indices[ii];

      if self.elements[idx] == Element::Hydrogen{
        if hydrons.len()==3{ return None;}
        hydrons.push(idx);
      }
    }


    if hydrons.len()!=3{ return None;}

    Some([hydrons[0], hydrons[1], hydrons[2] ])
  }
//------------------------------------------------------------------------------
  fn set_connection_indices(&mut self){

    for ipdb in 0..self.number {
      let n_serials = self.connections[ipdb].serials.len();
      for index in 0..n_serials{
        let serial = self.connections[ipdb].serials[index];
          let pdb_index = self.find_index(serial)
            .expect("Could not find index.");
          self.connections[ipdb].indices.push(pdb_index);
      }
    }

  }
//------------------------------------------------------------------------------
  fn get_centroid(&self, serials: &[u32]) -> Vector3 {

    let mut centroid = Vector3::new();
    for id in serials {
      let idx = self.find_index(*id).expect("Cannot find atom.");
      centroid = &centroid + &self.coordinates[idx];
    }

    centroid.scale(1.0/serials.len() as f64)
  }
//------------------------------------------------------------------------------
  fn find_exchange_groups(&self) -> Vec::<ExchangeGroup> {

    let n_max: usize = self.number/5;

    let mut exchange_groups = Vec::<ExchangeGroup>::with_capacity(n_max);

    for ii in 0..self.number{

      let hydrogens = self.get_methyl_hydrogen_indices(ii);
       if let Some([h0,h1,h2]) = hydrogens {
         let r_carbon = self.coordinates[ii].clone();
         let r_h0 = self.coordinates[h0].clone();
         let r_h1 = self.coordinates[h1].clone();
         let r_h2 = self.coordinates[h2].clone();


         if self.elements[ii] == Element::Carbon{
           exchange_groups.push(ExchangeGroup::Methyl(
                 C3Rotor::from(r_carbon, r_h0, r_h1, r_h2, [h0,h1,h2]) ));

         }else if self.elements[ii] == Element::Nitrogen{
           exchange_groups.push(ExchangeGroup::PrimaryAmonium(
                 C3Rotor::from(r_carbon, r_h0, r_h1, r_h2, [h0,h1,h2]) ));
         }
       }
    }

    exchange_groups
  }

//------------------------------------------------------------------------------
pub fn minimize_absolute_difference_for_step( x: f64, x0: f64, step: f64 )
 -> f64 {

   assert!(step > 0.0);

   let mut x1 = x;

   loop{

     if (x1 - x0).abs() <= (x1 + step - x0).abs()
       && (x1 - x0).abs() <= (x1 - step - x0).abs(){
         break;
       }
     else if (x1 - x0).abs() > (x1 + step - x0).abs(){
      x1 += step;
     }
     else if (x1 - x0).abs() > (x1 - step - x0).abs(){
      x1 -= step;
     }
   }

  x1
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb::*;
  use crate::physical_constants::*;

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

  //----------------------------------------------------------------------------
  #[test]
  fn test_find_exchange_groups(){
    let filename = "./assets/TEMPO.pdb";
    let mut pdb = read_pdb(filename).expect("Could not read pdb file.");

    let r0 = pdb.get_centroid(&vec![28,29]);
    pdb.set_origin(&r0);


    assert_eq!(pdb.exchange_groups.len(),4);

    let distances = [3.0352,   2.9653,   3.0415,   2.8101];
    for ii in 0..pdb.exchange_groups.len(){
      if let ExchangeGroup::Methyl(methyl) = &pdb.exchange_groups[ii] { 
        let r = methyl.center().magnitude();
        assert!( (r-distances[ii]).abs()<1e-4 )
      }else{
        panic!("There should only be methyls for exchange groups.")
      }
    }
  }
  
  #[test]
  fn test_minimize_absolute_difference_for_step(){

    assert_eq!(minimize_absolute_difference_for_step(1.2,1.3,1.0),1.2);
    assert!((minimize_absolute_difference_for_step(2.2,1.3,1.0)-1.2).abs()
        < 1e-12);
    assert!((minimize_absolute_difference_for_step(-2.8,1.3,1.0)-1.2).abs()
        < 1e-12);

    assert!((minimize_absolute_difference_for_step(1.0,1.5,1.0)-1.0).abs()
        < 1e-12);
  }  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


