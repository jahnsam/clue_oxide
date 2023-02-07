use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::physical_constants::{PI,Element,Isotope};
use crate::space_3d::Vector3D;
use crate::structure::particle::Particle;
use crate::structure::Structure;
use crate::vec_funcs;

use substring::Substring;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function reads the contents of a PDB file into a Vec<Structure>.
/// Each model within the PDB is saved as a separate Structure.
/// Only the infomation directly in the PDB is set.  
/// To get additional information, such as on spins and methyl groups,
/// feed the output of this function to build_primary_structure(). 
///
///```
///# use clue_oxide::structure::pdb::parse_pdb;
/// let filename = "./assets/TEMPO.pdb";
/// let file = std::fs::read_to_string(filename).unwrap();
/// let structures = parse_pdb(&file).unwrap();
///```
///
pub fn parse_pdb(file: &str) -> Result< Vec::<Structure>, CluEError>
{

  // Crystallographic Information 
  let cell_offsets = parse_pdb_crystal(&file)?;

  
  // Models
  let (model_indices,endmdl_indices) = parse_pdb_models(&file)?;
  let n_structs = model_indices.len();


  // Connections
  let mut connections_vec = Vec::<AdjacencyList>::new();


  let mut structures = Vec::<Structure>::with_capacity(n_structs);


  for (ii, mdl) in model_indices.iter().enumerate(){

    // Atoms
    let bath_particles 
      = parse_pdb_atoms(&file.substring(*mdl, endmdl_indices[ii] ))?;

    // Connections
    if connections_vec.is_empty(){
      let connections = parse_pdb_connections(&file,&bath_particles)?;
      connections_vec.push(connections);
    }

    structures.push(Structure::new(
        bath_particles,
        connections_vec[0].clone(),
        cell_offsets.clone(),
        ));
  }

  Ok(structures)
}
//------------------------------------------------------------------------------
fn parse_pdb_atoms(file: &str) -> Result< Vec::<Particle>, CluEError> {


  let mut atom_indices: Vec::<usize> 
    = file.match_indices("ATOM").map(|(idx,_substr)| idx).collect();
  let mut hetatm_indices: Vec::<usize> 
    = file.match_indices("HETATM").map(|(idx,_substr)| idx).collect();
  atom_indices.append(&mut hetatm_indices);
  atom_indices = vec_funcs::unique(atom_indices);
  atom_indices.sort();


  // Read atoms.
  let n_atoms = atom_indices.len();
  let mut bath_particles = Vec::<Particle>::with_capacity(n_atoms);

  for idx in atom_indices.iter() {
    let line = get_line(file.substring(*idx,file.len()));
    let particle = parse_atom(line)?;
    bath_particles.push(particle)
  }


  Ok(bath_particles)
}
//------------------------------------------------------------------------------
fn parse_pdb_connections(file: &str,bath_particles: &[Particle]) 
  -> Result< AdjacencyList, CluEError> 
{

  //let n_atoms = file.matches("ATOM").count() + file.matches("HETATM").count();
  let n_atoms = bath_particles.len();

  // Read connections.
  let conect_indices: Vec::<usize> 
    = file.match_indices("CONECT").map(|(idx,_substr)| idx).collect();

  let mut connections = AdjacencyList::with_capacity(n_atoms);
  for line_idx in conect_indices.iter() {
    let line = get_line(file.substring(*line_idx,file.len()));
    match parse_connections(line){
   
      Ok(serials) => {

        let serial0 = vec![serials[0]];
        let index0 = serials_to_indices(serial0,&bath_particles);
        let index0 = index0[0];

        let indices = serials_to_indices(serials,&bath_particles);


        for idx in indices.iter(){
          connections.connect(index0,*idx);
        }
      },

      Err(error) => return Err(error),
    }
  }

  Ok(connections)
}
//------------------------------------------------------------------------------
fn parse_pdb_crystal(file: &str)
  -> Result< Vec::<Vector3D>, CluEError> 
{
  let cryst1_indices: Vec::<usize> 
    = file.match_indices("CRYST1").map(|(idx,_substr)| idx).collect();
  
  if !cryst1_indices.is_empty() {
    let line = get_line(file.substring(cryst1_indices[0],file.len()));
    return parse_crystal_line(line);
  }
    Ok(Vec::<Vector3D>::new())
  
}
//------------------------------------------------------------------------------
fn parse_pdb_models(file: &str)
  -> Result< (Vec::<usize>,Vec::<usize>), CluEError> 
{
  let mut model_indices: Vec::<usize> 
    = file.match_indices("MODEL").map(|(idx,_substr)| idx).collect();

  let mut endmdl_indices: Vec::<usize> 
    = file.match_indices("ENDMDL").map(|(idx,_substr)| idx).collect();

  assert_eq!(model_indices.len(), endmdl_indices.len());
  if model_indices.is_empty(){
    model_indices.push(0);
    endmdl_indices.push(file.len());
  }

  Ok((model_indices,endmdl_indices))
}
//------------------------------------------------------------------------------
fn get_line<'a>(file_slice: &'a str) -> &'a str {
  
  for (ii,ch) in file_slice.chars().enumerate() {

    if ch == '\n' {
      return &file_slice[0..ii+1];
    }
  }

  file_slice
}
//------------------------------------------------------------------------------
/*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/
fn parse_atom(line: &str) -> Result<Particle,CluEError>{

  let serial_result = line[6..=10].trim().parse::<u32>();
  let x_coor = line[30..=37].trim().parse::<f64>();
  let y_coor = line[38..=45].trim().parse::<f64>();
  let z_coor = line[46..=53].trim().parse::<f64>();
  let element = Element::from(line[76..=77].trim())?;
  let residue = Some(line[17..=19].trim().to_string());
  let residue_sequence_number = line[22..=25].trim().parse::<u32>();
 


  let coordinates: Vector3D;
  if let (Ok(x),Ok(y), Ok(z)) = (x_coor, y_coor, z_coor) {
    coordinates = Vector3D::from([x,y,z]);
  }else{
    return Err(CluEError::CannotParseLine(line.to_string()));
  }

  let serial: u32;
  if let Ok(ser) = serial_result{
    serial = ser;
  }else{
    return Err(CluEError::CannotParseLine(line.to_string()));
  }

  Ok(Particle{
      element,
      isotope: Isotope::most_common_for(&element),
      coordinates,
      serial: Some(serial),
      residue,
      residue_sequence_number: residue_sequence_number.ok(),
      //exchange_group: None,
      //exchange_group_id: None,
      })
}
//------------------------------------------------------------------------------
/*
COLUMNS       DATA  TYPE      FIELD        DEFINITION
-------------------------------------------------------------------------
 1 -  6        Record name    "CONECT"
 7 - 11        Integer        serial       Atom  serial number
12 - 16        Integer        serial       Serial number of bonded atom
17 - 21        Integer        serial       Serial  number of bonded atom
22 - 26        Integer        serial       Serial number of bonded atom
27 - 31        Integer        serial       Serial number of bonded atom
*/
fn parse_connections(conect_line: &str) -> Result<Vec::<u32>,CluEError>{
  let line = conect_line[6..].trim_end();

  let num = line.len();
  let n = 5;
  let mut serials = Vec::<u32>::with_capacity(num/n);
  let mut idx = 0;
  while idx < num {
    let idx_end = usize::min(idx+n,num);
    let connect: u32; 
    if let Ok(serial) = line[idx .. idx_end].trim().parse::<u32>(){
      connect = serial;
    }else{
      return Err(CluEError::CannotParseLine(conect_line.to_string() ));
    }
    serials.push(connect);
    idx += n;
  }

  Ok(serials)

}
//------------------------------------------------------------------------------
fn serials_to_indices(serials: Vec::<u32>,bath_particles: &[Particle])
 -> Vec::<usize>
{
  let mut indices = Vec::<usize>::with_capacity(serials.len());
  for (idx,particle) in bath_particles.iter().enumerate(){
    if let Some(serial) = (*particle).serial{
      if serials.contains(&serial){
        indices.push(idx);
      }
    }
  }
  indices
}
//------------------------------------------------------------------------------
/*
COLUMNS       DATA  TYPE    FIELD          DEFINITION
-------------------------------------------------------------
 1 -  6       Record name   "CRYST1"
 7 - 15       Real(9.3)     a              a (Angstroms).
16 - 24       Real(9.3)     b              b (Angstroms).
25 - 33       Real(9.3)     c              c (Angstroms).
34 - 40       Real(7.2)     alpha          alpha (degrees).
41 - 47       Real(7.2)     beta           beta (degrees).
48 - 54       Real(7.2)     gamma          gamma (degrees).
56 - 66       LString       sGroup         Space  group.
67 - 70       Integer       z              Z value.
*/
fn parse_crystal_line(line: &str) 
  -> Result<Vec::<Vector3D>,CluEError>
{
  if let (Ok(a),Ok(b),Ok(c),Ok(mut alpha),Ok(mut beta),Ok(mut gamma)) = (
      line[6..=14].trim().parse::<f64>(),
      line[15..=23].trim().parse::<f64>(),
      line[24..=32].trim().parse::<f64>(),
      line[33..=39].trim().parse::<f64>(),
      line[40..=46].trim().parse::<f64>(),
      line[47..=53].trim().parse::<f64>() ){

    alpha *= PI/180.0;
    beta *= PI/180.0;
    gamma *= PI/180.0;

    let a_vec = Vector3D::from([a, 0.0, 0.0]);
    let b_vec = Vector3D::from([b*f64::cos(gamma), b*f64::sin(gamma), 0.0]);
    let cx = c * f64::cos(beta);
    let cy = c*(f64::cos(alpha) - f64::cos(beta)*f64::cos(gamma))/f64::sin(gamma);
    let cz = c * f64::sqrt( 1.0 -cx*cx - cy*cy);
    let c_vec = Vector3D::from([cx, cy, cz]);
  

    return Ok(vec![a_vec,b_vec,c_vec]);
  }
  Err(CluEError::CannotParseLine(line.to_string() ))
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests {
  use super::*;
  #[test]
  fn test_parse_pdb(){
    let filename = "./assets/TEMPO.pdb";
    let file = std::fs::read_to_string(filename).unwrap();
    let structures = parse_pdb(&file).unwrap();

    assert_eq!(structures.len(),1);
    assert_eq!(structures[0].bath_particles.len(),29);
    assert_eq!(structures[0].bath_particles[27].element,Element::Nitrogen);
    assert_eq!(structures[0].bath_particles[28].element,Element::Oxygen);
    assert!(structures[0].connections.are_connected(27,28));

    assert_eq!(structures[0].bath_particles[27].serial, Some(28));
    assert_eq!(structures[0].bath_particles[27].residue,
        Some("TEM".to_string()));
    assert_eq!(structures[0].bath_particles[27].residue_sequence_number, 
        Some(1));
    assert_eq!(structures[0].bath_particles[27].coordinates.x(), 36.440);
    assert_eq!(structures[0].bath_particles[27].coordinates.y(), 36.900);
    assert_eq!(structures[0].bath_particles[27].coordinates.z(), 37.100);

    assert_eq!(structures[0].cell_offsets[0].x(),72.5676);
    assert_eq!(structures[0].cell_offsets[0].y(),0.0);
    assert_eq!(structures[0].cell_offsets[0].z(),0.0);

    assert!((structures[0].cell_offsets[1].x()-0.0).abs() < 1e12);
    assert_eq!(structures[0].cell_offsets[1].y(),72.5676);
    assert!((structures[0].cell_offsets[1].z()-0.0).abs() < 1e12);

    assert!((structures[0].cell_offsets[2].x()-0.0).abs() < 1e12);
    assert!((structures[0].cell_offsets[2].y()-0.0).abs() < 1e12);
    assert_eq!(structures[0].cell_offsets[2].z(),72.5676);


    //         O
    //   3HC   N  CH3
    //  3HC-C    C-CH3
    //    2HC    CH2 
    //        CH2

    assert!(structures[0].connections.are_connected(0,1)); 
    assert!(structures[0].connections.are_connected(0,5));
    assert!(structures[0].connections.are_connected(0,9));
    assert!(structures[0].connections.are_connected(0,27));
    assert!(structures[0].connections.are_connected(1,2));
    assert!(structures[0].connections.are_connected(1,3));
    assert!(structures[0].connections.are_connected(1,4));
    assert!(structures[0].connections.are_connected(5,6));
    assert!(structures[0].connections.are_connected(5,7));
    assert!(structures[0].connections.are_connected(5,8));
    assert!(structures[0].connections.are_connected(9,10));
    assert!(structures[0].connections.are_connected(9,11));
    assert!(structures[0].connections.are_connected(9,12));
    assert!(structures[0].connections.are_connected(12,13));
    assert!(structures[0].connections.are_connected(12,14));
    assert!(structures[0].connections.are_connected(12,15));
    assert!(structures[0].connections.are_connected(15,16));
    assert!(structures[0].connections.are_connected(15,17));
    assert!(structures[0].connections.are_connected(15,18));
    assert!(structures[0].connections.are_connected(18,19));
    assert!(structures[0].connections.are_connected(18,23));
    assert!(structures[0].connections.are_connected(18,27));
    assert!(structures[0].connections.are_connected(19,20));
    assert!(structures[0].connections.are_connected(19,21));
    assert!(structures[0].connections.are_connected(19,22));
    assert!(structures[0].connections.are_connected(23,24));
    assert!(structures[0].connections.are_connected(23,25));
    assert!(structures[0].connections.are_connected(23,26));
    assert!(structures[0].connections.are_connected(27,28));

    assert!(!structures[0].connections.are_connected(0,28));
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
