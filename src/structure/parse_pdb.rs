use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::physical_constants::{ANGSTROM,PI,Element,Isotope};
use crate::structure::{Structure, particle::Particle};
use crate::space_3d::Vector3D;

use std::io::BufRead;
use std::collections::HashMap;

pub fn parse_pdb(filename: &str) -> Result< Vec::<Structure>, CluEError>{

   let (n_atoms, n_models) = count_pdb_atoms(filename)?;
   
   let (bath_particles_sets, serial_to_index) 
     = parse_atoms(filename, n_atoms,n_models)?;
   let mut connections = AdjacencyList::with_capacity(n_atoms);
   let mut cell_offsets = Vec::<Vector3D>::new();
   
   let Ok(file) = std::fs::File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };
   let lines = std::io::BufReader::new(file).lines();


   for line_result in lines {
     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };
     if line.contains("CONECT"){
       let a = parse_connections(&line);
     }else if line.contains("CRYST1"){
       let a = parse_crystal_line(&line)?;
     }
   }
    

  let mut structures = Vec::<Structure>::with_capacity(n_models);
  for bath_particles in bath_particles_sets{
    structures.push(Structure::new(
      bath_particles,
      connections.clone(),
      cell_offsets.clone())
    )
  }
  Ok(structures)
}
//------------------------------------------------------------------------------
fn parse_atoms(filename: &str,n_atoms: usize, n_models: usize) 
  -> Result<(Vec::<Vec::<Particle>>, HashMap::<u32,usize>),CluEError>
{
   let mut serial_to_index = HashMap::<u32,usize>::new(); 
   let mut bath_particles = Vec::<Vec::<Particle>>::with_capacity(n_models);
   for ii in 0..n_models{
     bath_particles.push( Vec::<Particle>::with_capacity(n_atoms));
   }

   let mut model_idx = 0;

   let Ok(file) = std::fs::File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };
   
   let lines = std::io::BufReader::new(file).lines();
   for line_result in lines {
     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };
     if line.contains("ATOM") || line.contains("HETATM") {
       let particle = parse_atom(&line)?;

       if model_idx == 0{
         if let Some(serial) = particle.serial{
           if !serial_to_index.contains_key(&serial){
             serial_to_index.insert(serial,bath_particles[0].len());
           }
         }
       }

       bath_particles[model_idx].push(particle);
     }else if line.contains("ENDMDL"){
       model_idx += 1;
     }
   }

   Ok((bath_particles,serial_to_index))
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

    let a = a*ANGSTROM;
    let b = b*ANGSTROM;
    let c = c*ANGSTROM;

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
    coordinates = Vector3D::from([x*ANGSTROM,y*ANGSTROM,z*ANGSTROM]);
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
fn count_pdb_atoms(filename: &str) -> Result<(usize,usize),CluEError> {

   let mut count = 0;
   let mut model_count = 0;
   let mut increment = 1;
   let Ok(file) = std::fs::File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };
   let lines = std::io::BufReader::new(file).lines();
   for line_result in lines {
     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };
     if line.contains("MODEL"){
       model_count += 1;
     }else if line.contains("ENDMDL"){
       increment = 0;
     }else if line.contains("ATOM") || line.contains("HETATM") {
       count += increment;
     }
   }

   if model_count==0{
     model_count += 1;
   }
   Ok((count,model_count))
}


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests {
  use super::*;
  #[test]
  fn test_parse_pdb(){
    let filename = "./assets/TEMPO_wat_gly_70A.pdb";
    let structures = parse_pdb(&filename).unwrap();
  }
}
