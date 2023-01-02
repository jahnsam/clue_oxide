use std::error::Error;
use std::fs;
use std::io::{self, BufRead};
use std::path::Path;

use crate::vector3::Vector3;
use crate::structure::exchange_groups::*;
use crate::physical_constants::*;

#[derive(Debug, Clone)]
pub struct PDB{
  number: usize,
  crystal: CrystalInfo,
  serials: Vec<u32>,
  residues: Vec<String>,
  elements: Vec<Element>,
  coordinates: Vec<Vector3>,
  connections: Vec<ConnectionInfo>,
  residue_sequence_numbers: Vec<u32>,
  exchange_groups: Vec::<ExchangeGroup>,
  exchange_group_ids: Vec::<Option<usize>>,
  solvent_exchangabilities: Vec::<bool>,
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn read_pdb(filename: &str) -> Result<PDB, Box<dyn Error> >{

    let n_atoms = count_pdb_atoms(filename)?;
    let mut pdb = PDB::new(n_atoms as usize);

    // Read lines.
    if let Ok(lines) = read_lines(filename) {
        // Consumes the iterator and returns an Option<String>.
        for line in lines.flatten() {
          let linetype = parse_line(&line[..]);

          match linetype {
            LineType::Coordinate(coor_info) =>{
              pdb.serials.push(coor_info.serial);
              pdb.residues.push(coor_info.residue);
              pdb.elements.push(coor_info.element);
              pdb.coordinates.push(
                    coor_info.coordinates,
              );
              pdb.residue_sequence_numbers.push(
                  coor_info.residue_sequence_number);

              pdb.exchange_group_ids.push(None);
            },
            LineType::Crystal(crystal) => {
              pdb.crystal = crystal;
            },
            LineType::Connection(connection_info) => {
              pdb.connections.push(connection_info);
            }
            _ => (),
          }
        }
    }

    pdb.set_connection_indices();

    // Set exchange groups.
    pdb.exchange_groups = pdb.find_exchange_groups();
    
    for (ii, egroup) in pdb.exchange_groups.iter().enumerate() {
      let indices = egroup.indices();
      for h in  indices{
        pdb.exchange_group_ids[h] = Some(ii);
      }
    }
    

  Ok(pdb)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl PDB{
  //----------------------------------------------------------------------------
  pub fn new(n_atoms: usize) -> PDB {
  
    PDB {
      number: n_atoms,
      crystal: CrystalInfo{ a: 0.0,
        b: 0.0,
        c: 0.0,
        alpha: 0.0,
        beta: 0.0,
        gamma: 0.0,
        cell_offsets: [Vector3::new(), Vector3::new(), Vector3::new() ],
        },
      elements: Vec::<Element>::with_capacity(n_atoms),
      serials: Vec::<u32>::with_capacity(n_atoms),
      residues: Vec::<String>::with_capacity(n_atoms),
      coordinates: Vec::<Vector3>::with_capacity(n_atoms),
      connections: Vec::<ConnectionInfo>::with_capacity(n_atoms),
      residue_sequence_numbers: Vec::<u32>::with_capacity(n_atoms),
      exchange_groups: Vec::<ExchangeGroup>::new(),
      exchange_group_ids: Vec::<Option<usize>>::with_capacity(n_atoms),
      solvent_exchangabilities: Vec::<bool>::with_capacity(n_atoms),
    }
  }
    
  //----------------------------------------------------------------------------
  pub fn print(&self) {
  
    println!("PDB {{");
    println!("  number: {}\n", self.number);
    println!("  {:?}",self.crystal);

    for ii in 0..self.number{
      println!("    {}  {:?}  {}  {}  {}", 
          self.serials[ii], self.elements[ii], 
          self.coordinates[ii].x, 
          self.coordinates[ii].y, 
          self.coordinates[ii].z);
    }
    println!("\n  connections serials:");
    for ii in 0..self.connections.len(){
      println!("    {:?}", self.connections[ii].serials);
    }
    println!("}}");
  }
  
  
  //----------------------------------------------------------------------------
  pub fn number(&self) -> usize{
    self.serials.len()
  }
  
  //----------------------------------------------------------------------------
  pub fn coordinates(&self,n: usize) -> Vector3{
    self.coordinates[n].clone()
  }
  //----------------------------------------------------------------------------
  pub fn serial(&self, n: usize) -> u32{
    self.serials[n]
  }
  //----------------------------------------------------------------------------
  pub fn residue(&self, n: usize) -> String{
    self.residues[n].clone()
  }
  //----------------------------------------------------------------------------
  pub fn element(&self, n: usize) -> Element {
    self.elements[n]
  }
  //----------------------------------------------------------------------------
  pub fn residue_sequence_number(&self, n: usize) -> u32{
    self.residue_sequence_numbers[n]
  }
  //----------------------------------------------------------------------------
  pub fn exchange_group(&self, n: usize) -> Option<&ExchangeGroup>{
    if let Some(idx) = self.exchange_group_ids[n]{
      return Some(&self.exchange_groups[idx]);
    }
    None
  }
  //----------------------------------------------------------------------------
  pub fn crystal_a(&self) -> f64 {self.crystal.a}
  pub fn crystal_b(&self) -> f64 {self.crystal.b}
  pub fn crystal_c(&self) -> f64 {self.crystal.c}
  pub fn crystal_alpha(&self) -> f64 {self.crystal.alpha}
  pub fn crystal_beta(&self) -> f64 {self.crystal.beta}
  pub fn crystal_gamma(&self) -> f64 {self.crystal.gamma}
  //----------------------------------------------------------------------------
  pub fn find_index(&self, serial: u32) -> Option<usize> {
  
    for (index, value) in self.serials.iter().enumerate(){
      if serial == *value {
        return Some(index);
      }
    }

    None
  }
  //----------------------------------------------------------------------------
  pub fn find_exchange_groups(&self) -> Vec::<ExchangeGroup> {
  
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
  //----------------------------------------------------------------------------
  pub fn get_centroid(&self, serials: &[u32]) -> Vector3 {
  
    let mut centroid = Vector3::new();
    for id in serials {
      let idx = self.find_index(*id).expect("Cannot find atom.");
      centroid = &centroid + &self.coordinates[idx];
    }

    centroid.scale(1.0/serials.len() as f64)
  }
  //----------------------------------------------------------------------------
  
  pub fn translate(&mut self, r: &Vector3){
    for ii in 0..self.number{
      self.coordinates[ii] = &self.coordinates[ii] + r;
    }
    for egroup in &mut self.exchange_groups{
      egroup.translate(r); 
    }
  }
  
  //----------------------------------------------------------------------------
  pub fn set_origin(&mut self, r: &Vector3){
    self.translate(&r.scale(-1.0));
  }
  //----------------------------------------------------------------------------
  fn get_connections_index_from_serial(&self, serial: u32) 
    -> Option<usize>
    {
    for ii in 0..self.connections.len(){
      if !self.connections[ii].serials.is_empty() 
        && self.connections[ii].serials[0] == serial{
          return Some(ii);
      }
      
    }
    None
  }
  //----------------------------------------------------------------------------
  
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
 
  
  
  
  //----------------------------------------------------------------------------
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
 //----------------------------------------------------------------------------- 
 pub fn reconnect_bonds(&mut self){

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
 //----------------------------------------------------------------------------- 
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#[derive(Debug, Clone)]
struct CrystalInfo{
  a: f64,
  b: f64,
  c: f64,
  alpha: f64,
  beta: f64,
  gamma: f64,
  cell_offsets: [Vector3; 3],
}

#[derive(Debug, Clone)]
struct CoordinateInfo {
  serial: u32,
  residue: String,
  coordinates: Vector3,
  element: Element,
  residue_sequence_number: u32,
}

#[derive(Debug, Clone)]
struct ConnectionInfo {
  serials: Vec<u32>,
  indices: Vec<usize>,
}

enum LineType {
  Crystal(CrystalInfo),
  Coordinate(CoordinateInfo),
  Connection(ConnectionInfo),
  Other,
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn count_pdb_atoms(filename: &str) -> Result<u32,Box<dyn Error>> {

    let mut count: u32 = 0;
    let lines = read_lines(filename)?;
        for line in lines.flatten() {
          if line.contains("ATOM") || line.contains("HETATM") {
            count += 1;
          }
        }
    Ok(count)
}
  //----------------------------------------------------------------------------
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<fs::File>>>
where P: AsRef<Path>, {
    let file = fs::File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
//------------------------------------------------------------------------------

pub fn parse_config(args: &[String]) -> Result<String,&'static str>{
  if args.len() < 2{
    return Err("Please provide an input file.");
  }
  Ok(args[1].clone())
}
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
fn read_coordinate_line(line: &str) -> LineType {

  let serial: u32 = line[6..=10].trim().parse().expect("no serial");
  let x: f64 = line[30..=37].trim().parse().expect("no x");
  let y: f64 = line[38..=45].trim().parse().expect("no y");
  let z: f64 = line[46..=53].trim().parse().expect("no z");
  let element = Element::from(line[76..=77].trim()).unwrap();
  let residue: String = line[17..=19].trim().parse().expect("no residue");
  let residue_sequence_number: u32 = line[22..=25].trim().parse()
    .expect("no residue sequence number");
 
  let coordinates = Vector3::from([x,y,z]);

  let coor_info = CoordinateInfo{
    serial,
    residue,
    coordinates,
    element,  
    residue_sequence_number,
  };

  LineType::Coordinate(coor_info)
}

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
fn read_crystal_line(line: &str) -> LineType {
  let a: f64 = line[6..=14].trim().parse().expect("no a");
  let b: f64 = line[15..=23].trim().parse().expect("no b");
  let c: f64 = line[24..=32].trim().parse().expect("no c");
  let alpha: f64 = line[33..=39].trim().parse().expect("no alpha");
  let beta: f64 = line[40..=46].trim().parse().expect("no beta");
  let gamma: f64 = line[47..=53].trim().parse().expect("no gamma");


  let a_vec = Vector3{x: a, y: 0.0, z: 0.0};
  let b_vec = Vector3{x: b*f64::cos(gamma), y: b*f64::sin(gamma), z: 0.0};
  let cx = c * f64::cos(beta);
  let cy = c*(f64::cos(alpha) - f64::cos(beta)*f64::cos(gamma))/f64::sin(gamma);
  let cz = c * f64::sqrt( 1.0 -cx*cx - cy*cy);
  let c_vec = Vector3{x: cx, y: cy, z: cz};
  
  let cryst_info = CrystalInfo{
  a,
  b,
  c,
  alpha,
  beta,
  gamma,
  cell_offsets: [a_vec, b_vec, c_vec],
  };

  LineType::Crystal(cryst_info)
}


/*
COLUMNS       DATA  TYPE      FIELD        DEFINITION
-------------------------------------------------------------------------
 1 -  6        Record name    "CONECT"
 7 - 11       Integer        serial       Atom  serial number
12 - 16        Integer        serial       Serial number of bonded atom
17 - 21        Integer        serial       Serial  number of bonded atom
22 - 26        Integer        serial       Serial number of bonded atom
27 - 31        Integer        serial       Serial number of bonded atom
*/
fn read_connetions_line(line: &str) -> LineType {
  let line = line[6..].trim_end();

  let num = line.len();
  let n = 5;
  let mut serials = Vec::<u32>::with_capacity(num/n);
  let mut idx = 0;
  while idx < num {
    let idx_end = usize::min(idx+n,num);
    let connect: u32 = line[idx .. idx_end].trim().parse()
      .unwrap_or_else(|_| panic!("Cannot parse {}--{} of line: '{}'.",
            idx,idx_end, line));
    serials.push(connect);
    idx += n;
  }

  let indices = Vec::<usize>::with_capacity(serials.len());

  LineType::Connection(ConnectionInfo{
    serials,
    indices,
  })

}
//----------------------------------------------------------------------------
fn parse_line(line: &str) -> LineType {

  if line.len() < 4 {
    return LineType::Other;
  }

  let model = line[0..6].trim();
  match model {
    "ATOM" | "HETATM" => read_coordinate_line(line),
    "CONECT" => read_connetions_line(line),
    "CRYST1" => read_crystal_line(line),
    _ => LineType::Other
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
mod tests {
  use super::*;
  #[test]
  fn test_read_pdb(){
    let filename = "./assets/TEMPO.pdb";
    let mut pdb = read_pdb(filename).expect("Could not read pdb file.");

    let r0 = pdb.get_centroid(&vec![28,29]);
    pdb.set_origin(&r0);
    
    let r_en0 = 7.1545e-01;
    let idx = pdb.find_index(28).unwrap();
    let r_en = pdb.coordinates[idx].magnitude();

    assert!((r_en-r_en0).abs()<1e-4);

    assert_eq!(pdb.number,29);
    assert_eq!(pdb.number,pdb.coordinates.len());
    assert_eq!(pdb.number,pdb.elements.len());
    assert_eq!(Element::Nitrogen,pdb.elements[27]);
    assert_eq!(pdb.number,pdb.serials.len());
    assert_eq!(pdb.number,pdb.residues.len());
    assert_eq!("TEM",pdb.residues[0] );
    assert_eq!(pdb.number,pdb.residue_sequence_numbers.len());
    assert_eq!(1 ,pdb.residue_sequence_numbers[0]);
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
