use std::fs::File;
use std::error::Error;
use std::fs;
use std::io::{self, BufRead};
use std::path::Path;


pub struct PDB{
  number: u32,
  a: f64,
  b: f64,
  c: f64,
  alpha: f64,
  beta: f64,
  gamma: f64,
  serials: Vec<u32>,
  elements: Vec<String>,
  x: Vec<f64>,
  y: Vec<f64>,
  z: Vec<f64>,
  connections: Vec<ConnectionInfo>,
}


struct CrystalInfo{
  a: f64,
  b: f64,
  c: f64,
  alpha: f64,
  beta: f64,
  gamma: f64,
}

struct CoordinateInfo {
  serial: u32,
  x: f64,
  y: f64,
  z: f64,
  element: String,
}

struct ConnectionInfo {
  serials: Vec<u32>,
  indices: Vec<usize>,
}

enum LineType {
  CrystalLine(CrystalInfo),
  CoordinateLine(CoordinateInfo),
  ConnectionLine(ConnectionInfo),
  OtherLine,
}



impl PDB{
  
  pub fn write(&self, filename: &str) -> Result<(), Box<dyn Error> > {
    let mut output = std::fs::File::create(filename)?;


    Ok(())
  }
  pub fn new(n_atoms: usize) -> PDB {
  
    PDB {
      number: n_atoms as u32,
      a: 0.0,
      b: 0.0,
      c: 0.0,
      alpha: 0.0,
      beta: 0.0,
      gamma: 0.0,
      elements: Vec::<String>::with_capacity(n_atoms),
      serials: Vec::<u32>::with_capacity(n_atoms),
      x: Vec::<f64>::with_capacity(n_atoms),
      y: Vec::<f64>::with_capacity(n_atoms),
      z: Vec::<f64>::with_capacity(n_atoms),
    }
  }
    
  pub fn print(&self) {
  
    println!("PDB {}", "{");
    println!("  number: {}", self.number);
    println!("  a : {}, b : {}, c : {}, alpha : {}, beta : {}, gamma : {},",
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma);
    println!("  serials : {:?}", self.serials);
    println!("  elements : {:?}", self.elements);
    println!("  x : {:?}", self.x);
    println!("  y : {:?}", self.x);
    println!("  z : {:?}", self.x);
    println!("{}", "}");
  }
  
  fn count_pdb_atoms(filename: &str) -> Result<u32,Box<dyn Error>> {
  
      let mut count: u32 = 0;
      let lines = read_lines(filename)?;
          for line in lines {
            if let Ok(ip) = line {
              if ip.contains("ATOM") || ip.contains("HETATM") {
                count += 1;
              }
            }
          }
      Ok(count)
  }
  
  //----------------------------------------------------------------------------
  pub fn read_pdb(filename: String) -> Result<PDB, Box<dyn Error> >{
      // Read lines.
      println!("Reading {}.", filename);
  
      let n_atoms = PDB::count_pdb_atoms(&filename)?;
      let mut pdb = PDB::new(n_atoms as usize);
  
      if let Ok(lines) = read_lines(filename) {
          // Consumes the iterator, returns an (Optional) String
          for line in lines {
              if let Ok(ip) = line {
                let linetype = PDB::parse_line(&ip[..]);
  
                match linetype {
                  LineType::CoordinateLine(coor_info) =>{
                    pdb.serials.push(coor_info.serial);
                    pdb.elements.push(coor_info.element);
                    pdb.x.push(coor_info.x);
                    pdb.y.push(coor_info.y);
                    pdb.z.push(coor_info.z);
                  },
                  LineType::CrystalLine(cryst_info) => {
                    pdb.a = cryst_info.a;
                    pdb.b = cryst_info.b;
                    pdb.c = cryst_info.c;
                    pdb.alpha = cryst_info.alpha;
                    pdb.beta = cryst_info.beta;
                    pdb.gamma = cryst_info.gamma;
                  },
                  _ => (),
                }
              }
          }
      }
  
      pdb.set_connection_indices()
  
    Ok(pdb)
  }
  
  //----------------------------------------------------------------------------
  pub fn number(&self) -> u32{
    self.x.len().try_into().unwrap()
  }
  
  //----------------------------------------------------------------------------
  pub fn pos(&self,n: usize) -> [f64; 3]{
    [self.x[n],self.y[n], self.z[n] ]
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
    let element = line[76..=77].trim();
  
    let coor_info = CoordinateInfo{
      serial,
      x,
      y,
      z,
      element: element.to_string(),  
    };
  
    LineType::CoordinateLine(coor_info)
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
  
    let cryst_info = CrystalInfo{
    a,
    b,
    c,
    alpha,
    beta,
    gamma,
    };
  
    LineType::CrystalLine(cryst_info)
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
    let line = line[11..].trim_end();
    
    let n_id:usize = (line.len() + 5 - 1)/5;
    
    let mut serials = Vec::<u32>::with_capacity(n_id);
    let mut indices = Vec::<usize>::with_capacity(n_id);

    let mut line_idx = 0;

    while line_idx < line.len(){

      let new_idx = std::cmp::min(line_idx+5, line.len());

      let serial_id: u32 = line[line_idx..new_idx-1]
        .trim().parse().expect("Could not parse CONECT line.");

      serials.push(serial_id);
      
      // Initial guess for the index.
      indices.push((serial_id-1) as usize);
      
      line_idx = new_idx;
    }

    let ids = Vec::<u32>::with_capacity(n_id);

    ConnectionInfo{
      serials:
      indices:
    }
  
  }
  //----------------------------------------------------------------------------
  fn set_connection_indices(&mut self){

    for ipdb in 0..self.number {
      for (index, serial) in self.connections[ipdb].serials.iter().enumerate(){
        self.connections.indics[index] = self.find_index(serial)
          .expect("Could not find index."); 
      }
    }

  }
  //----------------------------------------------------------------------------
  pub fn find_index(&self, serial: u32) -> Option<usize> {
  
    for (index, value) in self.serial.iter().enumerate(){
      if serial == value {
        return Some(index);
      }
    }

    None
  }
  //----------------------------------------------------------------------------
  fn parse_line(line: &str) -> LineType {
  
    if line.len() < 4 {
      return LineType::OtherLine;
    }
  
    let model = line[0..6].trim();
    match model {
      "ATOM" | "HETATM" => PDB::read_coordinate_line(&line),
      "CONECT" => LineType::OtherLine,
      "CRYST1" => PDB::read_crystal_line(&line),
      _ => LineType::OtherLine
    }
  
  }
  
  
  
}
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<fs::File>>>
where P: AsRef<Path>, {
    let file = fs::File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn parse_config(args: &[String]) -> Result<String,&'static str>{
  if args.len() < 2{
    return Err("Please provide an input file.");
  }
  Ok(args[1].clone())
}
