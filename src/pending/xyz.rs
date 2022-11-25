
use std::fs;
use std::fs::File;
use std::io::{Error, Write};

pub struct XYZ{
  number: usize,
  elements: Vec<String>,
  x: Vec<f64>,
  y: Vec<f64>,
  z: Vec<f64>,
}

enum LineType{
 Number(usize),
 Element,
 Other,
}

impl XYZ {

pub fn print(&self) {
  println!("xyz = {}","{");
  println!("  n : {},", self.number);
  println!("  elements : {:?},", self.elements);
  println!("  x : {:?},", self.x);
  println!("  y : {:?},", self.y);
  println!("  z : {:?},", self.z);
  println!("{};","}");
}

pub fn write(&self, filename: &str) -> Result<(), Error> {
  let mut output = File::create(filename)?;
  let line = format!("{}\n",self.number);
  write!(output, "{}\n", line);

  for ii in 0..self.number {
    let line = format!("{}       {}       {}       {}\n",
        self.elements[ii], self.x[ii], self.y[ii], self.z[ii]);
    write!(output, "{}\n", line);
  }
  Ok(())
}

pub fn read(filename: &str) -> XYZ{
  println!("Reading {}.",filename);

  let lines = fs::read_to_string(filename)
    .expect("Could not read file.");
  let lines = lines.lines();


  let mut number = 0;
  let mut elements = Vec::<String>::new();
  let mut x = Vec::<f64>::new();
  let mut y = Vec::<f64>::new();
  let mut z = Vec::<f64>::new();

  for line in lines {
    if let LineType::Number(nn) = read_line_type(&line){
      number = nn;
    };
    let words = line.split_whitespace();
    let mut counter = 3;
    for word in words {
      let line_type = read_line_type(&word);
      match line_type {
        LineType::Element  => {
        counter = 0;
        elements.push(word.to_string());
        continue;
      },
        _ => (),
      }

      match counter {
        0 => {x.push(word.trim().parse().expect("No x."));
          counter += 1;},
        1 => {y.push(word.trim().parse().expect("No y."));
          counter += 1;},
        2 => {z.push(word.trim().parse().expect("No z."));
          counter += 1;},
        _ => {break;}
      } 
 
    }

  }
  XYZ {
    number,
    elements,
    x,
    y,
    z,
  }
}

}


fn read_line_type(word: &str) -> LineType {
  match word {
 "H"|"D"|"T"                                                                   |"He"|
"Li"|"Be"                                             | "B"| "C"| "N"| "O"| "F"|"Ne"|
"Na"|"Mg"                                             |"Al"|"Si"| "P"|"Se"|"Cl"|"Ar"|
 "K"|"Ca"|"Sc"|"Ti"| "V"|"Cr"|"Mn"|"Fe"|"Co"|"Cu"|"Zn"|"Ga"|"Ge"|"As"|"Se"|"Br"|"Kr"|  
"Rb"|"Sr"|"Zr"|"Nb"|"Mo"|"Tc"|"Ru"|"Rh"|"Pd"|"Ag"|"Cd"|"In"|"Sn"|"Sb"|"Te"| "I"|"Xe"
    => {
        LineType::Element
    },
      _ => {
        if let Ok(n) = word.trim().parse(){
          LineType::Number(n)
        } else {
          LineType::Other
        }
      },
  }
}

