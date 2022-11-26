//extern crate matrix_oxide;
//use matrix_oxide as mox;
//use num_complex::Complex;


// Numbers
//pub const I: Complex<f64> = Complex::<f64>{re: 0.0, im: 1.0};
//pub const ONE: Complex<f64> = Complex::<f64>{re: 1.0, im: 0.0};
//pub const I: mox::Z64 = mox::I;
//pub const ONE: mox::Z64 = mox::Z64{re: 1.0, im: 0.0};

// 2021: https://introcs.cs.princeton.edu/java/data/pi-10million.txt 
pub const PI: f64 = 3.141592653589793;

// Units
pub const METER: f64 = 1.0;
pub const SECOND: f64 = 1.0;
pub const JOULE: f64 = 1.0;
pub const TESLA: f64 = 1.0;
pub const KELVIN: f64 = 1.0;
pub const KILOGRAM: f64 = JOULE/METER/METER*SECOND*SECOND;
pub const HERTZ: f64 = 1.0/SECOND;
pub const ANGSTROM: f64 = METER*1e-10;
pub const AMPERE: f64 = KILOGRAM/SECOND/SECOND/TESLA;
pub const NEWTON: f64 = KILOGRAM*METER/SECOND/SECOND;
pub const COULOMB: f64 = AMPERE*SECOND; 
pub const VOLT: f64 = JOULE/COULOMB;

// 2022 https://physics.nist.gov/cgi-bin/cuu/Value?h
pub const PLANCK: f64 = 6.62607015e-34*JOULE/HERTZ;

// 2021: https://physics.nist.gov/cgi-bin/cuu/Value?hbar 
pub const HBAR: f64 = 1.054571817e-34 * JOULE * SECOND;

// 2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?mub
pub const MUB:f64 = 9.2740100783e-24 * JOULE/TESLA;

// 2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?mun
pub const MUN: f64 = 5.0507837461e-27 * JOULE/TESLA;

//2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?mu0
pub const MU0: f64 = 1.25663706212e-6 * NEWTON/AMPERE/AMPERE;

//2022 https://physics.nist.gov/cgi-bin/cuu/Value?k
pub const BOLTZMANN: f64 = 1.380649e-23 * JOULE/KELVIN;

// 2022 https://physics.nist.gov/cgi-bin/cuu/Value?na
pub const AVOGADRO: f64 = 6.02214076e23;

// 2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?gem 
pub const ELECTRON_G: f64 = -2.00231930436256;

// 2022 https://physics.nist.gov/cgi-bin/cuu/Value?e
pub const ELEMENTARY_CHARGE: f64 = 1.602176634e-19 * COULOMB;

pub const ELECTRON_VOLT: f64 = ELEMENTARY_CHARGE * VOLT;


// Conversions
pub const JOULES_TO_HERTZ: f64 = 1.0/PLANCK * HERTZ/JOULE;
pub const RAD_PER_S_TO_HZ: f64 = 0.5/PI;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum Element{
   Electron,                                                              
   Hydrogen,                                                              
   Deuterium,                                                             
   Helium,                                                                
   Lithium,                                                               
   Beryllium,                                                             
   Boron,                                                                 
   Carbon,                                                                
   Nitrogen,                                                              
   Oxygen,                                                               
   Fluorine,                                                             
   Neon,                                                                 
   Sodium,                                                               
   Magnesium,                                                            
   Aluminium,                                                            
   Silicon,                                                              
   Phosphorus,                                                           
   Sulfur,                                                               
   Chlorine,                                                             
   Argon,
}

impl Element{
  pub fn from(element: &str) -> Option<Element> {

    match element {
       "e" => return Some(Element::Electron),
       "H" | "D" | "T" => return Some(Element::Hydrogen),
      "He" => return Some(Element::Helium),
      "Li" => return Some(Element::Lithium),
      "Be" => return Some(Element::Beryllium),
       "B" => return Some(Element::Boron),
       "C" => return Some(Element::Carbon),
       "N" => return Some(Element::Nitrogen),
       "O" => return Some(Element::Oxygen),
       "F" => return Some(Element::Fluorine),
      "Ne" => return Some(Element::Neon),
      "Na" => return Some(Element::Sodium),
      "Al" => return Some(Element::Aluminium),
      "Si" => return Some(Element::Silicon),
       "P" => return Some(Element::Phosphorus),
       "S" => return Some(Element::Sulfur),
      "Cl" => return Some(Element::Chlorine),
      "Ar" => return Some(Element::Argon),
      _ => return None,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(PartialEq, Debug, Clone, Copy)]
pub enum Isotope{
  Electron,                                                                  
  Hydrogen1,                                                                
  Hydrogen2,                                                                
  Hydrogen3,                                                                
  Helium3,                                                                  
  Helium4,                                                                  
  Lithium6,                                                                 
  Lithium7,                                                                 
  Beryllium9,                                                               
  Boron10,                                                                  
  Boron11,                                                                  
  Carbon12,                                                                 
  Carbon13,                                                                 
  Nitrogen14,                                                               
  Nitrogen15,                                                               
  Oxygen16,                                                                 
  Oxygen17,                                                                 
  Fluorine19,                                                               
  Neon20,                                                                   
  Neon21,                                                                   
  Sodium22,                                                                 
  Sodium23,                                                                 
  Magnesium24,                                                              
  Magnesium25,                                                              
  Aluminium27,                                                              
  Silicon28,                                                                
  Silicon29,                                                                
  Phosphorus31,                                                             
  Sulfur32,                                                                 
  Sulfur33,                                                                 
  Chlorine35,                                                               
  Chlorine36,                                                               
  Chlorine37,                                                               
  Argon39,                                                                  
  Argon40,  
}


impl Isotope{

  pub fn g_value(&self) -> f64{

    match self{
      // https://www.physics.nist.gov/cgi-bin/cuu/Value?gem
      Isotope::Electron => return -2.00231930436256, 
      Isotope::Hydrogen1 => return 5.5869,
      Isotope::Hydrogen2 => return 0.857438,
      Isotope::Hydrogen3 => return 5.95799,
      Isotope::Helium3 => return -4.255,
      Isotope::Helium4 => return 0.0,
      Isotope::Lithium6 => return 0.822047,
      Isotope::Lithium7 => return 2.17095,
      Isotope::Beryllium9 => return -0.78495,
      Isotope::Boron10 => return 0.600215,
      Isotope::Boron11 => return 1.79243,
      Isotope::Carbon12 => return 0.0,
      Isotope::Carbon13 => return 1.40482,
      Isotope::Nitrogen14 => return 0.403761,
      Isotope::Nitrogen15 => return -0.566378,
      Isotope::Oxygen16 => return 0.0,
      Isotope::Oxygen17 => return -0.757516,
      Isotope::Fluorine19 => return 5.25774,
      Isotope::Neon20 => return 0.0,
      Isotope::Neon21 => return -0.441198,
      Isotope::Sodium22 => return 0.582,
      Isotope::Sodium23 => return 1.47835,
      Isotope::Magnesium24 => return 0.0,
      Isotope::Magnesium25 => return -0.34218,
      Isotope::Aluminium27 => return 1.4566,
      Isotope::Silicon28 => return 0.0,
      Isotope::Silicon29 => return -1.11058,
      Isotope::Phosphorus31 => return 2.2632,
      Isotope::Sulfur32 => return 0.0,
      Isotope::Sulfur33 => return 0.429214,
      Isotope::Chlorine35 => return 0.547916,
      Isotope::Chlorine36 => return 0.642735,
      Isotope::Chlorine37 => return 0.456082,
      Isotope::Argon39 => return -0.4537,
      Isotope::Argon40 => return 0.0, 
    }
  }
  //----------------------------------------------------------------------------

  pub fn spin_multiplicity(&self) -> usize {
    match self {
         Isotope::Electron => return 2,
         Isotope::Hydrogen1 => return 2,
         Isotope::Hydrogen2 => return 3,
         Isotope::Hydrogen3 => return 2,
         Isotope::Helium3 => return 2,
         Isotope::Helium4 => return 1,
         Isotope::Lithium6 => return 3,
         Isotope::Lithium7 => return 4,
         Isotope::Beryllium9 => return 4,
         Isotope::Boron10 => return 7,
         Isotope::Boron11 => return 4,
         Isotope::Carbon12 => return 1,
         Isotope::Carbon13 => return 2,
         Isotope::Nitrogen14 => return 3,
         Isotope::Nitrogen15 => return 2,
         Isotope::Oxygen16 => return 1,
         Isotope::Oxygen17 => return 6,
         Isotope::Fluorine19 => return 2,
         Isotope::Neon20 => return 1,
         Isotope::Neon21 => return 4,
         Isotope::Sodium22 => return 7,
         Isotope::Sodium23 => return 4,
         Isotope::Magnesium24 => return 1,
         Isotope::Magnesium25 => return 6,
         Isotope::Aluminium27 => return 6,
         Isotope::Silicon28 => return 1,
         Isotope::Silicon29 => return 2,
         Isotope::Phosphorus31 => return 2,
         Isotope::Sulfur32 => return 1,
         Isotope::Sulfur33 => return 4,
         Isotope::Chlorine35 => return 5,
         Isotope::Chlorine36 => return 5,
         Isotope::Chlorine37 => return 4,
         Isotope::Argon39 => return 8,
         Isotope::Argon40 => return 1,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
