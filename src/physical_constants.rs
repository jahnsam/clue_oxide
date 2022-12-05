use strum_macros::EnumIter;
//extern crate matrix_oxide;
//use matrix_oxide as mox;
use num_complex::Complex;


// Numbers
pub const I: Complex<f64> = Complex::<f64>{re: 0.0, im: 1.0};
pub const ONE: Complex<f64> = Complex::<f64>{re: 1.0, im: 0.0};
//pub const I: mox::Z64 = mox::I;
//pub const ONE: mox::Z64 = mox::Z64{re: 1.0, im: 0.0};

// 2021: https://introcs.cs.princeton.edu/java/data/pi-10million.txt 
//pub const PI: f64 = 3.141592653589793;
pub const PI: f64 = std::f64::consts::PI;

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
/// Specifies a particle type, such as a chemical element or an electron spin.
#[derive(PartialEq, Debug, Clone, Copy, EnumIter)]
pub enum Element{
   Electron,                                                              
   Hydrogen,                                                              
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
       "e" => Some(Element::Electron),
       "H" | "D" | "T" => Some(Element::Hydrogen),
      "He" => Some(Element::Helium),
      "Li" => Some(Element::Lithium),
      "Be" => Some(Element::Beryllium),
       "B" => Some(Element::Boron),
       "C" => Some(Element::Carbon),
       "N" => Some(Element::Nitrogen),
       "O" => Some(Element::Oxygen),
       "F" => Some(Element::Fluorine),
      "Ne" => Some(Element::Neon),
      "Na" => Some(Element::Sodium),
      "Al" => Some(Element::Aluminium),
      "Si" => Some(Element::Silicon),
       "P" => Some(Element::Phosphorus),
       "S" => Some(Element::Sulfur),
      "Cl" => Some(Element::Chlorine),
      "Ar" => Some(Element::Argon),
      _ => None,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// Specifies the particle type and specific varient.  
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
      Isotope::Electron => -2.00231930436256, 
      Isotope::Hydrogen1 => 5.5869,
      Isotope::Hydrogen2 => 0.857438,
      Isotope::Hydrogen3 => 5.95799,
      Isotope::Helium3 => -4.255,
      Isotope::Helium4 => 0.0,
      Isotope::Lithium6 => 0.822047,
      Isotope::Lithium7 => 2.17095,
      Isotope::Beryllium9 => -0.78495,
      Isotope::Boron10 => 0.600215,
      Isotope::Boron11 => 1.79243,
      Isotope::Carbon12 => 0.0,
      Isotope::Carbon13 => 1.40482,
      Isotope::Nitrogen14 => 0.403761,
      Isotope::Nitrogen15 => -0.566378,
      Isotope::Oxygen16 => 0.0,
      Isotope::Oxygen17 => -0.757516,
      Isotope::Fluorine19 => 5.25774,
      Isotope::Neon20 => 0.0,
      Isotope::Neon21 => -0.441198,
      Isotope::Sodium22 => 0.582,
      Isotope::Sodium23 => 1.47835,
      Isotope::Magnesium24 => 0.0,
      Isotope::Magnesium25 => -0.34218,
      Isotope::Aluminium27 => 1.4566,
      Isotope::Silicon28 => 0.0,
      Isotope::Silicon29 => -1.11058,
      Isotope::Phosphorus31 => 2.2632,
      Isotope::Sulfur32 => 0.0,
      Isotope::Sulfur33 => 0.429214,
      Isotope::Chlorine35 => 0.547916,
      Isotope::Chlorine36 => 0.642735,
      Isotope::Chlorine37 => 0.456082,
      Isotope::Argon39 => -0.4537,
      Isotope::Argon40 => 0.0, 
    }
  }
  //----------------------------------------------------------------------------

  pub fn spin_multiplicity(&self) -> usize {
    match self {
         Isotope::Electron => 2,
         Isotope::Hydrogen1 => 2,
         Isotope::Hydrogen2 => 3,
         Isotope::Hydrogen3 => 2,
         Isotope::Helium3 => 2,
         Isotope::Helium4 => 1,
         Isotope::Lithium6 => 3,
         Isotope::Lithium7 => 4,
         Isotope::Beryllium9 => 4,
         Isotope::Boron10 => 7,
         Isotope::Boron11 => 4,
         Isotope::Carbon12 => 1,
         Isotope::Carbon13 => 2,
         Isotope::Nitrogen14 => 3,
         Isotope::Nitrogen15 => 2,
         Isotope::Oxygen16 => 1,
         Isotope::Oxygen17 => 6,
         Isotope::Fluorine19 => 2,
         Isotope::Neon20 => 1,
         Isotope::Neon21 => 4,
         Isotope::Sodium22 => 7,
         Isotope::Sodium23 => 4,
         Isotope::Magnesium24 => 1,
         Isotope::Magnesium25 => 6,
         Isotope::Aluminium27 => 6,
         Isotope::Silicon28 => 1,
         Isotope::Silicon29 => 2,
         Isotope::Phosphorus31 => 2,
         Isotope::Sulfur32 => 1,
         Isotope::Sulfur33 => 4,
         Isotope::Chlorine35 => 5,
         Isotope::Chlorine36 => 5,
         Isotope::Chlorine37 => 4,
         Isotope::Argon39 => 8,
         Isotope::Argon40 => 1,
    }
  }
  //----------------------------------------------------------------------------
  pub fn most_common_for(element: &Element) -> Isotope{

    match element{
      Element::Electron => Isotope::Electron,
      Element::Hydrogen => Isotope::Hydrogen1,
      //Element::Deuterium => Isotope::Hydrogen2,
      Element::Helium => Isotope::Helium4,
      Element::Lithium => Isotope::Lithium7,
      Element::Beryllium => Isotope::Beryllium9,
      Element::Boron => Isotope::Boron11,
      Element::Carbon => Isotope::Carbon12,
      Element::Nitrogen => Isotope::Nitrogen14,
      Element::Oxygen => Isotope::Oxygen16,
      Element::Fluorine => Isotope::Fluorine19,
      Element::Neon => Isotope::Neon20,
      Element::Sodium => Isotope::Sodium23,
      Element::Magnesium => Isotope::Magnesium24,
      Element::Aluminium => Isotope::Aluminium27,
      Element::Silicon => Isotope::Silicon28,
      Element::Phosphorus => Isotope::Phosphorus31,
      Element::Sulfur => Isotope::Sulfur32,
      Element::Chlorine => Isotope::Chlorine35,
      Element::Argon => Isotope::Argon40,
    }
  
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
