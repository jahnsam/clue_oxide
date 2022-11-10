
#[derive(PartialEq)]
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
       "H" => return Some(Element::Hydrogen),
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
