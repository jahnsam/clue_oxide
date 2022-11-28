use std::fmt;

#[derive(Debug,Clone)]
pub enum CluEError{
  NoCentralSpinCoor,
}

impl fmt::Display for CluEError{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{
      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),
    }
  }
}
