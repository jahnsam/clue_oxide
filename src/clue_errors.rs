use std::fmt;

#[derive(Debug,Clone)]
pub enum CluEError{
  InvalidConfigFile(String),
  //InvalidToken(usize,String),
  NoCentralSpinCoor,
  UnmatchedBlockComment(usize),
}

impl fmt::Display for CluEError{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{
      CluEError::InvalidConfigFile(filename) => write!(f,
          "cannot not read config file \"{}\"", filename),

      //CluEError::InvalidToken(line_number,err_token) => write!(f,
      //    "{}:  invalid token \"{}\"",line_number, err_token),

      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),

      CluEError::UnmatchedBlockComment(line_number) => write!(f,
          "{}: unmatched \"*/\"", line_number),
    }
  }
}
