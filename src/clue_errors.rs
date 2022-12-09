use std::fmt;

#[derive(PartialEq,Debug,Clone)]
pub enum CluEError{
  CannotConvertToFloat(usize,String),
  InvalidConfigFile(String),
  //InvalidToken(usize,String),
  NoCentralSpinCoor,
  NoRelationalOperators(usize),
  NotAnOperator(usize,String),
  UnmatchedBlockComment(usize),
  TooManyRelationalOperators(usize),
}

impl fmt::Display for CluEError{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{
      CluEError::CannotConvertToFloat(line_number, token) => write!(f,
          "{}: cannot convert \"{}\" to type float", line_number,token),

      CluEError::InvalidConfigFile(filename) => write!(f,
          "cannot not read config file \"{}\"", filename),

      //CluEError::InvalidToken(line_number,err_token) => write!(f,
      //    "{}:  invalid token \"{}\"",line_number, err_token),

      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),
      
      CluEError::NoRelationalOperators(line_number) => write!(f,
          "{}: no relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::NotAnOperator(line_number, token) => write!(f,
          "{}: cannot interpret \"{}\" as an operator", line_number,token),

      CluEError::UnmatchedBlockComment(line_number) => write!(f,
          "{}: unmatched \"*/\"", line_number),

      CluEError::TooManyRelationalOperators(line_number) => write!(f,
          "{}: too many relational operators (=, <, >, in, ...), are present", 
          line_number),
    }
  }
}
