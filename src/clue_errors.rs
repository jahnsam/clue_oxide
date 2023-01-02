use std::fmt;

#[derive(PartialEq,Debug,Clone)]
pub enum CluEError{
  CannotAddTokens,
  CannotCombineTokens(usize),
  CannotDivTokens,
  CannotMulTokens,
  CannotPowTokens,
  CannotSubTokens,
  CannotConvertToFloat(usize,String),
  CannotConvertToVector(usize),
  CannotParseElement(String),
  CannotParseLine(String),
  ConfigModeNotRecognized(String),
  EmptyVector(usize),
  IndexOutOfBounds(usize,usize,usize),
  InvalidConfigFile(String),
  InvalidToken(usize,String),
  ModeAttributeWrongBrackets,
  ModeAttributeWrongOption(String),
  ModeAttributeWrongSharp,
  NoCentralSpinCoor,
  NoClustersOfSize(usize),
  NoRelationalOperators(usize),
  NotAnOperator(usize,String),
  OptionAlreadySet(usize,String),
  UnmatchedBlockComment(usize),
  UnmatchedDelimiter(usize),
  TooManyRelationalOperators(usize),
  WrongVectorLength(usize,usize,usize)
}

impl fmt::Display for CluEError{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{

      CluEError::CannotAddTokens => write!(f,
          "cannot add tokens meaningfully"),

      CluEError::CannotCombineTokens(line_number) => write!(f,
          "{}: cannot combine tokens meaningfully", line_number),

      CluEError::CannotDivTokens => write!(f,
          "cannot divide tokens meaningfully"),

      CluEError::CannotMulTokens => write!(f,
          "cannot multiply tokens meaningfully"),

      CluEError::CannotParseElement(element) => write!(f,
          "cannot parse \"{}\" as an element", element),

      CluEError::CannotParseLine(line) => write!(f,
          "cannot parse line \"{}\"", line),

      CluEError::CannotSubTokens => write!(f,
          "cannot subtract tokens meaningfully"),

      CluEError::CannotPowTokens => write!(f,
          "cannot do token^token meaningfully"),

      CluEError::CannotConvertToFloat(line_number, token) => write!(f,
          "{}: cannot convert \"{}\" to type float", line_number,token),

      CluEError::CannotConvertToVector(line_number) => write!(f,
          "{}: cannot find vector", line_number),

      CluEError::ConfigModeNotRecognized(mode) => write!(f,
          "#[{}] is not recognized",mode),

      CluEError::EmptyVector(line_number) => write!(f,
          "{}: supplied vector is emptry", line_number),

      CluEError::IndexOutOfBounds(line_number,idx, len) => write!(f,
          "{}: cannot access element {} from array of length {}", 
          line_number, idx, len),

      CluEError::InvalidConfigFile(filename) => write!(f,
          "cannot not read config file \"{}\"", filename),

      CluEError::InvalidToken(line_number,err_token) => write!(f,
          "{}: invalid token \"{}\"",line_number, err_token),

      CluEError::ModeAttributeWrongBrackets => write!(f,
          "modes details should with square brackets"),

      CluEError::ModeAttributeWrongOption(mode)=> write!(f,
          "#[{}...] contains an invalid option",mode),

      CluEError::ModeAttributeWrongSharp => write!(f,
          "modes are specified with a single '#' at the start"),

      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),
      
      CluEError::NoClustersOfSize(size) => write!(f,
          "cannot find any clusters of size {}", size),

      CluEError::NoRelationalOperators(line_number) => write!(f,
          "{}: no relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::NotAnOperator(line_number, token) => write!(f,
          "{}: cannot interpret \"{}\" as an operator", line_number,token),

      CluEError::OptionAlreadySet(line_number,err_token) => write!(f,
          "{}: \"{}\" has already been set",line_number, err_token),

      CluEError::UnmatchedBlockComment(line_number) => write!(f,
          "{}: unmatched \"*/\"", line_number),

      CluEError::UnmatchedDelimiter(line_number) => write!(f,
          "{}: unmatched delimiter", line_number),

      CluEError::TooManyRelationalOperators(line_number) => write!(f,
          "{}: too many relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::WrongVectorLength(line_number, expected, actual) => write!(f,
          "{}: expected vector of length {}, but recieved a length of {}", 
          line_number, expected,actual),
    }
  }
}
