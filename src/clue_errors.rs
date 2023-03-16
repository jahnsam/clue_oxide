use std::fmt;

#[derive(PartialEq,Debug,Clone)]
pub enum CluEError{
  AtomDoesNotSpecifyElement(usize),
  CannotAddTokens,
  CannotCombineTokens(usize),
  CannotConvertSerialToIndex(u32),
  CannotDivTokens,
  CannotMulTokens,
  CannotPowTokens,
  CannotSampleBinomialDistribution(usize,f64),
  CannotSubTokens,
  CannontAugmentFilter(usize,String),
  CannotConvertToFloat(usize,String),
  CannotConvertToVector(usize),
  CannotOpenFile(String),
  CannotParseElement(String),
  CannotParseLine(String),
  ConfigModeNotRecognized(String),
  EmptyVector(usize),
  ExpectedEquality(usize),
  ExpectedFloatRHS(usize),
  ExpectedVecOfNFloatsRHS(usize,usize),
  IndexOutOfBounds(usize,usize,usize),
  InvalidConfigFile(String),
  InvalidToken(usize,String),
  ModeAttributeWrongBrackets,
  ModeAttributeWrongOption(String),
  ModeAttributeWrongSharp,
  MultipleCosubstitutionGroups(usize),
  NoCentralSpinCoor,
  NoClustersOfSize(usize),
  NoRadius,
  NoRelationalOperators(usize),
  NoRHS(usize),
  NotAnOperator(usize,String),
  OptionAlreadySet(usize,String),
  UnmatchedBlockComment(usize),
  UnmatchedDelimiter(usize),
  TooManyRelationalOperators(usize),
  UnassignedCosubstitutionGroup(usize),
  UnrecognizedOption(String),
  WrongVectorLength(usize,usize,usize)
}

impl fmt::Display for CluEError{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{

      CluEError::AtomDoesNotSpecifyElement(serial) => write!(f,
          "atom {} does not specify an element",serial),

      CluEError::CannotAddTokens => write!(f,
          "cannot add tokens meaningfully"),

      CluEError::CannontAugmentFilter(index,secondary_filter) => write!(f,
          "cannot use secondary filter \"{}\" with particle {}",
         secondary_filter,index),

      CluEError::CannotCombineTokens(line_number) => write!(f,
          "{}: cannot combine tokens meaningfully", line_number),

      CluEError::CannotConvertSerialToIndex(serial) => write!(f,
          "cannot convert serial id, {}, to an index",serial),

      CluEError::CannotDivTokens => write!(f,
          "cannot divide tokens meaningfully"),

      CluEError::CannotMulTokens => write!(f,
          "cannot multiply tokens meaningfully"),

      CluEError::CannotOpenFile(file) => write!(f,
          "cannot open \"{}\"", file),

      CluEError::CannotParseElement(element) => write!(f,
          "cannot parse \"{}\" as an element", element),

      CluEError::CannotParseLine(line) => write!(f,
          "cannot parse line \"{}\"", line),

      CluEError::CannotSampleBinomialDistribution(n,p) => write!(f,
          "cannot sample from the binomial distribution B(n={},p={})",
          n,p),

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

      CluEError::ExpectedEquality(line_number) => write!(f,
          "{}: expected an equaliy",line_number),

      CluEError::ExpectedFloatRHS(line_number) => write!(f,
          "{}: expected a float on the right hand side",line_number),

      CluEError::ExpectedVecOfNFloatsRHS(line_number,n) => write!(f,
          "{}: expected a vector of {} floats on the right hand side",
          line_number,n),

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

      CluEError::MultipleCosubstitutionGroups(index)=> write!(f,
          "particle {} is assigned to more than one cosubstitution group",
          index),

      CluEError::ModeAttributeWrongSharp => write!(f,
          "modes are specified with a single '#' at the start"),

      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),
      
      CluEError::NoClustersOfSize(size) => write!(f,
          "cannot find any clusters of size {}", size),

      CluEError::NoRadius => write!(f,
          "system radius not set"),

      CluEError::NoRelationalOperators(line_number) => write!(f,
          "{}: no relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::NotAnOperator(line_number, token) => write!(f,
          "{}: cannot interpret \"{}\" as an operator", line_number,token),

      CluEError::NoRHS(line_number) => write!(f,
          "{}: cannot read right hand side",line_number),

      CluEError::OptionAlreadySet(line_number,err_token) => write!(f,
          "{}: \"{}\" has already been set",line_number, err_token),

      CluEError::UnmatchedBlockComment(line_number) => write!(f,
          "{}: unmatched \"*/\"", line_number),

      CluEError::UnmatchedDelimiter(line_number) => write!(f,
          "{}: unmatched delimiter", line_number),

      CluEError::TooManyRelationalOperators(line_number) => write!(f,
          "{}: too many relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::UnassignedCosubstitutionGroup(index)=> write!(f,
          "particle {} cannot be assigned to a cosubstitution group",
          index),

      CluEError::UnrecognizedOption(option) => write!(f,
          "unrecognized option \"{}\"",option),

      CluEError::WrongVectorLength(line_number, expected, actual) => write!(f,
          "{}: expected vector of length {}, but recieved a length of {}", 
          line_number, expected,actual),
    }
  }
}
