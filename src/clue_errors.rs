use std::fmt;

#[derive(PartialEq,Debug,Clone)]
pub enum CluEError{
  AllSignalsNotSameLength(String),
  AtomDoesNotSpecifyElement(usize),
  AllVectorsNotSameLength(String),
  CannotAddTokens,
  CannontAugmentFilter(usize,String),
  CannotCombineTokens(usize),
  CannotConvertToFloat(usize,String),
  CannotConvertSerialToIndex(u32),
  CannotConvertToVector(usize),
  CannotDiagonalizeHamiltonian(String),
  CannotDivTokens,
  CannotFindCellID(usize),
  CannotFindSpinOp(String),
  CannotInferEigenvalues(usize),
  CannotMulTokens,
  CannotOpenFile(String),
  CannotParseElement(String),
  CannotParseLine(String),
  CannotParseIsotope(String),
  CannotParseRHS(usize),
  CannotParseSecondaryParticleFilter(String),
  CannotPowTokens,
  CannotSampleBinomialDistribution(usize,f64),
  CannotSetExchangeCoupling(usize),
  CannotSubTokens,
  CannotWriteFile(String),
  ClusterHasNoSignal(String),
  ConfigModeNotRecognized(String),
  EmptyVector(usize),
  ExpectedEquality(usize),
  ExpectedFloatRHS(usize),
  ExpectedIntRHS(usize),
  ExpectedNonNegativeIntRHS(usize),
  ExpectedVecOfNFloatsRHS(usize,usize),
  ExpectedNumber(usize),
  IncorrectNumberOfAxes(usize,usize),
  InorrectNumberOfCellOffsets(usize,usize),
  IndexOutOfBounds(usize,usize,usize),
  InvalidArgument(usize,String),
  InvalidAxes,
  InvalidClusterPartitionKey,
  InvalidConfigFile(String),
  InvalidPulseSequence(usize),
  InvalidToken(usize,String),
  LenghMismatchTimepointsIncrements(usize,usize),
  MissingFilter(String),
  MissingFilterArgument(usize,String),
  MissingFilterLabel(usize),
  MissingHeader(String),
  MissingProperties(String),
  MissingPropertiesLabel(usize),
  ModeAttributeWrongBrackets,
  ModeAttributeWrongOption(String),
  ModeAttributeWrongSharp,
  MultipleCosubstitutionGroups(usize),
  NoArgument(usize),
  NoClusterBatchSize,
  NoCentralSpin,
  NoCentralSpinCoor,
  NoCentralSpinIdentity,
  NoCentralSpinTransition,
  NoClusterMethod,
  NoClustersOfSize(usize),
  NoDensityMatrixMethod,
  NoHyperfineSpecifier(String,String),
  NoInputFile,
  NoLoadGeometry,
  NoMagneticField,
  NoModelIndex,
  NoMaxClusterSize,
  NoPulseSequence,
  NoQuadrupoleSpecifier(String,String),
  NoRadius,
  NoRelationalOperators(usize),
  NoRHS(usize),
  NoSpinOpForClusterSize(usize,usize),
  NoSpinOpWithMultiplicity(usize),
  NoStructureFile,
  NotAnOperator(usize,String),
  NotAProperSubset(String,String),
  NoTemperature,
  NoTensorValues,
  NoTimeAxis,
  NoTimeIncrements,
  NoTimepoints,
  OptionAlreadySet(usize,String),
  TensorNotSet(usize),
  TooManyRelationalOperators(usize),
  TooManyRHSArguments(usize),
  UnavailableSpinOp(usize,usize),
  UnassignedCosubstitutionGroup(usize),
  UnequalLengths(String,usize,String,usize),
  UnmatchedBlockComment(usize),
  UnmatchedDelimiter(usize),
  UnrecognizedOption(String),
  UnrecognizedVectorSpecifier(String),
  VectorSpecifierDoesNotSpecifyUniqueVector(String),
  WrongClusterSizeForAnalyticCCE(usize),
  WrongNumberOfAxes(usize,usize),
  WrongVectorLength(usize,usize,usize)
}

impl fmt::Display for CluEError{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{

      CluEError::AllSignalsNotSameLength(filename) => write!(f,
          "for \"{}\",signals must all have the same length", filename),

      CluEError::AtomDoesNotSpecifyElement(serial) => write!(f,
          "atom {} does not specify an element",serial),

      CluEError::AllVectorsNotSameLength(filename) => write!(f,
          "for \"{}\",signals must all have the same length", filename),

      CluEError::CannotAddTokens => write!(f,
          "cannot add tokens meaningfully"),

      CluEError::CannontAugmentFilter(index,secondary_filter) => write!(f,
          "cannot use secondary filter \"{}\" with particle {}",
         secondary_filter,index),

      CluEError::CannotCombineTokens(line_number) => write!(f,
          "line {}, cannot combine tokens meaningfully", line_number),

      CluEError::CannotConvertSerialToIndex(serial) => write!(f,
          "cannot convert serial id, {}, to an index",serial),

      CluEError::CannotDiagonalizeHamiltonian(matrix) => write!(f,
          "cannot diagonalize Hamiltonian,\n {}",matrix),

      CluEError::CannotDivTokens => write!(f,
          "cannot divide tokens meaningfully"),

      CluEError::CannotFindCellID(idx) => write!(f,
          "cannot determine cell id for particle {}",idx),

      CluEError::CannotFindSpinOp(sop) => write!(f,
          "cannot find spin operator \"{}\"",sop),

      CluEError::CannotInferEigenvalues(line_number) => write!(f,
          "line {}, cannot infer eigenvalues from input",line_number),
      
      CluEError::CannotMulTokens => write!(f,
          "cannot multiply tokens meaningfully"),

      CluEError::CannotOpenFile(file) => write!(f,
          "cannot open \"{}\"", file),

      CluEError::CannotParseElement(element) => write!(f,
          "cannot parse \"{}\" as an element", element),

      CluEError::CannotParseLine(line) => write!(f,
          "cannot parse line \"{}\"", line),

      CluEError::CannotParseIsotope(isotope) => write!(f,
          "cannot parse \"{}\" as an isotope", isotope),

      CluEError::CannotParseRHS(line_number) => write!(f,
          "line {}, cannot parse line righ hand side", line_number),

      CluEError::CannotParseSecondaryParticleFilter(filter) => write!(f,
          "cannot parse secondary particle filter \"{}\"", filter),

      CluEError::CannotSampleBinomialDistribution(n,p) => write!(f,
          "cannot sample from the binomial distribution B(n={},p={})",
          n,p),

      CluEError::CannotSetExchangeCoupling(line_number) => write!(f,
          "line {}, cannot set exchange coupling",line_number),

      CluEError::CannotSubTokens => write!(f,
          "cannot subtract tokens meaningfully"),

      CluEError::CannotPowTokens => write!(f,
          "cannot do token^token meaningfully"),

      CluEError::CannotConvertToFloat(line_number, token) => write!(f,
          "line {}, cannot convert \"{}\" to type float", line_number,token),

      CluEError::CannotConvertToVector(line_number) => write!(f,
          "line {}, cannot find vector", line_number),

      CluEError::ClusterHasNoSignal(cluster) => write!(f,
          "expected cluster {} to have a signal, but found none", cluster),

      CluEError::CannotWriteFile(file) => write!(f,
          "cannot write to \"{}\"", file),

      CluEError::ConfigModeNotRecognized(mode) => write!(f,
          "#[{}] is not recognized",mode),

      CluEError::EmptyVector(line_number) => write!(f,
          "line {}, supplied vector is emptry", line_number),

      CluEError::ExpectedEquality(line_number) => write!(f,
          "line {}, expected an equaliy",line_number),

      CluEError::ExpectedFloatRHS(line_number) => write!(f,
          "line {}, expected a float on the right hand side",line_number),

      CluEError::ExpectedIntRHS(line_number) => write!(f,
          "line {}, expected an integer on the right hand side",line_number),

      CluEError::ExpectedNonNegativeIntRHS(line_number) => write!(f,
          "line {}, expected a non-negative integer on the right hand side",
          line_number),

      CluEError::ExpectedNumber(line_number) => write!(f,
          "line {}, expected a number",line_number),

      CluEError::ExpectedVecOfNFloatsRHS(line_number,n) => write!(f,
          "line {}, expected a vector of {} floats on the right hand side",
          line_number,n),

      CluEError::IncorrectNumberOfAxes(n,n_ref)=> write!(f,
          "expected {} axes, but {} were provided",n_ref, n),

      CluEError::InorrectNumberOfCellOffsets(n,n_ref)=> write!(f,
          "expected {} cell offsets, but {} were provided",n_ref, n),

      CluEError::IndexOutOfBounds(line_number,idx, len) => write!(f,
          "line {}, cannot access element {} from array of length {}", 
          line_number, idx, len),

      CluEError::InvalidAxes => write!(f,
          "invalid axes"),

      CluEError::InvalidArgument(line_number,expected_arg) => write!(f,
          "line {}, argument should be a(n) {}",line_number, expected_arg),

      CluEError::InvalidClusterPartitionKey => write!(f,
          "invalid cluster partition key"),

      CluEError::InvalidConfigFile(filename) => write!(f,
          "cannot not read config file \"{}\"", filename),

      CluEError::InvalidPulseSequence(line_number) => write!(f,
          "line {}, invalid pulse sequence",line_number),

      CluEError::InvalidToken(line_number,err_token) => write!(f,
          "line {}, invalid token \"{}\"",line_number, err_token),

      CluEError::LenghMismatchTimepointsIncrements(n_dts,dts) => write!(f,
          "there are {} timepoint specifications, but {} time increments",
          n_dts,dts),

      CluEError::MissingFilter(label) => write!(f,
          "no filter with label \"{}\"",label),

      CluEError::MissingFilterArgument(line_number,fn_name) => write!(f,
          "line {}, missing filter argument in \"{}()\"",line_number,fn_name),

      CluEError::MissingFilterLabel(line_number) => write!(f,
          "line {}, missing label in filter",line_number),

      CluEError::MissingHeader(filename) => write!(f,
          "in \"{}\", every entry must have a header",filename),

      CluEError::MissingProperties(label) => write!(f,
          "no properties with label \"{}\"",label),

      CluEError::MissingPropertiesLabel(line_number) => write!(f,
          "line {}, missing label in properties",line_number),

      CluEError::ModeAttributeWrongBrackets => write!(f,
          "modes details should with square brackets"),

      CluEError::ModeAttributeWrongOption(mode)=> write!(f,
          "#[{}...] contains an invalid option",mode),

      CluEError::MultipleCosubstitutionGroups(index)=> write!(f,
          "particle {} is assigned to more than one cosubstitution group",
          index),

      CluEError::ModeAttributeWrongSharp => write!(f,
          "modes are specified with a single '#' at the start"),

      CluEError::NoArgument(line_number) => write!(f,
          "line {}, expected a function argument",line_number),

      CluEError::NoClusterBatchSize => write!(f,
          "batch size for clusters not specified"),

      CluEError::NoCentralSpin => write!(f,
          "no detected spin defined"),
      
      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),
      
      CluEError::NoCentralSpinIdentity => write!(f,
          "the detected spin's identity was not defined"),
      
      CluEError::NoCentralSpinTransition => write!(f,
          "the detected spin transition was not defined"),
      
      CluEError::NoClusterMethod => write!(f,
          "no cluster method specified"),
      
      CluEError::NoClustersOfSize(size) => write!(f,
          "cannot find any clusters of size {}", size),

      CluEError::NoDensityMatrixMethod=> write!(f,
          "no density matrix method specified"),
      
      CluEError::NoInputFile => write!(f,
          "no input file"),
      
      CluEError::NoHyperfineSpecifier(label,isotope) => write!(f,
          "no hyperfine specifier found for {} {}",label,isotope),
      
      CluEError::NoLoadGeometry => write!(f,
          "geometry for loading in the system was not defined"),
      
      CluEError::NoMagneticField => write!(f,
          "please specify the applied magnetic field"),

      CluEError::NoModelIndex => write!(f,
          "PDB model not selected"),

      CluEError::NoMaxClusterSize => write!(f,
          "maximum cluster size not set"),

      CluEError::NoPulseSequence => write!(f,
          "no pulse sequence is defined"),
      
      CluEError::NoQuadrupoleSpecifier(label,isotope) => write!(f,
          "no quadrupole specifier found for {} {}",label,isotope),
      
      CluEError::NoRadius => write!(f,
          "system radius not set"),

      CluEError::NoRelationalOperators(line_number) => write!(f,
          "line {}, no relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::NotAnOperator(line_number, token) => write!(f,
          "line {}, cannot interpret \"{}\" as an operator", line_number,token),

      CluEError::NotAProperSubset(subcluster, cluster) => write!(f,
          "{} is not a proper subset of {}", subcluster, cluster),

      CluEError::NoTemperature => write!(f,
          "no temperature specified"),
      
      CluEError::NoTensorValues => write!(f,
          "no values specified for tensor"),
      
      CluEError::NoTimeAxis => write!(f,
          "the time-axis has not been built"),
      
      CluEError::NoTimeIncrements => write!(f,
          "no time increments defined"),
      
      CluEError::NoTimepoints => write!(f,
          "please specify how many timepoints there are for each increment"),
      
      CluEError::NoRHS(line_number) => write!(f,
          "line {}, cannot read right hand side",line_number),

      CluEError::NoSpinOpForClusterSize(cluster_size,max_size) => write!(f,
          "no spin operators for clusters of size {} are built,\
          only sizes upto {}",cluster_size,max_size),

      CluEError::NoSpinOpWithMultiplicity(spin_multiplicity) => write!(f,
          "no spin-{} operators are built", 
          (*spin_multiplicity as f64 - 1.0)/2.0),

      CluEError::NoStructureFile => write!(f,
          "no structure file defined"),

      CluEError::OptionAlreadySet(line_number,err_token) => write!(f,
          "line {}, \"{}\" has already been set",line_number, err_token),

      CluEError::UnmatchedBlockComment(line_number) => write!(f,
          "line {}, unmatched \"*/\"", line_number),

      CluEError::UnmatchedDelimiter(line_number) => write!(f,
          "line {}, unmatched delimiter", line_number),

      CluEError::TensorNotSet(index) => write!(f,
          "no tensor set for particle {}",index),

      CluEError::TooManyRelationalOperators(line_number) => write!(f,
          "line {}, too many relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::TooManyRHSArguments(line_number) => write!(f,
          "line {}, too many arguments on the right hand side", 
          line_number),

      CluEError::UnassignedCosubstitutionGroup(index)=> write!(f,
          "particle {} cannot be assigned to a cosubstitution group",
          index),

      CluEError::UnavailableSpinOp(op_pos,n_ops) => write!(f,
          "spin operator index, {}, exceeds vector length of {}",
          op_pos,n_ops),

      CluEError::UnequalLengths(vec0,len0,vec1,len1) => write!(f,
          "unequal lengths, {} has length {}, but {} has length {}",
          vec0,len0,vec1,len1),

      CluEError::UnrecognizedOption(option) => write!(f,
          "unrecognized option \"{}\"",option),

      CluEError::UnrecognizedVectorSpecifier(specifier) => write!(f,
          "unrecognized vector specifier \"{}\"",specifier),

      CluEError::VectorSpecifierDoesNotSpecifyUniqueVector(specifier) 
        => write!(f,
          "vector specifier does not specify unique vector \"{}\"",specifier),

      CluEError::WrongClusterSizeForAnalyticCCE(given_size) => write!(f,
          "analytic 2-CCE cannot work with clusters of size {}",given_size),

      CluEError::WrongNumberOfAxes(num_axes, expected_num) => write!(f,
          "{} axes were provided, but {} are expected",num_axes, expected_num),

      CluEError::WrongVectorLength(line_number, expected, actual) => write!(f,
          "line {}, expected vector of length {}, but recieved a length of {}", 
          line_number, expected,actual),
    }
  }
}
