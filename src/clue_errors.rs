use std::fmt;

#[derive(PartialEq,Debug,Clone)]
pub enum CluEError{
  AllSignalsNotSameLength(String),
  AtomDoesNotSpecifyElement(usize),
  AllVectorsNotSameLength(String),
  BondsAreNotDefined,
  CannotAddPointToGrid(usize,usize),
  CannotAddTokens,
  CannontAugmentFilter(usize,String),
  CannotCombineTokens(usize),
  CannotConvertToFloat(usize,String),
  CannotConvertSerialToIndex(u32),
  CannotConvertToVector(usize),
  CannotCreateDir(String),
  CannotDiagonalizeHamiltonian(String),
  CannotDivTokens,
  CannotFindCellID(usize),
  CannotFindSpinOp(String),
  CannotFindRefIndexFromBathIndex(usize),
  CannotFindRefIndexFromNthActive(usize),
  CannotInferEigenvalues(usize),
  CannotMatchVertexToIndex(usize),
  CannotMulTokens,
  CannotOpenFile(String),
  CannotParseElement(String),
  CannotParseLine(String),
  CannotParseIsotope(String),
  CannotParseRHS(usize),
  CannotParseSecondaryParticleFilter(String),
  CannotPowTokens,
  CannotReadGrid(String),
  CannotSampleBinomialDistribution(usize,f64),
  CannotSetExchangeCoupling(usize),
  CannotSubTokens,
  CannotTakeTrace(String),
  CannotWriteFile(String),
  ClusterHasNoSignal(String),
  ConfigModeNotRecognized(String),
  EmptyVector(usize),
  ExpectedClusterSetWithNSizes(usize,usize),
  ExpectedEquality(usize),
  ExpectedBoolRHS(usize),
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
  InvalidGeometry(usize,String),
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
  NANTensorBathDipoleDipole(usize,String,usize,String),
  NANTensorBathZeeman(usize,String),
  NANTensorDetectedZeeman,
  NANTensorExchangeCoupling(usize,String,usize,String),
  NANTensorHyperfine(usize,String),
  NANTensorQuadrupole(usize,String),
  NoArgument(usize),
  NoCentralSpin,
  NoCentralSpinCoor,
  NoCentralSpinIdentity,
  NoCentralSpinTransition,
  NoClashDistancePBC,
  NoClusterBatchSize,
  NoClusterMethod,
  NoClustersOfSize(usize),
  NoDensityMatrixMethod,
  NoDetectedSpinIdentity,
  NoDetectedSpinMultiplicity,
  NoGMatrixSpecifier,
  NoGMatrixValues,
  NoHyperfineSpecifier(String,String),
  NoInputFile,
  NoLoadGeometry,
  NoMagneticField,
  NoModelIndex,
  NoMaxClusterSize,
  NoNeighborCutoffDistance,
  NoOrientationGrid,
  NoPulseSequence,
  NoQuadrupoleSpecifier(String,String),
  NoRadius,
  NoRelationalOperators(usize),
  NoRemovePartialMethyls,
  NoRHS(usize),
  NoSpinOpForClusterSize(usize,usize),
  NoSpinOpWithMultiplicity(usize),
  NoStructureFile,
  NotA3DVector(usize),
  NotALebedevGrid(usize),
  NotAnOperator(usize,String),
  NotAProperSubset(String,String),
  NoTemperature,
  NoTensorValues,
  NoTimeAxis,
  NoTimeIncrements,
  NoTimepoints,
  OptionAlreadySet(usize,String),
  ParticlesClash(usize,String,usize,String,f64,f64),
  SaveNameEmpty,
  SaveNameNotSet,
  SecondaryFilterRequiresAnIndex(String),
  SpinPropertiesNeedsALabel,
  SpinPropertiesNeedsAnIsotope(String),
  StructurePropertiesNeedsALabel,
  StructurePropertiesDoesNotNeedAnIsotope(String),
  TensorNotSet(usize),
  TooFewRHSArguments(usize),
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
  WrongOrientationGridDim(usize,usize,usize),
  WrongProbabilityDistributionDim(usize,usize,usize),
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

      CluEError:: BondsAreNotDefined => write!(f,
          "no chemical bonds are established"),

      CluEError::CannotAddPointToGrid(point_dim, grid_dim) => write!(f,
          "cannot add {}D point t0 {}D grid",point_dim, grid_dim),

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

      CluEError::CannotCreateDir(path) => write!(f,
          "cannot create directory \"{}\"",path),

      CluEError::CannotDivTokens => write!(f,
          "cannot divide tokens meaningfully"),

      CluEError::CannotFindCellID(idx) => write!(f,
          "cannot determine cell id for particle {}",idx),

      CluEError::CannotFindSpinOp(sop) => write!(f,
          "cannot find spin operator \"{}\"",sop),

      CluEError::CannotFindRefIndexFromBathIndex(bath_index) => write!(f,
          "cannot find reference index for bath index \"{}\"", bath_index),

      CluEError::CannotFindRefIndexFromNthActive(n) => write!(f,
          "cannot find reference index for the {}th active particle", n),

      CluEError::CannotInferEigenvalues(line_number) => write!(f,
          "line {}, cannot infer eigenvalues from input",line_number),
      
      CluEError::CannotMatchVertexToIndex(vertex) => write!(f,
          "cannot match vertex {} to an index",vertex),
      
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
          "line {}, cannot parse line right hand side", line_number),

      CluEError::CannotParseSecondaryParticleFilter(filter) => write!(f,
          "cannot parse secondary particle filter \"{}\"", filter),

      CluEError::CannotSampleBinomialDistribution(n,p) => write!(f,
          "cannot sample from the binomial distribution B(n={},p={})",
          n,p),

      CluEError::CannotSetExchangeCoupling(line_number) => write!(f,
          "line {}, cannot set exchange coupling",line_number),

      CluEError::CannotSubTokens => write!(f,
          "cannot subtract tokens meaningfully"),

      CluEError::CannotTakeTrace(matrix) => write!(f,
          "cannot take trace of \"{}\"",matrix),

      CluEError::CannotPowTokens => write!(f,
          "cannot do token^token meaningfully"),

      CluEError::CannotReadGrid(filename) => write!(f,
          "cannot read grid from \"{}\": \
the grid should be specified as a csv file with one column per dimension, \
followed by a column for the weights", filename),

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

      CluEError::ExpectedClusterSetWithNSizes(n_exp,n_act) => write!(f,
          "expected a cluster set with {} sizes, but got {} sizes",
          n_exp, n_act),

      CluEError::ExpectedEquality(line_number) => write!(f,
          "line {}, expected an equaliy",line_number),

      CluEError::ExpectedBoolRHS(line_number) => write!(f,
          "line {}, expected a bool on the right hand side",line_number),

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

      CluEError::InvalidGeometry(line_number,rhs) => write!(f,
          "line {}, argument \"{}\" is not a valid load_geometry",
          line_number, rhs),

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

      CluEError::NANTensorBathDipoleDipole(idx0,isotope0,idx1,isotope1) 
        => write!(f,
          "dipole-dipole tensor for particles {} {} and {} {} contains NANs",
          idx0,isotope0,idx1,isotope1),

      CluEError::NANTensorBathZeeman(idx,isotope) => write!(f,
          "zeeman tensor for particle {} {} contains NANs",idx,isotope),

      CluEError::NANTensorDetectedZeeman => write!(f,
          "Zeeman tensor for the detected particle contains NANs"),

      CluEError::NANTensorExchangeCoupling(idx0,isotope0,idx1,isotope1) 
        => write!(f,
          "exchange tensor for particles {} {} and {} {} contains NANs",
          idx0,isotope0,idx1,isotope1),

      CluEError::NANTensorHyperfine(idx,isotope) => write!(f,
          "hyperfine tensor for particle {} {} contains NANs",idx,isotope),

      CluEError::NANTensorQuadrupole(idx,isotope) => write!(f,
          "electric quadrupole tensor for particle {} {} contains NANs",
          idx,isotope),

      CluEError::NoArgument(line_number) => write!(f,
          "line {}, expected a function argument",line_number),

      CluEError::NoCentralSpin => write!(f,
          "no detected spin defined"),
      
      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),
      
      CluEError::NoCentralSpinIdentity => write!(f,
          "the detected spin's identity was not defined"),
      
      CluEError::NoCentralSpinTransition => write!(f,
          "the detected spin transition was not defined"),

      CluEError::NoClashDistancePBC => write!(f,
          "clash_distance_pbc is not defined"),

      CluEError::NoClusterBatchSize => write!(f,
          "batch size for clusters not specified"),

      CluEError::NoClusterMethod => write!(f,
          "no cluster method specified"),
      
      CluEError::NoClustersOfSize(size) => write!(f,
          "cannot find any clusters of size {}", size),
      
      CluEError::NoDensityMatrixMethod=> write!(f,
          "no density matrix method specified"),
      
      CluEError::NoDetectedSpinIdentity => write!(f,
          "detected_spin_identity is not set"),

      CluEError::NoDetectedSpinMultiplicity => write!(f,
          "detected_spin_identity is not set"),

      CluEError::NoGMatrixSpecifier => write!(f,
          "no g-matrix specifier"),
      
      CluEError::NoGMatrixValues => write!(f,
          "g-matrix values are not set"),
      
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

      CluEError::NoNeighborCutoffDistance => write!(f,
          "neighbor_cutoff_distance is not defined"),
      
      CluEError::NoOrientationGrid => write!(f,
          "orientation_grid is not defined"),
      
      CluEError::NoPulseSequence => write!(f,
          "no pulse sequence is defined"),
      
      CluEError::NoQuadrupoleSpecifier(label,isotope) => write!(f,
          "no quadrupole specifier found for {} {}",label,isotope),
      
      CluEError::NoRadius => write!(f,
          "system radius not set"),

      CluEError::NoRelationalOperators(line_number) => write!(f,
          "line {}, no relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::NoRemovePartialMethyls => write!(f,
          "remove_partial_methyls: bool is not set"),
      
      CluEError::NotA3DVector(dim) => write!(f,
          "vector is {}-dimensional, not 3-dimensional", 
          dim),

      CluEError::NotALebedevGrid(n_ori) => write!(f,
          "\"{}\" is not a valid Lebedev grid; please choose n_ori in \
{{6, 14, 26, 38, 50, 74, 86, 110, 146, 170, \
194,  230, 266, 302, 350, 434, 590, 770, 974, 1202, \
1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,5294, 5810}}",n_ori),

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

      CluEError::ParticlesClash(idx0,elmt0,idx1,elmt1,r,r_clash) => write!(f,
          "particle {} {} and particle {} {} are {} Å apart, closer than the
clash distance of {} Å",idx0,elmt0,idx1,elmt1,r,r_clash),

      CluEError::OptionAlreadySet(line_number,err_token) => write!(f,
          "line {}, \"{}\" has already been set",line_number, err_token),

      CluEError::UnmatchedBlockComment(line_number) => write!(f,
          "line {}, unmatched \"*/\"", line_number),

      CluEError::UnmatchedDelimiter(line_number) => write!(f,
          "line {}, unmatched delimiter", line_number),

      CluEError::SaveNameEmpty => write!(f,
          "save name is empty"),

      CluEError::SaveNameNotSet => write!(f,
          "save name not set"),

      CluEError::SecondaryFilterRequiresAnIndex(filter) => write!(f,
          "secondary particle filter, \"{}\", requires a particle index",
          filter),

      CluEError::SpinPropertiesNeedsALabel => write!(f,
          "spin_properties requires a label to be set: 
#[spin_properties(label = LABEL, isotope = ISOTOPE)]"),

      CluEError::SpinPropertiesNeedsAnIsotope(label) => write!(f,
          "spin_properties requires an isotope to be set: 
#[spin_properties(label = {}, isotope = ISOTOPE)]",label),

      CluEError::StructurePropertiesNeedsALabel => write!(f,
          "structure_properties requires a label to be set: 
#[structure_properties(label = LABEL)]"),

      CluEError::StructurePropertiesDoesNotNeedAnIsotope(label) => write!(f,
          "structure_properties does not require an isotope to be set: 
#[spin_properties(label = {})]",label),

      CluEError::TensorNotSet(index) => write!(f,
          "no tensor set for particle {}",index),

      CluEError::TooFewRHSArguments(line_number) => write!(f,
          "line {}, too few arguments on the right hand side", 
          line_number),

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

      CluEError::WrongOrientationGridDim(line_number, expected, actual) 
        => write!(f,
          "line {}, cannot add a {}D, point to a {}D grid", 
          line_number, expected,actual),

      CluEError::WrongProbabilityDistributionDim(line_number, 
          expected, actual) => write!(f,
          "line {}, expected probability distribution with {} dimensions, \
but found {} dimensions", line_number, expected,actual),

      CluEError::WrongVectorLength(line_number, expected, actual) => write!(f,
          "line {}, expected vector of length {}, but recieved a length of {}", 
          line_number, expected,actual),
    }
  }
}
