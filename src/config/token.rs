use crate::clue_errors::*;
use crate::config::token_stream::split_on_token;
use crate::isotopes::Isotope;
use std::ops::{Add,Sub,Mul,Div};

use std::fmt;
//------------------------------------------------------------------------------
/// This enum specifies the tokens that user input text files get converted to.
// After adding a new token to this list, the token should also be added to 
// fmt(), identify_token(), and test_identify_token() below.
#[derive(PartialEq, Debug, Clone)]
pub enum Token{
 Active,
 ApplyPBC,
 ApproxThermal,
 //AtomName,
 AtomNames,
 Bang,
 Ball,
 BondedElements,                                               
 BondedIndices,                                                     
 BlockCommentEnd, 
 BlockCommentStart, 
 CarrPurcell,
 CCE,
 Centroid,
 CentroidOverSerials,
 ClashDistancePBC,
 ClusterBatchSize,
 ClusterDensityMatrix,
 ClusterMethod,
 Clusters,
 ClustersFile,
 Colon,
 Comma,
 Config,
 Cosubstitute,
 Cube,
 CurlyBracketClose,
 CurlyBracketOpen,
 DetectedSpinGMatrix,
 DetectedSpinGX,
 DetectedSpinGY,
 DetectedSpinGZ,
 DetectedSpinIdentity,
 DetectedSpinPosition,
 DetectedSpinMultiplicity,
 DetectedSpinTransition,
 Diff,
 Distance,
 DoubleQuote,
 ElectricQuadrupoleCoupling,
 ElectricQuadrupoleX,
 ElectricQuadrupoleY,
 ElectricQuadrupoleZ,
 //Element,
 Elements,                                                    
 EOL,
 Equals,
 Extracells,
 ExchangeGroupsAndParticles,
 False,
 Filter,
 Float(f64),
 GCCE,
 GMatrix,
 GreaterThan,
 GreaterThanEqualTo,
 Group,
 GX,
 GY,
 GZ,
 Hat,
 Hahn,
 HyperfineCoupling,
 HyperfineX,
 HyperfineY,
 HyperfineZ,
 Identity,
 In,
 Indices,                                                      
 InputStructureFile,
 Int(i32),
 Isotope,
 IsotopeAbundances,
 Label,
 Lebedev,
 LessThan,
 LessThanEqualTo,
 LineComment, 
 LoadGeometry,
 MagneticField,
 MaxCellSize,
 MaxClusterSize,
 MaxSpinOrder,
 MinCellSize,
 Minus,
 Mode(ModeAttribute),
 NeighborCutoffDeltaHyperfine,
 NeighborCutoffCoupling,
 NeighborCutoffDipolePerpendicular,
 NeighborCutoffDistance,
 NeighborCutoff3SpinHahnModDepth,
 NeighborCutoff3SpinHahnTaylor4,
 Not,
 NotEqual,
 NotIn,
 NumberSystemInstances,
 NumberTimepoints,
 OrientationGrid,
 Particles,
 PartitioningMethod,
 Path,
 ParenthesisClose,
 ParenthesisOpen,
 PrimaryCell,
 Plus,
 PulseSequence,
 R2CCE,
 Radius,
 Random,
 ReadCSV,
 Residues,                                                     
 ResSeqNums,                                                   
 RNGSeed,
 SaveDir,
 Semicolon,
 Serials,                                                      
 Set,
 Sharp,
 Slash,
 //Sphere,
 SpinProperties,
 Spin,
 SquareBracketClose,
 SquareBracketOpen,
 StructureProperties,
 Structures,
 SystemName,
 Temperature,
 Tensors,
 Thermal,
 TimeIncrements,
 Times,
 True,
 TunnelSplitting,
 Type,
 UnitOfClustering,
 UserInputValue(String),
 Vector,
 VectorF64(Vec::<f64>),
 VectorI32(Vec::<i32>),
 VectorString(Vec::<String>),
 VoidProbability,
 Whitespace,
 WriteAuxiliarySignals,
 WriteBath,
 WriteClusters,
 WriteExchangeGroups,
 WriteInfo,
 WriteMethylPartitions,
 WriteOrientationSignals,
 WriteSansSpinSignals,
 WriteStructurePDB,
 WriteTensors,
}
impl fmt::Display for Token{
  /// This function translates tokens to strings.
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{
      Token::Active => write!(f,"active"), 
      Token::ApplyPBC => write!(f,"apply_pbc"), 
      Token::ApproxThermal => write!(f,"approx_thermal"),
      //Token::AtomName => write!(f,"atom_name"),
      Token::AtomNames => write!(f,"atom_names"),
      Token::Bang => write!(f,"!"), 
      Token::Ball => write!(f,"ball"), 
      Token::BondedIndices => write!(f,"bonded_indices"),
      Token::BondedElements => write!(f,"bonded_elements"),
      Token::BlockCommentEnd => write!(f,"*/"), 
      Token::BlockCommentStart => write!(f,"/*"), 
      Token::CarrPurcell => write!(f,"cp"),
      Token::CCE => write!(f,"cce"),
      Token::Centroid => write!(f,"centroid"),
      Token::CentroidOverSerials => write!(f,"centroid_over_serials"),
      Token::ClashDistancePBC => write!(f,"clash_distance_pbc"),
      Token::ClusterBatchSize => write!(f,"cluster_batch_size"),
      Token::ClusterDensityMatrix => write!(f,"cluster_density_matrix"),
      Token::ClusterMethod => write!(f,"cluster_method"),
      Token::Clusters => write!(f,"clusters"),
      Token::ClustersFile => write!(f,"input_clusters_file"),
      Token::Colon => write!(f,":"),
      Token::Comma => write!(f,","),
      Token::Config => write!(f,"config"),
      Token::Cosubstitute => write!(f,"cosubstitute"),
      Token::Cube => write!(f,"cube"),
      Token::CurlyBracketClose => write!(f,"}}"),
      Token::CurlyBracketOpen => write!(f,"{{"),
      Token::DetectedSpinGMatrix => write!(f,"detected_spin_g_matrix"),
      Token::DetectedSpinGX => write!(f,"detected_spin_g_x"),
      Token::DetectedSpinGY => write!(f,"detected_spin_g_y"),
      Token::DetectedSpinGZ => write!(f,"detected_spin_g_z"),
      Token::DetectedSpinIdentity => write!(f,"detected_spin_identity"),
      Token::DetectedSpinMultiplicity => write!(f,"detected_spin_multiplicity"),
      Token::DetectedSpinPosition => write!(f,"detected_spin_position"),
      Token::DetectedSpinTransition => write!(f,"detected_spin_transition"),
      Token::Diff => write!(f,"diff"),
      Token::Distance => write!(f,"distance"),
      Token::DoubleQuote => write!(f,"\""),
      Token::ElectricQuadrupoleCoupling  
        => write!(f,"electric_quadrupole_coupling"),
      Token::ElectricQuadrupoleX  => write!(f,"electric_quadrupole_x"),
      Token::ElectricQuadrupoleY  => write!(f,"electric_quadrupole_y"),
      Token::ElectricQuadrupoleZ  => write!(f,"electric_quadrupole_z"),
      //Token::Element => write!(f,"element"),
      Token::Elements => write!(f,"elements"),
      Token::EOL => writeln!(f),
      Token::Equals => write!(f,"="),
      Token::ExchangeGroupsAndParticles 
        => write!(f,"exchange_groups_and_particles"),
      Token::Extracells => write!(f,"extracells"),
      Token::False => writeln!(f,"false"),
      Token::Filter => write!(f,"filter"),
      Token::Float(x) => write!(f,"{}",x),
      Token::GCCE => write!(f,"gcce"),
      Token::GMatrix => write!(f,"g_matrix"),
      Token::GreaterThan => write!(f,">"),
      Token::GreaterThanEqualTo => write!(f,">="),
      Token::Group => write!(f,"group"),
      Token::GX => write!(f,"g_x"),
      Token::GY => write!(f,"g_y"),
      Token::GZ => write!(f,"g_z"),
      Token::Hahn => write!(f,"hahn"),
      Token::Hat => write!(f,"^"),
      Token::HyperfineCoupling => write!(f,"hyperfine_coupling"),
      Token::HyperfineX => write!(f,"hyperfine_x"),
      Token::HyperfineY => write!(f,"hyperfine_y"),
      Token::HyperfineZ => write!(f,"hyperfine_z"),
      Token::Identity => write!(f,"identity"),
      Token::In => write!(f,"in"),
      Token::Indices => write!(f,"indices"),
      Token::InputStructureFile => write!(f,"input_structure_file"),
      Token::Int(n) => write!(f,"{}",n),
      Token::Isotope => write!(f,"isotope"),
      Token::IsotopeAbundances => write!(f,"isotope_abundances"),
      Token::Label => write!(f,"label"),
      Token::Lebedev => write!(f,"lebedev"),
      Token::LessThan => write!(f,"<"),
      Token::LessThanEqualTo => write!(f,"<="),
      Token::LineComment => write!(f,"//"), 
      Token::LoadGeometry => write!(f,"load_geometry"),
      Token::MagneticField => write!(f,"magnetic_field"),
      Token::MaxCellSize => write!(f,"max_cell_size"),
      Token::MaxClusterSize => write!(f,"max_cluster_size"),
      Token::MaxSpinOrder => write!(f,"max_spin_order"),
      Token::MinCellSize => write!(f,"min_cell_size"),
      Token::Minus => write!(f,"-"),
      Token::Mode(mode) => write!(f,"{}",mode),
      Token::NeighborCutoffDeltaHyperfine 
        => write!(f,"neighbor_cutoff_delta_hyperfine"),
      Token::NeighborCutoffCoupling 
        => write!(f,"neighbor_cutoff_coupling"),
      Token::NeighborCutoffDipolePerpendicular 
        => write!(f,"neighbor_cutoff_dipole_perpendicular"),
      Token::NeighborCutoffDistance 
        => write!(f,"neighbor_cutoff_distance"),
      Token::NeighborCutoff3SpinHahnModDepth 
        => write!(f,"neighbor_cutoff_3_spin_hahn_mod_depth"),
      Token::NeighborCutoff3SpinHahnTaylor4 
        => write!(f,"neighbor_cutoff_3_spin_hahn_taylor_4"),
      Token::Not => write!(f,"not"),
      Token::NotEqual => write!(f,"!="),
      Token::NotIn => write!(f,"not in"),
      Token::NumberSystemInstances => write!(f,"number_system_instances"),
      Token::NumberTimepoints => write!(f,"number_timepoints"),
      Token::OrientationGrid => write!(f,"orientation_grid"),
      Token::Particles => write!(f,"particles"),
      Token::PartitioningMethod => write!(f,"partition_method"),
      Token::Path => write!(f,"path"),
      Token::ParenthesisClose => write!(f,")"),
      Token::ParenthesisOpen => write!(f,"("),
      Token::PrimaryCell => write!(f,"primary_cell"),
      Token::Plus => write!(f,"+"),
      Token::PulseSequence => write!(f,"pulse_sequence"),
      Token::R2CCE => write!(f,"r2cce"),
      Token::Radius => write!(f,"radius"),
      Token::Random => write!(f,"random"),
      Token::ReadCSV => write!(f,"read_csv"),
      Token::Residues => write!(f,"residues"),
      Token::ResSeqNums => write!(f,"residue_sequence_numbers"),
      Token::RNGSeed => write!(f,"rng_seed"),
      Token::SaveDir => write!(f,"save_dir"),
      Token::Semicolon => write!(f,";"),
      Token::Serials => write!(f,"serials"),
      Token::Set => write!(f,"set"),
      Token::Sharp => write!(f,"#"),
      Token::Slash => write!(f,"/"),
      Token::SpinProperties => write!(f,"spin_properties"),
      Token::Spin => write!(f,"spin"),
      //Token::Sphere => write!(f,"sphere"),
      Token::SquareBracketClose => write!(f,"]"),
      Token::SquareBracketOpen => write!(f,"["),
      Token::StructureProperties => write!(f,"structure_properties"),
      Token::Structures => write!(f,"structures"),
      Token::SystemName => write!(f,"system_name"),
      Token::Temperature => write!(f,"temperature"),
      Token::Tensors => write!(f,"tensors"),
      Token::Thermal => write!(f,"thermal"),
      Token::TimeIncrements => write!(f,"time_increments"),
      Token::Times => write!(f,"*"),
      Token::True => writeln!(f,"true"),
      Token::TunnelSplitting => write!(f,"tunnel_splitting"),
      Token::Type => write!(f,"type"),
      Token::UnitOfClustering => write!(f,"unit_of_clustering"),
      Token::UserInputValue(string) => write!(f,"{}",string),
      Token::Vector => write!(f,"vector"), 
      Token::VectorF64(v) => write!(f,"{:?}",v), 
      Token::VectorI32(v) => write!(f,"{:?}",v), 
      Token::VectorString(v) => write!(f,"{:?}",v), 
      Token::VoidProbability => write!(f,"void_probability"),
      Token::Whitespace => write!(f," "),
      Token::WriteAuxiliarySignals => write!(f,"write_auxiliary_signals"),
      Token::WriteBath => write!(f,"write_bath"),
      Token::WriteClusters => write!(f,"write_clusters"),
      Token::WriteExchangeGroups => write!(f,"write_exchange_groups"),
      Token::WriteInfo => write!(f,"write_info"),
      Token::WriteMethylPartitions => write!(f,"write_methyl_partitions"),
      Token::WriteOrientationSignals => write!(f,"write_orientation_signals"),
      Token::WriteSansSpinSignals => write!(f,"write_sans_spin_signals"),
      Token::WriteStructurePDB => write!(f,"write_structure_pdb"),
      Token::WriteTensors => write!(f,"write_tensors"),
    }
  }
}
//------------------------------------------------------------------------------
/// This function translates strings to tokens.
pub fn identify_token(word: &str) -> Option<Token>{
  match word{
    "active" => Some(Token::Active),
    "apply_pbc" => Some(Token::ApplyPBC),
    "approx_thermal" => Some(Token::ApproxThermal),
    //"atom_name" => Some(Token::AtomName),
    "atom_names" => Some(Token::AtomNames),
    "!" => Some(Token::Bang),
    "ball" => Some(Token::Ball),
    "bonded_indices" => Some(Token::BondedIndices),
    "bonded_elements" => Some(Token::BondedElements),
    "*/" => Some(Token::BlockCommentEnd),
    "/*" => Some(Token::BlockCommentStart),
    "cp" => Some(Token::CarrPurcell),
    "cce" => Some(Token::CCE),
    "centroid" => Some(Token::Centroid),
    "centroid_over_serials" => Some(Token::CentroidOverSerials),
    "clash_distance_pbc" => Some(Token::ClashDistancePBC),
    "cluster_batch_size" => Some(Token::ClusterBatchSize),
    "cluster_density_matrix" => Some(Token::ClusterDensityMatrix),
    "cluster_method" => Some(Token::ClusterMethod),
    "clusters" => Some(Token::Clusters),
    "input_clusters_file" => Some(Token::ClustersFile),
    ":" => Some(Token::Colon),
    "," => Some(Token::Comma),
    "config" => Some(Token::Config),
    "cosubstitute" => Some(Token::Cosubstitute),
    "cube" => Some(Token::Cube),
    "}" => Some(Token::CurlyBracketClose),
    "{" => Some(Token::CurlyBracketOpen),
    "detected_spin_g_matrix" => Some(Token::DetectedSpinGMatrix),
    "detected_spin_g_x" => Some(Token::DetectedSpinGX),
    "detected_spin_g_y" => Some(Token::DetectedSpinGY),
    "detected_spin_g_z" => Some(Token::DetectedSpinGZ),
    "detected_spin_identity" => Some(Token::DetectedSpinIdentity),
    "detected_spin_multiplicity" => Some(Token::DetectedSpinMultiplicity),
    "detected_spin_position" => Some(Token::DetectedSpinPosition),
    "detected_spin_transition" => Some(Token::DetectedSpinTransition),
    "diff" => Some(Token::Diff),
    "distance" => Some(Token::Distance),
    "\"" => Some(Token::DoubleQuote),
    "electric_quadrupole_coupling" => Some(Token::ElectricQuadrupoleCoupling),
    "electric_quadrupole_x" => Some(Token::ElectricQuadrupoleX),
    "electric_quadrupole_y" => Some(Token::ElectricQuadrupoleY),
    "electric_quadrupole_z" => Some(Token::ElectricQuadrupoleZ),
    //"element" => Some(Token::Element),
    "elements" => Some(Token::Elements), 
    "\n" => Some(Token::EOL),
    "=" => Some(Token::Equals),
    "exchange_groups_and_particles" => Some(Token::ExchangeGroupsAndParticles),
    "extracells" => Some(Token::Extracells),
    "false" => Some(Token::False),
    "filter" => Some(Token::Filter),
    "gcce" => Some(Token::GCCE),
    "g_matrix" => Some(Token::GMatrix),
    ">" => Some(Token::GreaterThan),
    ">=" => Some(Token::GreaterThanEqualTo),
    "group" => Some(Token::Group),
    "g_x" => Some(Token::GX),
    "g_y" => Some(Token::GY),
    "g_z" => Some(Token::GZ),
    "identity" => Some(Token::Identity),
    "in" => Some(Token::In),
    "input_structure_file" => Some(Token::InputStructureFile),
    "indices" => Some(Token::Indices), 
    "isotope" => Some(Token::Isotope),
    "isotope_abundances" => Some(Token::IsotopeAbundances),
    "label" => Some(Token::Label),
    "load_geometry" => Some(Token::LoadGeometry),
    "lebedev" => Some(Token::Lebedev),
    "hahn" => Some(Token::Hahn),
    "^" => Some(Token::Hat),
    "hyperfine_coupling" => Some(Token::HyperfineCoupling),
    "hyperfine_x" => Some(Token::HyperfineX),
    "hyperfine_y" => Some(Token::HyperfineY),
    "hyperfine_z" => Some(Token::HyperfineZ),
    "<" => Some(Token::LessThan),
    "<=" => Some(Token::LessThanEqualTo),
    "//" => Some(Token::LineComment),
    "magnetic_field" => Some(Token::MagneticField),
    "max_cell_size" => Some(Token::MaxCellSize),
    "max_cluster_size" => Some(Token::MaxClusterSize),
    "max_spin_order" => Some(Token::MaxSpinOrder),
    "min_cell_size" => Some(Token::MinCellSize),
    "-" => Some(Token::Minus),
    "neighbor_cutoff_delta_hyperfine" 
      => Some(Token::NeighborCutoffDeltaHyperfine),
    "neighbor_cutoff_coupling" 
      => Some(Token::NeighborCutoffCoupling),
    "neighbor_cutoff_dipole_perpendicular" 
      => Some(Token::NeighborCutoffDipolePerpendicular),
    "neighbor_cutoff_distance" 
      => Some(Token::NeighborCutoffDistance),
    "neighbor_cutoff_3_spin_hahn_mod_depth" 
      => Some(Token::NeighborCutoff3SpinHahnModDepth),
    "neighbor_cutoff_3_spin_hahn_taylor_4" 
      => Some(Token::NeighborCutoff3SpinHahnTaylor4),
    "not" => Some(Token::Not),
    "!=" => Some(Token::NotEqual),
    "not in" => Some(Token::NotIn),
    "number_system_instances" => Some(Token::NumberSystemInstances),
    "number_timepoints" => Some(Token::NumberTimepoints),
    "orientation_grid" => Some(Token::OrientationGrid),
    "particles" => Some(Token::Particles),
    "partitioning_method" => Some(Token::PartitioningMethod),
    "path" => Some(Token::Path),
    ")" => Some(Token::ParenthesisClose),
    "(" => Some(Token::ParenthesisOpen),
    "+" => Some(Token::Plus),
    "primary_cell" => Some(Token::PrimaryCell),
    "pulse_sequence" => Some(Token::PulseSequence),
    "r2cce" => Some(Token::R2CCE),
    "radius" => Some(Token::Radius),
    "random" => Some(Token::Random),
    "read_csv" => Some(Token::ReadCSV),
    "residues" => Some(Token::Residues),
    "residue_sequence_numbers" => Some(Token::ResSeqNums),
    "rng_seed" => Some(Token::RNGSeed),
    "save_dir" => Some(Token::SaveDir),
    ";" => Some(Token::Semicolon),
    "serials" => Some(Token::Serials), 
    "set" => Some(Token::Set), 
    "#" => Some(Token::Sharp),
    "/" => Some(Token::Slash),
    "spin_properties" => Some(Token::SpinProperties),
    //"sphere" => Some(Token::Sphere),
    "spin" => Some(Token::Spin),
    "]" => Some(Token::SquareBracketClose),
    "[" => Some(Token::SquareBracketOpen),
    "structure_properties" => Some(Token::StructureProperties),
    "structures" => Some(Token::Structures),
    "system_name" => Some(Token::SystemName),
    "temperature" => Some(Token::Temperature),
    "tensors" => Some(Token::Tensors),
    "thermal" => Some(Token::Thermal),
    "*" => Some(Token::Times),
    "true" => Some(Token::True),
    "time_increments" => Some(Token::TimeIncrements),
    "tunnel_splitting" => Some(Token::TunnelSplitting),
    "type" => Some(Token::Type),
    "unit_of_clustering" => Some(Token::UnitOfClustering),
    "vector" => Some(Token::Vector),
    "void_probability" => Some(Token::VoidProbability),
    " " => Some(Token::Whitespace),
    "write_auxiliary_signals" => Some(Token::WriteAuxiliarySignals),
    "write_bath" => Some(Token::WriteBath),
    "write_clusters" => Some(Token::WriteClusters),
    "write_exchange_groups" => Some(Token::WriteExchangeGroups),
    "write_info" => Some(Token::WriteInfo),
    "write_methyl_partitions" => Some(Token::WriteMethylPartitions),
    "write_orientation_signals" => Some(Token::WriteOrientationSignals),
    "write_sans_spin_signals" => Some(Token::WriteSansSpinSignals),
    "write_structure_pdb" => Some(Token::WriteStructurePDB),
    "write_tensors" => Some(Token::WriteTensors),
    _ => None
  }
}

impl Add for Token{
  
  type Output = Result<Token,CluEError>;
  /// This function implements "+" for tokens.
  /// The function will return `Err` if the tokens cannot be added.
  fn add(self, rhs: Token) -> Result<Token,CluEError>{
    
    match (self, rhs){
      (Token::Float(x), Token::Float(y)) => Ok(Token::Float(x+y)),
      (Token::Float(x), Token::Int(b)) => Ok(Token::Float(x + b as f64)),
      (Token::Float(x), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| x+y).collect()  )),
      (Token::Float(x), Token::VectorI32(v)) => 
        Ok(Token::VectorF64( v.iter().map(|b| x + *b as f64).collect()  )),

      (Token::Int(a), Token::Float(y)) => Ok(Token::Float(a as f64 + y)),
      (Token::Int(a), Token::Int(b)) => Ok(Token::Int(a+b)),
      (Token::Int(a), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| a as f64 + y).collect()  )),
      (Token::Int(a), Token::VectorI32(v)) => 
        Ok(Token::VectorI32( v.iter().map(|b| a+b).collect()  )),
      
      (Token::VectorF64(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x+y).collect()  )),
      (Token::VectorF64(u), Token::Int(b)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x+b as f64).collect()  )),
      (Token::VectorF64(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotAddTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] + v[ii]); }
        Ok(Token::VectorF64(w))
      }, 
      (Token::VectorF64(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotAddTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] + v[ii] as f64); }
        Ok(Token::VectorF64(w))
      }, 

      (Token::VectorI32(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|a| *a as f64 +y).collect()  )),
      (Token::VectorI32(u), Token::Int(b)) => 
        Ok(Token::VectorI32( u.iter().map(|a| a+b).collect()  )),
      (Token::VectorI32(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotAddTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] as f64 + v[ii]); }
        Ok(Token::VectorF64(w))
      },
      (Token::VectorI32(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotAddTokens);}
        let mut w = Vec::<i32>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] + v[ii]); }
        Ok(Token::VectorI32(w))
      }, 
      (_,_) => Err(CluEError::CannotAddTokens),
    }
  }
}

impl Sub for Token{
  
  type Output = Result<Token,CluEError>;
  /// This function implements "-" for tokens.
  /// The function will return `Err` if the tokens cannot be subtracted.
  fn sub(self, rhs: Token) -> Result<Token,CluEError>{
    
    match (self, rhs){
      (Token::Float(x), Token::Float(y)) => Ok(Token::Float(x-y)),
      (Token::Float(x), Token::Int(b)) => Ok(Token::Float(x - b as f64)),
      (Token::Float(x), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| x-y).collect()  )),
      (Token::Float(x), Token::VectorI32(v)) => 
        Ok(Token::VectorF64( v.iter().map(|b| x - *b as f64).collect()  )),

      (Token::Int(a), Token::Float(y)) => Ok(Token::Float(a as f64 - y)),
      (Token::Int(a), Token::Int(b)) => Ok(Token::Int(a-b)),
      (Token::Int(a), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| a as f64 - y).collect()  )),
      (Token::Int(a), Token::VectorI32(v)) => 
        Ok(Token::VectorI32( v.iter().map(|b| a-b).collect()  )),
      
      (Token::VectorF64(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x-y).collect()  )),
      (Token::VectorF64(u), Token::Int(b)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x-b as f64).collect()  )),
      (Token::VectorF64(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotSubTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] - v[ii]); }
        Ok(Token::VectorF64(w))
      }, 
      (Token::VectorF64(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotSubTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] - v[ii] as f64); }
        Ok(Token::VectorF64(w))
      }, 

      (Token::VectorI32(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|a| *a as f64 -y).collect()  )),
      (Token::VectorI32(u), Token::Int(b)) => 
        Ok(Token::VectorI32( u.iter().map(|a| a-b).collect()  )),
      (Token::VectorI32(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotSubTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] as f64 - v[ii]); }
        Ok(Token::VectorF64(w))
      },
      (Token::VectorI32(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotSubTokens);}
        let mut w = Vec::<i32>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] - v[ii]); }
        Ok(Token::VectorI32(w))
      }, 
      (_,_) => Err(CluEError::CannotSubTokens),
    }
  }
}

impl Mul for Token{
  
  type Output = Result<Token,CluEError>;
  /// This function implements "*" for tokens.
  /// The function will return `Err` if the tokens cannot be multiplied.
  fn mul(self, rhs: Token) -> Result<Token,CluEError>{
    
    match (self, rhs){
      (Token::Float(x), Token::Float(y)) => Ok(Token::Float(x*y)),
      (Token::Float(x), Token::Int(b)) => Ok(Token::Float(x * b as f64)),
      (Token::Float(x), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| x*y).collect()  )),
      (Token::Float(x), Token::VectorI32(v)) => 
        Ok(Token::VectorF64( v.iter().map(|b| x * *b as f64).collect()  )),

      (Token::Int(a), Token::Float(y)) => Ok(Token::Float(a as f64 * y)),
      (Token::Int(a), Token::Int(b)) => Ok(Token::Int(a*b)),
      (Token::Int(a), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| a as f64 * y).collect()  )),
      (Token::Int(a), Token::VectorI32(v)) => 
        Ok(Token::VectorI32( v.iter().map(|b| a*b).collect()  )),
      
      (Token::VectorF64(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x*y).collect()  )),
      (Token::VectorF64(u), Token::Int(b)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x*b as f64).collect()  )),
      (Token::VectorF64(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotMulTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] * v[ii]); }
        Ok(Token::VectorF64(w))
      }, 
      (Token::VectorF64(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotMulTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] * v[ii] as f64); }
        Ok(Token::VectorF64(w))
      }, 

      (Token::VectorI32(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|a| *a as f64 *y).collect()  )),
      (Token::VectorI32(u), Token::Int(b)) => 
        Ok(Token::VectorI32( u.iter().map(|a| a*b).collect()  )),
      (Token::VectorI32(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotMulTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] as f64 * v[ii]); }
        Ok(Token::VectorF64(w))
      },
      (Token::VectorI32(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotMulTokens);}
        let mut w = Vec::<i32>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] * v[ii]); }
        Ok(Token::VectorI32(w))
      }, 
      (_,_) => Err(CluEError::CannotMulTokens),
    }
  }
}

impl Div for Token{
  
  type Output = Result<Token,CluEError>;
  /// This function implements "/" for tokens.
  /// The function will return `Err` if the tokens cannot be divided.
  fn div(self, rhs: Token) -> Result<Token,CluEError>{
    
    match (self, rhs){
      (Token::Float(x), Token::Float(y)) => Ok(Token::Float(x/y)),
      (Token::Float(x), Token::Int(b)) => Ok(Token::Float(x / b as f64)),
      (Token::Float(x), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| x/y).collect()  )),
      (Token::Float(x), Token::VectorI32(v)) => 
        Ok(Token::VectorF64( v.iter().map(|b| x / *b as f64).collect()  )),

      (Token::Int(a), Token::Float(y)) => Ok(Token::Float(a as f64 / y)),
      (Token::Int(a), Token::Int(b)) => Ok(Token::Int(a/b)),
      (Token::Int(a), Token::VectorF64(v)) => 
        Ok(Token::VectorF64( v.iter().map(|y| a as f64 / y).collect()  )),
      (Token::Int(a), Token::VectorI32(v)) => 
        Ok(Token::VectorI32( v.iter().map(|b| a/b).collect()  )),
      
      (Token::VectorF64(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x/y).collect()  )),
      (Token::VectorF64(u), Token::Int(b)) => 
        Ok(Token::VectorF64( u.iter().map(|x| x/b as f64).collect()  )),
      (Token::VectorF64(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotDivTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] / v[ii]); }
        Ok(Token::VectorF64(w))
      }, 
      (Token::VectorF64(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotDivTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] / v[ii] as f64); }
        Ok(Token::VectorF64(w))
      }, 

      (Token::VectorI32(u), Token::Float(y)) => 
        Ok(Token::VectorF64( u.iter().map(|a| *a as f64 /y).collect()  )),
      (Token::VectorI32(u), Token::Int(b)) => 
        Ok(Token::VectorI32( u.iter().map(|a| a/b).collect()  )),
      (Token::VectorI32(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotDivTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] as f64 / v[ii]); }
        Ok(Token::VectorF64(w))
      },
      (Token::VectorI32(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotDivTokens);}
        let mut w = Vec::<i32>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] / v[ii]); }
        Ok(Token::VectorI32(w))
      }, 
      (_,_) => Err(CluEError::CannotDivTokens),
    }
  }
}

impl Token{
  /// This function allows one token to be raised to the power of another
  /// token.  The output is `Err` if the tokens do not allow for exponentiation.
  pub fn pow(self, rhs: Token) -> Result<Token,CluEError>{

    match (self, rhs){
      (Token::Float(x), Token::Float(y)) => Ok(Token::Float(x.powf(y) )),
      (Token::Float(x), Token::Int(b)) => Ok(Token::Float(x.powi(b) )),
      (Token::Float(x), Token::VectorF64(v)) =>
        Ok(Token::VectorF64( v.iter().map(|y| x.powf(*y) ).collect()  )),
      (Token::Float(x), Token::VectorI32(v)) =>
        Ok(Token::VectorF64( v.iter().map(|b| x.powi(*b)).collect()  )),

      (Token::Int(a), Token::Float(y)) => Ok(Token::Float((a as f64).powf(y) )),
      (Token::Int(a), Token::Int(b)) =>{ 
        if b < 0{ return Ok(Token::Float((a as f64).powi(b) )); }
        Ok(Token::Int(a.pow(b as u32) )) 
      },
      (Token::Int(a), Token::VectorF64(v)) =>
        Ok(Token::VectorF64( v.iter().map(|y| (a as f64).powf(*y) ).collect())),
      (Token::Int(a), Token::VectorI32(v)) =>{
        let mut any_neg = false;
        for b in v.iter(){
          if *b<0{
            any_neg = true;
            break;
          }
        }
        if any_neg{
          return Ok(Token::VectorF64( v.iter().map(
                |b| (a as f64).powi( *b )).collect()  ));
        }
        Ok(Token::VectorI32( v.iter().map(
                |b| a.pow((*b) as u32) ).collect()  ))
      },

      (Token::VectorF64(u), Token::Float(y)) =>
        Ok(Token::VectorF64( u.iter().map(|x| x.powf(y) ).collect()  )),
      (Token::VectorF64(u), Token::Int(b)) =>
        Ok(Token::VectorF64( u.iter().map(|x| x.powi(b) ).collect()  )),
      (Token::VectorF64(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotPowTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii] .powf( v[ii]) ); }
        Ok(Token::VectorF64(w))
      },
      (Token::VectorF64(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotPowTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push(u[ii].powi(v[ii]) ); }
        Ok(Token::VectorF64(w))
      },

      (Token::VectorI32(u), Token::Float(y)) =>
        Ok(Token::VectorF64( u.iter().map(|a| (*a as f64).powf(y) ).collect()  )),
      (Token::VectorI32(u), Token::Int(b)) =>
      {
        if b < 0{
          return Ok(Token::VectorF64( u.iter().map(
                |a| (*a as f64).powi(b) ).collect()  ));
        }
        Ok(Token::VectorI32( u.iter().map(
                |a| a.pow(b as u32) ).collect()  ))
      },
      (Token::VectorI32(u), Token::VectorF64(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotPowTokens);}
        let mut w = Vec::<f64>::with_capacity(u.len());
        for ii in 0..u.len(){ w.push((u[ii] as f64).powf(v[ii]) ); }
        Ok(Token::VectorF64(w))
      },
      (Token::VectorI32(u), Token::VectorI32(v)) => {
        if u.len() != v.len(){ return Err(CluEError::CannotPowTokens);}

        let mut any_neg = false;
        for b in v.iter(){
          if *b<0{
            any_neg = true;
            break;
          }
        }
        if any_neg{
          let mut w = Vec::<f64>::with_capacity(u.len());
          for ii in 0..u.len(){ 
            w.push((u[ii] as f64).powi(v[ii]) ); 
          }
          return Ok(Token::VectorF64(w))
        }

        let mut w = Vec::<i32>::with_capacity(u.len());
        for ii in 0..u.len(){ 
          w.push(u[ii].pow(v[ii] as u32) ); 
        }
        Ok(Token::VectorI32(w))
      },
      (_,_) => Err(CluEError::CannotPowTokens),

    }
  }
}
//------------------------------------------------------------------------------
/// This function answers if the token specifies a relationship.
pub fn is_relational_operator(token: &Token) -> bool{
  matches!( token,
    Token::Equals | Token::GreaterThan | Token::GreaterThanEqualTo
    | Token::In | Token::LessThan | Token::LessThanEqualTo | Token::NotEqual 
    | Token::NotIn
    )
}
//------------------------------------------------------------------------------
/// This function answers whether or not the token is "(", "[" or "{".
pub fn is_opening_delimiter(token: &Token) -> bool{
  matches!(token,
    Token::ParenthesisOpen | Token::SquareBracketOpen | Token::CurlyBracketOpen
    )
}
//------------------------------------------------------------------------------
/// This function answers whether or not the token is ")", "]" or "}".
pub fn is_closing_delimiter(token: &Token) -> bool{
  matches!(token,
    Token::ParenthesisClose | Token::SquareBracketClose 
    | Token::CurlyBracketClose
    )
}
//------------------------------------------------------------------------------
/// This function outputs a list of indices where the target token can be found
/// within a list of tokens.
pub fn find_token(target: &Token, tokens: &[Token]) -> Vec::<usize>{

  let mut out = Vec::<usize>::with_capacity(count_token(target,tokens) );

  for (ii, token) in tokens.iter().enumerate(){
    if *token == *target{
      out.push(ii);
    }
  }
  out
}
//------------------------------------------------------------------------------
/// This function outputs a the number of occurances of a target token
/// within a list of tokens.
pub fn count_token(target: &Token, tokens: &[Token]) -> usize{

  let mut counter = 0;

  for token in tokens.iter(){
    if token == target{
      counter += 1;
    }
  }
  counter   
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// 'ModeAttribute' keeps tracks of the current parsing mode:
/// "\#\[config\]", "\#\[group(...)\]", etc.
/// The other fields hold potential arguments of the mode header.
#[derive(PartialEq, Debug, Clone)]
pub struct ModeAttribute{
  pub mode:  ConfigMode,
  pub isotope: Option<Isotope>,
  pub label: Option<String>,
  pub path: Option<String>,
}

impl Default for ModeAttribute{
  // This function implements the `default()` for `ModeAttribute`.
  fn default() -> Self{
    ModeAttribute{
      mode:  ConfigMode::Config,
      isotope: None,
      label: None,
      path: None,
    }
  }
}

impl ModeAttribute{

  /// This function implements the `new()` for `ModeAttribute`.
  pub fn new() -> Self{
    ModeAttribute::default()
  }

  /// This function tries to interpret a list of tokens as a `ModeAttribute`.
  /// The output will be `Err` if it cannot.
  pub fn from(tokens: Vec::<Token>) -> Result<Self,CluEError>
  {

    // Locate "#"s within tokens.
    let sharp_indices = find_token(&Token::Sharp,&tokens);
    if sharp_indices.len() != 1 || sharp_indices[0] != 0{
      return Err(CluEError::ModeAttributeWrongSharp);
    } 

    // Check formating.
    let idx = sharp_indices[0];
    if tokens[idx+1] != Token::SquareBracketOpen
    && tokens[tokens.len()-1] != Token::SquareBracketClose{
      return Err(CluEError::ModeAttributeWrongBrackets);
    }

    // Determine which mode.
    let mode = ConfigMode::from(tokens[idx+2].clone())?; 

    let args = if tokens.len() > 4{
      split_on_token(tokens[idx+4..tokens.len()-2].to_vec(), 
        Token::Comma)
    }else{
      Vec::<Vec::<Token>>::new()
    };

    let mut allow_positional_argument = true;

    // Find label identifier.
    let label: Option<String>;
    let label_indices = find_token(&Token::Label,&tokens);
    
    if label_indices.is_empty(){
      if !args.is_empty(){
        label = Some(args[0][0].to_string());
      }else{
        label = None;
      }
    } else if label_indices.len() == 1 
      && tokens[label_indices[0]+1] == Token::Equals{
      
      allow_positional_argument = false;

      if let Token::UserInputValue(value)=&tokens[label_indices[0]+2]{
        label = Some(value.clone());
      }else{
        return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
      }
    } else {
      return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
    }

    // Find isotope identifier.
    let isotope: Option<Isotope>;
    let isotope_indices = find_token(&Token::Isotope, &tokens);

    if isotope_indices.is_empty(){
      if args.len() >= 2{
        if !allow_positional_argument{
          return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
        }
        isotope = Some(Isotope::from(&args[1][0].to_string())?);
      }else{
        isotope = None;
      }
    } else if isotope_indices.len() == 1 
      && tokens[isotope_indices[0]+1] == Token::Equals{
      
       allow_positional_argument = false;

       if  let Token::UserInputValue(value)=&tokens[isotope_indices[0]+2]{
        let isotope_value = Isotope::from(value)?;
        isotope = Some(isotope_value);
      }else{
        return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
      }
       
    } else {
      return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
    }

    // Find path identifier.
    let path: Option<String>;
    let path_indices = find_token(&Token::Path, &tokens);

    if path_indices.is_empty(){
      if args.len() >= 3{
        if !allow_positional_argument{
          return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
        }
        path = Some(args[2][0].to_string());
      }else{
        path = None;
      }
    } else if path_indices.len() == 1 
      && tokens[path_indices[0]+1] == Token::Equals{
      
       //allow_positional_argument = false;

       if  let Token::UserInputValue(value)=&tokens[path_indices[0]+2]{
        path = Some(value.clone());
      }else{
        return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
      }
       
    } else {
      return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
    }

    Ok(ModeAttribute{ mode,isotope, label, path})
  }
  //----------------------------------------------------------------------------
}
impl fmt::Display for ModeAttribute{
  // This function translates `ModeAttribute` to a string.
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
  
    let mode_str = self.mode.to_string();

    let label_str: String;
    if let Some(value) = &self.label{
      label_str = format!(", label = {}",value);
    }else{
      label_str = "".to_string();
    }

    let path_str: String;
    if let Some(value) = &self.path{
      path_str = format!(", path = {}",value);
    }else{
      path_str = "".to_string();
    }

    write!(f,"#[{}{}{}]",mode_str,label_str,path_str)

  }
}

/// This enum list the possible configuration modes.
#[derive(PartialEq, Debug, Clone)]
pub enum ConfigMode{
  Config,
  Filter,
  SpinProperties,
  StructureProperties,
}

impl ConfigMode{
  /// This function tries to interpret a list of tokens as a `ConfigMode`.
  /// The output will be `Err` if it cannot.
  pub fn from(token: Token) -> Result<Self,CluEError>{

    match token{
      Token::Config => Ok(ConfigMode::Config),
      Token::Filter | Token::Group => Ok(ConfigMode::Filter),
      Token::SpinProperties => Ok(ConfigMode::SpinProperties),
      Token::StructureProperties => Ok(ConfigMode::StructureProperties),
      _ => Err(CluEError::ConfigModeNotRecognized(token.to_string())),
    }
  }
  //----------------------------------------------------------------------------
}
impl fmt::Display for ConfigMode{
  // This function translates a `ConfigMode` to a string.
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{
      ConfigMode::Config => write!(f,"config"),
      ConfigMode::Filter => write!(f,"filter"),
      ConfigMode::SpinProperties => write!(f,"spin_properties"),
      ConfigMode::StructureProperties => write!(f,"structure_properties"),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  
  #[test]
  fn test_identify_token(){
    assert_eq!(identify_token("active"), Some(Token::Active));
    assert_eq!(identify_token("apply_pbc"), Some(Token::ApplyPBC));
    assert_eq!(identify_token("approx_thermal"), Some(Token::ApproxThermal));
    //assert_eq!(identify_token("atom_name"), Some(Token::AtomName));
    assert_eq!(identify_token("atom_names"), Some(Token::AtomNames));
    assert_eq!(identify_token("!"), Some(Token::Bang));
    assert_eq!(identify_token("ball"), Some(Token::Ball));
    assert_eq!(identify_token("*/"), Some(Token::BlockCommentEnd));
    assert_eq!(identify_token("/*"),Some( Token::BlockCommentStart));
    assert_eq!(identify_token("cp"),Some( Token::CarrPurcell));
    assert_eq!(identify_token("cce"),Some( Token::CCE));
    assert_eq!(identify_token("centroid"),
        Some(Token::Centroid));
    assert_eq!(identify_token("centroid_over_serials"),
        Some(Token::CentroidOverSerials));
    assert_eq!(identify_token("clash_distance_pbc"),
        Some(Token::ClashDistancePBC));
    assert_eq!(identify_token("cluster_batch_size"),
        Some(Token::ClusterBatchSize));
    assert_eq!(identify_token("cluster_density_matrix"),
        Some(Token::ClusterDensityMatrix));
    assert_eq!(identify_token("cluster_method"),Some(Token::ClusterMethod));
    assert_eq!(identify_token("clusters"),Some( Token::Clusters));
    assert_eq!(identify_token("input_clusters_file"),Some( Token::ClustersFile));
    assert_eq!(identify_token(":"),Some( Token::Colon));
    assert_eq!(identify_token(","),Some( Token::Comma));
    assert_eq!(identify_token("config"),Some( Token::Config));
    assert_eq!(identify_token("cosubstitute"),Some( Token::Cosubstitute));
    assert_eq!(identify_token("cube"),Some( Token::Cube));
    assert_eq!(identify_token("}"),Some( Token::CurlyBracketClose));
    assert_eq!(identify_token("{"),Some( Token::CurlyBracketOpen));
    assert_eq!(identify_token("detected_spin_g_matrix"),
        Some(Token::DetectedSpinGMatrix));
    assert_eq!(identify_token("detected_spin_g_x"),Some(Token::DetectedSpinGX));
    assert_eq!(identify_token("detected_spin_g_y"),Some(Token::DetectedSpinGY));
    assert_eq!(identify_token("detected_spin_g_z"),Some(Token::DetectedSpinGZ));
    assert_eq!(identify_token("detected_spin_multiplicity"),
        Some( Token::DetectedSpinMultiplicity));
    assert_eq!(identify_token("detected_spin_identity"),
        Some(Token::DetectedSpinIdentity));
    assert_eq!(identify_token("detected_spin_position"),
        Some(Token::DetectedSpinPosition));
    assert_eq!(identify_token("detected_spin_transition"),
        Some(Token::DetectedSpinTransition));
    assert_eq!(identify_token("diff"),Some( Token::Diff));
    assert_eq!(identify_token("distance"),Some( Token::Distance));
    assert_eq!(identify_token("\""),Some( Token::DoubleQuote));
    assert_eq!(identify_token("electric_quadrupole_coupling"),
        Some(Token::ElectricQuadrupoleCoupling));
    assert_eq!(identify_token("electric_quadrupole_x"),
        Some(Token::ElectricQuadrupoleX));
    assert_eq!(identify_token("electric_quadrupole_y"),
        Some(Token::ElectricQuadrupoleY));
    assert_eq!(identify_token("electric_quadrupole_z"),
        Some(Token::ElectricQuadrupoleZ));
    //assert_eq!(identify_token("element"),Some( Token::Element));
    assert_eq!(identify_token("\n"),Some( Token::EOL));
    assert_eq!(identify_token("="),Some( Token::Equals));
    assert_eq!(identify_token("exchange_groups_and_particles"),
        Some( Token::ExchangeGroupsAndParticles));
    assert_eq!(identify_token("extracells"), Some( Token::Extracells));
    assert_eq!(identify_token("false"),Some( Token::False));
    assert_eq!(identify_token("filter"),Some( Token::Filter));
    assert_eq!(identify_token("gcce"),Some( Token::GCCE));
    assert_eq!(identify_token("g_matrix"),Some( Token::GMatrix));
    assert_eq!(identify_token(">"),Some( Token::GreaterThan));
    assert_eq!(identify_token(">="),Some( Token::GreaterThanEqualTo));
    assert_eq!(identify_token("group"),Some( Token::Group));
    assert_eq!(identify_token("g_x"),Some( Token::GX));
    assert_eq!(identify_token("g_y"),Some( Token::GY));
    assert_eq!(identify_token("g_z"),Some( Token::GZ));
    assert_eq!(identify_token("hahn"),Some( Token::Hahn));
    assert_eq!(identify_token("^"),Some( Token::Hat));
    assert_eq!(identify_token("hyperfine_coupling"),
        Some(Token::HyperfineCoupling));
    assert_eq!(identify_token("hyperfine_x"),Some( Token::HyperfineX));
    assert_eq!(identify_token("hyperfine_y"),Some( Token::HyperfineY));
    assert_eq!(identify_token("hyperfine_z"),Some( Token::HyperfineZ));
    assert_eq!(identify_token("identity"),Some( Token::Identity));
    assert_eq!(identify_token("in"),Some( Token::In));
    assert_eq!(identify_token("input_structure_file"),
        Some(Token::InputStructureFile));
    assert_eq!(identify_token("isotope"),Some( Token::Isotope));
    assert_eq!(identify_token("isotope_abundances"),
        Some( Token::IsotopeAbundances));
    assert_eq!(identify_token("label"),Some( Token::Label));
    assert_eq!(identify_token("lebedev"),Some( Token::Lebedev));
    assert_eq!(identify_token("load_geometry"),Some( Token::LoadGeometry));
    assert_eq!(identify_token("<"),Some( Token::LessThan));
    assert_eq!(identify_token("<="),Some( Token::LessThanEqualTo));
    assert_eq!(identify_token("//"),Some( Token::LineComment));
    assert_eq!(identify_token("magnetic_field"),Some( Token::MagneticField));
    assert_eq!(identify_token("max_cell_size"),
        Some(Token::MaxCellSize));
    assert_eq!(identify_token("max_cluster_size"),
        Some(Token::MaxClusterSize));
    assert_eq!(identify_token("max_spin_order"),
        Some(Token::MaxSpinOrder));
    assert_eq!(identify_token("min_cell_size"),
        Some(Token::MinCellSize));
    assert_eq!(identify_token("-"),Some( Token::Minus));
    assert_eq!(identify_token("neighbor_cutoff_delta_hyperfine"),
        Some(Token::NeighborCutoffDeltaHyperfine));
    assert_eq!(identify_token("neighbor_cutoff_coupling"),
        Some(Token::NeighborCutoffCoupling));
    assert_eq!(identify_token("neighbor_cutoff_dipole_perpendicular"),
        Some(Token::NeighborCutoffDipolePerpendicular));
    assert_eq!(identify_token("neighbor_cutoff_distance"),
        Some(Token::NeighborCutoffDistance));
    assert_eq!(identify_token("neighbor_cutoff_3_spin_hahn_mod_depth"),
        Some(Token::NeighborCutoff3SpinHahnModDepth));
    assert_eq!(identify_token("neighbor_cutoff_3_spin_hahn_taylor_4"),
        Some(Token::NeighborCutoff3SpinHahnTaylor4));
    assert_eq!(identify_token("not"),Some( Token::Not));
    assert_eq!(identify_token("!="),Some( Token::NotEqual));
    assert_eq!(identify_token("not in"),Some( Token::NotIn));
    assert_eq!(identify_token("number_system_instances"),
        Some(Token::NumberSystemInstances));
    assert_eq!(identify_token("number_timepoints"),
        Some(Token::NumberTimepoints));
    assert_eq!(identify_token("orientation_grid"),Some(Token::OrientationGrid));
    assert_eq!(identify_token("particles"), Some( Token::Particles));
    assert_eq!(identify_token("partitioning_method"),
        Some( Token::PartitioningMethod));
    assert_eq!(identify_token("path"),Some( Token::Path));
    assert_eq!(identify_token(")"),Some( Token::ParenthesisClose));
    assert_eq!(identify_token("("),Some( Token::ParenthesisOpen));
    assert_eq!(identify_token("+"),Some( Token::Plus));
    assert_eq!(identify_token("primary_cell"), Some( Token::PrimaryCell));
    assert_eq!(identify_token("pulse_sequence"),Some( Token::PulseSequence));
    assert_eq!(identify_token("radius"),Some( Token::Radius));
    assert_eq!(identify_token("random"),Some( Token::Random));
    assert_eq!(identify_token("read_csv"), Some(Token::ReadCSV));
    assert_eq!(identify_token("rng_seed"),Some( Token::RNGSeed));
    assert_eq!(identify_token("save_dir"),Some( Token::SaveDir));
    assert_eq!(identify_token(";"),Some( Token::Semicolon));
    assert_eq!(identify_token("set"),Some( Token::Set));
    assert_eq!(identify_token("#"),Some( Token::Sharp));
    assert_eq!(identify_token("/"),Some( Token::Slash));
    //assert_eq!(identify_token("sphere"),Some( Token::Sphere));
    assert_eq!(identify_token("spin_properties"),
        Some(Token::SpinProperties));
    assert_eq!(identify_token("spin"),Some( Token::Spin));
    assert_eq!(identify_token("]"),Some( Token::SquareBracketClose));
    assert_eq!(identify_token("["),Some( Token::SquareBracketOpen));
    assert_eq!(identify_token("structure_properties"),
        Some(Token::StructureProperties));
    assert_eq!(identify_token("structures"),Some( Token::Structures));
    assert_eq!(identify_token("system_name"),Some( Token::SystemName));
    assert_eq!(identify_token("temperature"),Some( Token::Temperature));
    assert_eq!(identify_token("tensors"),Some( Token::Tensors));
    assert_eq!(identify_token("time_increments"),
        Some(Token::TimeIncrements));
    assert_eq!(identify_token("thermal"),Some( Token::Thermal));
    assert_eq!(identify_token("*"),Some( Token::Times));
    assert_eq!(identify_token("true"),Some( Token::True));
    assert_eq!(identify_token("tunnel_splitting"),
        Some(Token::TunnelSplitting));
    assert_eq!(identify_token("type"),Some( Token::Type));
    assert_eq!(identify_token("unit_of_clustering"),
        Some( Token::UnitOfClustering));
    assert_eq!(identify_token("vector"),Some( Token::Vector));
    assert_eq!(identify_token("void_probability"),
        Some( Token::VoidProbability));
    assert_eq!(identify_token(" "),Some( Token::Whitespace));
    assert_eq!(identify_token("write_auxiliary_signals"),
        Some(Token::WriteAuxiliarySignals));
    assert_eq!(identify_token("write_bath"),
        Some(Token::WriteBath));
    assert_eq!(identify_token("write_clusters"),
        Some(Token::WriteClusters));
    assert_eq!(identify_token("write_exchange_groups"),
        Some(Token::WriteExchangeGroups));
    assert_eq!(identify_token("write_info"),Some(Token::WriteInfo));
    assert_eq!(identify_token("write_methyl_partitions"),
        Some(Token::WriteMethylPartitions));
    assert_eq!(identify_token("write_orientation_signals"),
        Some(Token::WriteOrientationSignals));
    assert_eq!(identify_token("write_sans_spin_signals"),
        Some(Token::WriteSansSpinSignals));
    assert_eq!(identify_token("write_structure_pdb"),
        Some(Token::WriteStructurePDB));
    assert_eq!(identify_token("write_tensors"),
        Some(Token::WriteTensors));
    assert_eq!(identify_token("indices"),Some( Token::Indices)); 
    assert_eq!(identify_token("elements"),Some( Token::Elements)); 
    assert_eq!(identify_token("serials"),Some( Token::Serials)); 
    assert_eq!(identify_token("set"),Some( Token::Set)); 
    assert_eq!(identify_token("r2cce"),Some( Token::R2CCE));
    assert_eq!(identify_token("residues"),Some( Token::Residues));
    assert_eq!(identify_token("residue_sequence_numbers"),
        Some(Token::ResSeqNums));
    assert_eq!(identify_token("bonded_indices"),Some( Token::BondedIndices));
    assert_eq!(identify_token("bonded_elements"),Some(Token::BondedElements));

  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


