use crate::clue_errors::*;
use crate::physical_constants::Isotope;
use std::ops::{Add,Sub,Mul,Div};

use std::fmt;
//------------------------------------------------------------------------------
#[derive(PartialEq, Debug, Clone)]
pub enum Token{
 Bang,
 BondedElements,                                               
 BondedIndices,                                                     
 BlockCommentEnd, 
 BlockCommentStart, 
 CCE,
 CarrPurcell,
 CentroidOverSerials,
 ClusterBatchSize,
 ClusterMethod,
 Clusters,
 Comma,
 Config,
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
 DoubleQuote,
 ElectricQuadrupoleCoupling,
 ElectricQuadrupoleX,
 ElectricQuadrupoleY,
 ElectricQuadrupoleZ,
 Element,
 Elements,                                                    
 EOL,
 Equals,
 False,
 Filter,
 Float(f64),
 GreaterThan,
 GreaterThanEqualTo,
 Hat,
 Hahn,
 HyperfineCoupling,
 HyperfineX,
 HyperfineY,
 HyperfineZ,
 In,
 Indices,                                                      
 InputStructureFile,
 Int(i32),
 Isotope,
 Label,
 Lebedev,
 LessThan,
 LessThanEqualTo,
 LineComment, 
 LoadGeometry,
 MagneticField,
 MaxClusterSize,
 Minus,
 Mode(ModeAttribute),
 NeighborCutoffDeltaHyperfine,
 NeighborCutoffDipoleDipole,
 NeighborCutoffDistance,
 NeighborCutoff3SpinHahnModDepth,
 NeighborCutoff3SpinHahnTaylor4,
 Not,
 NotEqual,
 NotIn,
 NumberTimepoints,
 Path,
 ParenthesisClose,
 ParenthesisOpen,
 Plus,
 PulseSequence,
 R2CCE,
 Radius,
 Random,
 ReadCSV,
 RemovePartialMethyls,
 Residue,
 Residues,                                                     
 ResSeqNums,                                                   
 RNGSeed,
 Semicolon,
 Serials,                                                      
 Sharp,
 Slash,
 Sphere,
 SpinProperties,
 Spins,
 SquareBracketClose,
 SquareBracketOpen,
 StructureProperties,
 Structures,
 Temperature,
 Tensors,
 TimeIncrements,
 Times,
 True,
 TunnelSplitting,
 Type,
 UserInputValue(String),
 VectorF64(Vec::<f64>),
 VectorI32(Vec::<i32>),
 VectorString(Vec::<String>),
 Whitespace,
 WriteAuxiliarySignals,
 WriteBath,
 WriteClusters,
 WriteInfo,
 WriteExchangeGroups,
 WriteOrientationSignals,
 WriteStructurePDB,
 WriteTensors,
}
impl fmt::Display for Token{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{
      Token::Bang => write!(f,"!"), 
      Token::BondedIndices => write!(f,"bonded_indices"),
      Token::BondedElements => write!(f,"bonded_elements"),
      Token::BlockCommentEnd => write!(f,"*/"), 
      Token::BlockCommentStart => write!(f,"/*"), 
      Token::CarrPurcell => write!(f,"cp"),
      Token::CCE => write!(f,"cce"),
      Token::CentroidOverSerials => write!(f,"centroid_over_serials"),
      Token::ClusterBatchSize => write!(f,"cluster_batch_size"),
      Token::ClusterMethod => write!(f,"cluster_method"),
      Token::Clusters => write!(f,"clusters"),
      Token::Comma => write!(f,","),
      Token::Config => write!(f,"config"),
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
      Token::DoubleQuote => write!(f,"\""),
      Token::ElectricQuadrupoleCoupling  
        => write!(f,"electric_quadrupole_coupling"),
      Token::ElectricQuadrupoleX  => write!(f,"electric_quadrupole_x"),
      Token::ElectricQuadrupoleY  => write!(f,"electric_quadrupole_y"),
      Token::ElectricQuadrupoleZ  => write!(f,"electric_quadrupole_z"),
      Token::Element => write!(f,"element"),
      Token::Elements => write!(f,"elements"),
      Token::EOL => writeln!(f),
      Token::Equals => write!(f,"="),
      Token::False => writeln!(f,"false"),
      Token::Filter => write!(f,"filter"),
      Token::Float(x) => write!(f,"{}",x),
      Token::GreaterThan => write!(f,">"),
      Token::GreaterThanEqualTo => write!(f,">="),
      Token::Hahn => write!(f,"hahn"),
      Token::Hat => write!(f,"^"),
      Token::HyperfineCoupling => write!(f,"hyperfine_coupling"),
      Token::HyperfineX => write!(f,"hyperfine_x"),
      Token::HyperfineY => write!(f,"hyperfine_y"),
      Token::HyperfineZ => write!(f,"hyperfine_z"),
      Token::In => write!(f,"in"),
      Token::Indices => write!(f,"indices"),
      Token::InputStructureFile => write!(f,"input_structure_file"),
      Token::Int(n) => write!(f,"{}",n),
      Token::Isotope => write!(f,"isotope"),
      Token::Label => write!(f,"label"),
      Token::Lebedev => write!(f,"lebedev"),
      Token::LessThan => write!(f,"<"),
      Token::LessThanEqualTo => write!(f,"<="),
      Token::LineComment => write!(f,"//"), 
      Token::LoadGeometry => write!(f,"load_geometry"),
      Token::MagneticField => write!(f,"magnetic_field"),
      Token::MaxClusterSize => write!(f,"max_cluster_size"),
      Token::Minus => write!(f,"-"),
      Token::Mode(mode) => write!(f,"{}",mode),
      Token::NeighborCutoffDeltaHyperfine 
        => write!(f,"neighbor_cutoff_delta_hyperfine"),
      Token::NeighborCutoffDipoleDipole 
        => write!(f,"neighbor_cutoff_dipole_dipole"),
      Token::NeighborCutoffDistance 
        => write!(f,"neighbor_cutoff_distance"),
      Token::NeighborCutoff3SpinHahnModDepth 
        => write!(f,"neighbor_cutoff_3_spin_hahn_mod_depth"),
      Token::NeighborCutoff3SpinHahnTaylor4 
        => write!(f,"neighbor_cutoff_3_spin_hahn_taylor_4"),
      Token::Not => write!(f,"not"),
      Token::NotEqual => write!(f,"!="),
      Token::NotIn => write!(f,"not in"),
      Token::NumberTimepoints => write!(f,"*"),
      Token::Path => write!(f,"path"),
      Token::ParenthesisClose => write!(f,")"),
      Token::ParenthesisOpen => write!(f,"("),
      Token::Plus => write!(f,"+"),
      Token::PulseSequence => write!(f,"pulse_sequence"),
      Token::R2CCE => write!(f,"r2cce"),
      Token::Radius => write!(f,"radius"),
      Token::Random => write!(f,"random"),
      Token::ReadCSV => write!(f,"read_csv"),
      Token::RemovePartialMethyls => write!(f,"remove_partial_methyls"),
      Token::Residue => write!(f,"residue"),
      Token::Residues => write!(f,"residues"),
      Token::ResSeqNums => write!(f,"residue_sequence_numbers"),
      Token::RNGSeed => write!(f,"rng_seed"),
      Token::Semicolon => write!(f,";"),
      Token::Serials => write!(f,"serials"),
      Token::Sharp => write!(f,"#"),
      Token::Slash => write!(f,"/"),
      Token::SpinProperties => write!(f,"spin_properties"),
      Token::Spins => write!(f,"spins"),
      Token::Sphere => write!(f,"sphere"),
      Token::SquareBracketClose => write!(f,"]"),
      Token::SquareBracketOpen => write!(f,"["),
      Token::StructureProperties => write!(f,"structure_properties"),
      Token::Structures => write!(f,"structures"),
      Token::Temperature => write!(f,"temperature"),
      Token::Tensors => write!(f,"tensors"),
      Token::TimeIncrements => write!(f,"time_increments"),
      Token::Times => write!(f,"*"),
      Token::True => writeln!(f,"true"),
      Token::TunnelSplitting => write!(f,"tunnel_splitting"),
      Token::Type => write!(f,"type"),
      Token::UserInputValue(string) => write!(f,"{}",string),
      Token::VectorF64(v) => write!(f,"{:?}",v), 
      Token::VectorI32(v) => write!(f,"{:?}",v), 
      Token::VectorString(v) => write!(f,"{:?}",v), 
      Token::Whitespace => write!(f," "),
      Token::WriteAuxiliarySignals => write!(f,"write_auxiliary_signals"),
      Token::WriteBath => write!(f,"write_bath"),
      Token::WriteClusters => write!(f,"write_clusters"),
      Token::WriteExchangeGroups => write!(f,"write_exchange_groups"),
      Token::WriteInfo => write!(f,"write_info"),
      Token::WriteOrientationSignals => write!(f,"write_orientation_signals"),
      Token::WriteStructurePDB => write!(f,"write_structure_pdb"),
      Token::WriteTensors => write!(f,"write_tensors"),
    }
  }
}
//------------------------------------------------------------------------------
pub fn identify_token(word: &str) -> Option<Token>{
  match word{
    "!" => Some(Token::Bang),
    "bonded_indices" => Some(Token::BondedIndices),
    "bonded_elements" => Some(Token::BondedElements),
    "*/" => Some(Token::BlockCommentEnd),
    "/*" => Some(Token::BlockCommentStart),
    "cp" => Some(Token::CarrPurcell),
    "cce" => Some(Token::CCE),
    "centroid_over_serials" => Some(Token::CentroidOverSerials),
    "cluster_batch_size" => Some(Token::ClusterBatchSize),
    "cluster_method" => Some(Token::ClusterMethod),
    "clusters" => Some(Token::Clusters),
    "," => Some(Token::Comma),
    "config" => Some(Token::Config),
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
    "\"" => Some(Token::DoubleQuote),
    "electric_quadrupole_coupling" => Some(Token::ElectricQuadrupoleCoupling),
    "electric_quadrupole_x" => Some(Token::ElectricQuadrupoleX),
    "electric_quadrupole_y" => Some(Token::ElectricQuadrupoleY),
    "electric_quadrupole_z" => Some(Token::ElectricQuadrupoleZ),
    "element" => Some(Token::Element),
    "elements" => Some(Token::Elements), 
    "\n" => Some(Token::EOL),
    "=" => Some(Token::Equals),
    "false" => Some(Token::False),
    "filter" => Some(Token::Filter),
    ">" => Some(Token::GreaterThan),
    ">=" => Some(Token::GreaterThanEqualTo),
    "in" => Some(Token::In),
    "input_structure_file" => Some(Token::InputStructureFile),
    "indices" => Some(Token::Indices), 
    "isotope" => Some(Token::Isotope),
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
    "max_cluster_size" => Some(Token::MaxClusterSize),
    "-" => Some(Token::Minus),
    "neighbor_cutoff_delta_hyperfine" 
      => Some(Token::NeighborCutoffDeltaHyperfine),
    "neighbor_cutoff_dipole_dipole" 
      => Some(Token::NeighborCutoffDipoleDipole),
    "neighbor_cutoff_distance" 
      => Some(Token::NeighborCutoffDistance),
    "neighbor_cutoff_3_spin_hahn_mod_depth" 
      => Some(Token::NeighborCutoff3SpinHahnModDepth),
    "neighbor_cutoff_3_spin_hahn_taylor_4" 
      => Some(Token::NeighborCutoff3SpinHahnTaylor4),
    "not" => Some(Token::Not),
    "!=" => Some(Token::NotEqual),
    "not in" => Some(Token::NotIn),
    "number_timepoints" => Some(Token::NumberTimepoints),
    "path" => Some(Token::Path),
    ")" => Some(Token::ParenthesisClose),
    "(" => Some(Token::ParenthesisOpen),
    "+" => Some(Token::Plus),
    "pulse_sequence" => Some(Token::PulseSequence),
    "r2cce" => Some(Token::R2CCE),
    "radius" => Some(Token::Radius),
    "random" => Some(Token::Random),
    "read_csv" => Some(Token::ReadCSV),
    "remove_partial_methyls" => Some(Token::RemovePartialMethyls),
    "residue" => Some(Token::Residue),
    "residues" => Some(Token::Residues),
    "residue_sequence_numbers" => Some(Token::ResSeqNums),
    "rng_seed" => Some(Token::RNGSeed),
    ";" => Some(Token::Semicolon),
    "serials" => Some(Token::Serials), 
    "#" => Some(Token::Sharp),
    "/" => Some(Token::Slash),
    "spin_properties" => Some(Token::SpinProperties),
    "sphere" => Some(Token::Sphere),
    "spins" => Some(Token::Spins),
    "]" => Some(Token::SquareBracketClose),
    "[" => Some(Token::SquareBracketOpen),
    "structure_properties" => Some(Token::StructureProperties),
    "structures" => Some(Token::Structures),
    "temperature" => Some(Token::Temperature),
    "tensors" => Some(Token::Tensors),
    "*" => Some(Token::Times),
    "true" => Some(Token::True),
    "time_increments" => Some(Token::TimeIncrements),
    "tunnel_splitting" => Some(Token::TunnelSplitting),
    "type" => Some(Token::Type),
    " " => Some(Token::Whitespace),
    "write_auxiliary_signals" => Some(Token::WriteAuxiliarySignals),
    "write_bath" => Some(Token::WriteBath),
    "write_clusters" => Some(Token::WriteClusters),
    "write_exchange_groups" => Some(Token::WriteExchangeGroups),
    "write_info" => Some(Token::WriteInfo),
    "write_orientation_signals" => Some(Token::WriteOrientationSignals),
    "write_structure_pdb" => Some(Token::WriteStructurePDB),
    "write_tensors" => Some(Token::WriteTensors),
    _ => None
  }
}

impl Add for Token{
  
  type Output = Result<Token,CluEError>;
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

//trait Pow{
//  fn pow(self, )
//}
//impl Pow for Token{

//  type Output = Result<Token,CluEError>;
impl Token{
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
pub fn is_relational_operator(token: &Token) -> bool{
  matches!( token,
    Token::Equals | Token::GreaterThan | Token::GreaterThanEqualTo
    | Token::In | Token::LessThan | Token::LessThanEqualTo | Token::NotEqual 
    | Token::NotIn
    )
}
//------------------------------------------------------------------------------
pub fn is_opening_delimiter(token: &Token) -> bool{
  matches!(token,
    Token::ParenthesisOpen | Token::SquareBracketOpen | Token::CurlyBracketOpen
    )
}
//------------------------------------------------------------------------------
pub fn is_closing_delimiter(token: &Token) -> bool{
  matches!(token,
    Token::ParenthesisClose | Token::SquareBracketClose 
    | Token::CurlyBracketClose
    )
}
//------------------------------------------------------------------------------
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
#[derive(PartialEq, Debug, Clone)]
pub struct ModeAttribute{
  pub mode:  ConfigMode,
  pub isotope: Option<Isotope>,
  pub label: Option<String>,
  pub path: Option<String>,
  //pub mode_type: Option<String>,
}

impl Default for ModeAttribute{
  fn default() -> Self{
    ModeAttribute{
      mode:  ConfigMode::Config,
      isotope: None,
      label: None,
      path: None,
      //mode_type: None,
    }
  }
}

impl ModeAttribute{

  pub fn new() -> Self{
    ModeAttribute::default()
  }

  pub fn from(tokens: Vec::<Token>) -> Result<Self,CluEError>
  {

    let sharp_indices = find_token(&Token::Sharp,&tokens);
    if sharp_indices.len() != 1 || sharp_indices[0] != 0{
      return Err(CluEError::ModeAttributeWrongSharp);
    } 

    let idx = sharp_indices[0];
    if tokens[idx+1] != Token::SquareBracketOpen
    && tokens[tokens.len()-1] != Token::SquareBracketClose{
      return Err(CluEError::ModeAttributeWrongBrackets);
    }

    let mode = ConfigMode::from(tokens[idx+2].clone())?; 

    let isotope: Option<Isotope>;
    let isotope_indices = find_token(&Token::Isotope, &tokens);

    if isotope_indices.is_empty(){
      isotope = None;
    } else if isotope_indices.len() == 1 
      && tokens[isotope_indices[0]+1] == Token::Equals{
       if  let Token::UserInputValue(value)=&tokens[isotope_indices[0]+2]{
        let isotope_value = Isotope::from(value)?;
        isotope = Some(isotope_value);
      }else{
        return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
      }
       
    } else {
      return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
    }

    let label: Option<String>;
    let label_indices = find_token(&Token::Label,&tokens);
    
    if label_indices.is_empty(){
      label = None;
    } else if label_indices.len() == 1 
      && tokens[label_indices[0]+1] == Token::Equals{
      if let Token::UserInputValue(value)=&tokens[label_indices[0]+2]{
        label = Some(value.clone());
      }else{
        return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
      }
    } else {
      return Err(CluEError::ModeAttributeWrongOption(mode.to_string()));
    }

    let path: Option<String>;
    let path_indices = find_token(&Token::Path, &tokens);

    if path_indices.is_empty(){
      path = None;
    } else if path_indices.len() == 1 
      && tokens[path_indices[0]+1] == Token::Equals{
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

#[derive(PartialEq, Debug, Clone)]
pub enum ConfigMode{
  Clusters,
  Config,
  Filter,
  SpinProperties,
  Spins,
  StructureProperties,
  Structures,
  Tensors,
}

impl ConfigMode{
  pub fn from(token: Token) -> Result<Self,CluEError>{

    match token{
      Token::Clusters => Ok(ConfigMode::Clusters),
      Token::Config => Ok(ConfigMode::Config),
      Token::Filter => Ok(ConfigMode::Filter),
      Token::SpinProperties => Ok(ConfigMode::SpinProperties),
      Token::Spins => Ok(ConfigMode::Spins),
      Token::StructureProperties => Ok(ConfigMode::StructureProperties),
      Token::Structures => Ok(ConfigMode::Structures),
      Token::Tensors => Ok(ConfigMode::Tensors),
      _ => Err(CluEError::ConfigModeNotRecognized(token.to_string())),
    }
  }
  //----------------------------------------------------------------------------
}
impl fmt::Display for ConfigMode{
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{
      ConfigMode::Clusters => write!(f,"clusters"),
      ConfigMode::Config => write!(f,"config"),
      ConfigMode::Filter => write!(f,"filter"),
      ConfigMode::SpinProperties => write!(f,"spin_properties"),
      ConfigMode::Spins => write!(f,"spins"),
      ConfigMode::StructureProperties => write!(f,"structure_properties"),
      ConfigMode::Structures => write!(f,"structures"),
      ConfigMode::Tensors => write!(f,"tensors"),
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
    assert_eq!(identify_token("!"), Some(Token::Bang));
    assert_eq!(identify_token("*/"), Some(Token::BlockCommentEnd));
    assert_eq!(identify_token("/*"),Some( Token::BlockCommentStart));
    assert_eq!(identify_token("cp"),Some( Token::CarrPurcell));
    assert_eq!(identify_token("cce"),Some( Token::CCE));
    assert_eq!(identify_token("centroid_over_serials"),
        Some(Token::CentroidOverSerials));
    assert_eq!(identify_token("cluster_batch_size"),
        Some(Token::ClusterBatchSize));
    assert_eq!(identify_token("cluster_method"),Some(Token::ClusterMethod));
    assert_eq!(identify_token("clusters"),Some( Token::Clusters));
    assert_eq!(identify_token(","),Some( Token::Comma));
    assert_eq!(identify_token("config"),Some( Token::Config));
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
    assert_eq!(identify_token("\""),Some( Token::DoubleQuote));
    assert_eq!(identify_token("electric_quadrupole_coupling"),
        Some(Token::ElectricQuadrupoleCoupling));
    assert_eq!(identify_token("electric_quadrupole_x"),
        Some(Token::ElectricQuadrupoleX));
    assert_eq!(identify_token("electric_quadrupole_y"),
        Some(Token::ElectricQuadrupoleY));
    assert_eq!(identify_token("electric_quadrupole_z"),
        Some(Token::ElectricQuadrupoleZ));
    assert_eq!(identify_token("element"),Some( Token::Element));
    assert_eq!(identify_token("\n"),Some( Token::EOL));
    assert_eq!(identify_token("="),Some( Token::Equals));
    assert_eq!(identify_token("false"),Some( Token::False));
    assert_eq!(identify_token("filter"),Some( Token::Filter));
    assert_eq!(identify_token(">"),Some( Token::GreaterThan));
    assert_eq!(identify_token(">="),Some( Token::GreaterThanEqualTo));
    assert_eq!(identify_token("hahn"),Some( Token::Hahn));
    assert_eq!(identify_token("^"),Some( Token::Hat));
    assert_eq!(identify_token("hyperfine_coupling"),
        Some(Token::HyperfineCoupling));
    assert_eq!(identify_token("hyperfine_x"),Some( Token::HyperfineX));
    assert_eq!(identify_token("hyperfine_y"),Some( Token::HyperfineY));
    assert_eq!(identify_token("hyperfine_z"),Some( Token::HyperfineZ));
    assert_eq!(identify_token("in"),Some( Token::In));
    assert_eq!(identify_token("input_structure_file"),
        Some(Token::InputStructureFile));
    assert_eq!(identify_token("isotope"),Some( Token::Isotope));
    assert_eq!(identify_token("label"),Some( Token::Label));
    assert_eq!(identify_token("lebedev"),Some( Token::Lebedev));
    assert_eq!(identify_token("load_geometry"),Some( Token::LoadGeometry));
    assert_eq!(identify_token("<"),Some( Token::LessThan));
    assert_eq!(identify_token("<="),Some( Token::LessThanEqualTo));
    assert_eq!(identify_token("//"),Some( Token::LineComment));
    assert_eq!(identify_token("magnetic_field"),Some( Token::MagneticField));
    assert_eq!(identify_token("max_cluster_size"),
        Some(Token::MaxClusterSize));
    assert_eq!(identify_token("-"),Some( Token::Minus));
    assert_eq!(identify_token("neighbor_cutoff_delta_hyperfine"),
        Some(Token::NeighborCutoffDeltaHyperfine));
    assert_eq!(identify_token("neighbor_cutoff_dipole_dipole"),
        Some(Token::NeighborCutoffDipoleDipole));
    assert_eq!(identify_token("neighbor_cutoff_distance"),
        Some(Token::NeighborCutoffDistance));
    assert_eq!(identify_token("neighbor_cutoff_3_spin_hahn_mod_depth"),
        Some(Token::NeighborCutoff3SpinHahnModDepth));
    assert_eq!(identify_token("neighbor_cutoff_3_spin_hahn_taylor_4"),
        Some(Token::NeighborCutoff3SpinHahnTaylor4));
    assert_eq!(identify_token("not"),Some( Token::Not));
    assert_eq!(identify_token("!="),Some( Token::NotEqual));
    assert_eq!(identify_token("not in"),Some( Token::NotIn));
    assert_eq!(identify_token("number_timepoints"),
        Some(Token::NumberTimepoints));
    assert_eq!(identify_token("path"),Some( Token::Path));
    assert_eq!(identify_token(")"),Some( Token::ParenthesisClose));
    assert_eq!(identify_token("("),Some( Token::ParenthesisOpen));
    assert_eq!(identify_token("+"),Some( Token::Plus));
    assert_eq!(identify_token("pulse_sequence"),Some( Token::PulseSequence));
    assert_eq!(identify_token("radius"),Some( Token::Radius));
    assert_eq!(identify_token("random"),Some( Token::Random));
    assert_eq!(identify_token("read_csv"), Some(Token::ReadCSV));
    assert_eq!(identify_token("remove_partial_methyls"),
        Some(Token::RemovePartialMethyls));
    assert_eq!(identify_token("residue"),Some( Token::Residue));
    assert_eq!(identify_token("rng_seed"),Some( Token::RNGSeed));
    assert_eq!(identify_token(";"),Some( Token::Semicolon));
    assert_eq!(identify_token("#"),Some( Token::Sharp));
    assert_eq!(identify_token("/"),Some( Token::Slash));
    assert_eq!(identify_token("sphere"),Some( Token::Sphere));
    assert_eq!(identify_token("spin_properties"),
        Some(Token::SpinProperties));
    assert_eq!(identify_token("spins"),Some( Token::Spins));
    assert_eq!(identify_token("]"),Some( Token::SquareBracketClose));
    assert_eq!(identify_token("["),Some( Token::SquareBracketOpen));
    assert_eq!(identify_token("structure_properties"),
        Some(Token::StructureProperties));
    assert_eq!(identify_token("structures"),Some( Token::Structures));
    assert_eq!(identify_token("temperature"),Some( Token::Temperature));
    assert_eq!(identify_token("tensors"),Some( Token::Tensors));
    assert_eq!(identify_token("time_increments"),
        Some(Token::TimeIncrements));
    assert_eq!(identify_token("*"),Some( Token::Times));
    assert_eq!(identify_token("true"),Some( Token::True));
    assert_eq!(identify_token("tunnel_splitting"),
        Some(Token::TunnelSplitting));
    assert_eq!(identify_token("type"),Some( Token::Type));
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
    assert_eq!(identify_token("write_orientation_signals"),
        Some(Token::WriteOrientationSignals));
    assert_eq!(identify_token("write_structure_pdb"),
        Some(Token::WriteStructurePDB));
    assert_eq!(identify_token("write_tensors"),
        Some(Token::WriteTensors));
    assert_eq!(identify_token("indices"),Some( Token::Indices)); 
    assert_eq!(identify_token("elements"),Some( Token::Elements)); 
    assert_eq!(identify_token("serials"),Some( Token::Serials)); 
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


