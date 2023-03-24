use crate::clue_errors::*;
use std::ops::{Add,Sub,Mul,Div};
//------------------------------------------------------------------------------
#[derive(PartialEq, Debug, Clone)]
pub enum Token{
 Bang,
 BlockCommentEnd, 
 BlockCommentStart, 
 Clusters,
 Comma,
 Config,
 CurlyBracketClose,
 CurlyBracketOpen,
 Element,
 EOL,
 Equals,
 Filter,
 Float(f64),
 GreaterThan,
 GreaterThanEqualTo,
 Hat,
 In,
 Int(i32),
 Label,
 LessThan,
 LessThanEqualTo,
 LineComment, 
 MagneticField,
 Minus,
 Mode(ModeAttribute),
 Not,
 NotEqual,
 NotIn,
 Path,
 ParenthesisClose,
 ParenthesisOpen,
 Plus,
 Radius,
 Residue,
 RNGSeed,
 Semicolon,
 Sharp,
 Slash,
 Spins,
 SquareBracketClose,
 SquareBracketOpen,
 Structures,
 Tensors,
 Times,
 TunnelSplitting,
 Type,
 UserInputValue(String),
 VectorF64(Vec::<f64>),
 VectorI32(Vec::<i32>),
 VectorString(Vec::<String>),
 Whitespace,
 WriteStructurePDB,
 Indices,                                                      
 //NotIndices,                                                   
 Elements,                                                    
 //NotElements,                                                  
 Serials,                                                      
 //NotSerials,                                                   
 Residues,                                                     
 //NotResidues,                                                  
 ResSeqNums,                                                   
 //NotResSeqNums,                                                
 BondedIndices,                                                     
 //NotBondedIndices,                                                  
 BondedElements,                                               
 //NotBondedElements,
}
impl Token{
  pub fn to_string(&self) -> String{
    match self{
      Token::Bang => "!".to_string(), 
      Token::BlockCommentEnd => "*/".to_string(), 
      Token::BlockCommentStart => "/*".to_string(), 
      Token::Clusters => "clusters".to_string(),
      Token::Comma => ",".to_string(),
      Token::Config => "config".to_string(),
      Token::CurlyBracketClose => "}".to_string(),
      Token::CurlyBracketOpen => "{".to_string(),
      Token::Element => "element".to_string(),
      Token::EOL => "\n".to_string(),
      Token::Equals => "=".to_string(),
      Token::Filter => "filter".to_string(),
      Token::Float(x) => x.to_string(),
      Token::GreaterThan => ">".to_string(),
      Token::GreaterThanEqualTo => ">=".to_string(),
      Token::Hat => "^".to_string(),
      Token::In => "in".to_string(),
      Token::Int(n) => n.to_string(),
      Token::Label => "label".to_string(),
      Token::LessThan => "<".to_string(),
      Token::LessThanEqualTo => "<=".to_string(),
      Token::LineComment => "//".to_string(), 
      Token::MagneticField => "magnetic_field".to_string(),
      Token::Minus => "-".to_string(),
      Token::Mode(mode) => mode.to_string(),
      Token::Not => "not".to_string(),
      Token::NotEqual => "!=".to_string(),
      Token::NotIn => "not in".to_string(),
      Token::Path => "path".to_string(),
      Token::ParenthesisClose => ")".to_string(),
      Token::ParenthesisOpen => "(".to_string(),
      Token::Plus => "+".to_string(),
      Token::Radius => "radius".to_string(),
      Token::Residue => "residue".to_string(),
      Token::RNGSeed => "rng_seed".to_string(),
      Token::Semicolon => ";".to_string(),
      Token::Sharp => "#".to_string(),
      Token::Slash => "/".to_string(),
      Token::Spins => "spins".to_string(),
      Token::SquareBracketClose => "]".to_string(),
      Token::SquareBracketOpen => "[".to_string(),
      Token::Structures => "structures".to_string(),
      Token::Tensors => "tensors".to_string(),
      Token::Times => "*".to_string(),
      Token::TunnelSplitting => "tunnel_splitting".to_string(),
      Token::Type => "type".to_string(),
      Token::UserInputValue(string) => (*string).clone(),
      Token::VectorF64(v) => format!("{:?}",v), 
      Token::VectorI32(v) => format!("{:?}",v), 
      Token::VectorString(v) => format!("{:?}",v), 
      Token::Whitespace => " ".to_string(),
      Token::WriteStructurePDB => "write_structure_pdb".to_string(),
      Token::Indices => "indices".to_string(),
      Token::Elements => "elements".to_string(),
      Token::Serials => "serials".to_string(),
      Token::Residues => "residues".to_string(),
      Token::ResSeqNums => "residue_sequence_numbers".to_string(),
      Token::BondedIndices => "bonded_indices".to_string(),
      Token::BondedElements => "bonded_elements".to_string(),
    }
  }
}
//------------------------------------------------------------------------------
pub fn identify_token(word: &str) -> Option<Token>{
  match word{
    "!" => Some(Token::Bang),
    "*/" => Some(Token::BlockCommentEnd),
    "/*" => Some(Token::BlockCommentStart),
    "clusters" => Some(Token::Clusters),
    "," => Some(Token::Comma),
    "config" => Some(Token::Config),
    "}" => Some(Token::CurlyBracketClose),
    "{" => Some(Token::CurlyBracketOpen),
    "element" => Some(Token::Element),
    "\n" => Some(Token::EOL),
    "=" => Some(Token::Equals),
    "filter" => Some(Token::Filter),
    ">" => Some(Token::GreaterThan),
    ">=" => Some(Token::GreaterThanEqualTo),
    "in" => Some(Token::In),
    "label" => Some(Token::Label),
    "^" => Some(Token::Hat),
    "<" => Some(Token::LessThan),
    "<=" => Some(Token::LessThanEqualTo),
    "//" => Some(Token::LineComment),
    "magnetic_field" => Some(Token::MagneticField),
    "-" => Some(Token::Minus),
    "not" => Some(Token::Not),
    "!=" => Some(Token::NotEqual),
    "not in" => Some(Token::NotIn),
    "path" => Some(Token::Path),
    ")" => Some(Token::ParenthesisClose),
    "(" => Some(Token::ParenthesisOpen),
    "+" => Some(Token::Plus),
    "radius" => Some(Token::Radius),
    "residue" => Some(Token::Residue),
    "rng_seed" => Some(Token::RNGSeed),
    ";" => Some(Token::Semicolon),
    "#" => Some(Token::Sharp),
    "/" => Some(Token::Slash),
    "spins" => Some(Token::Spins),
    "]" => Some(Token::SquareBracketClose),
    "[" => Some(Token::SquareBracketOpen),
    "structures" => Some(Token::Structures),
    "tensors" => Some(Token::Tensors),
    "*" => Some(Token::Times),
    "tunnel_splitting" => Some(Token::TunnelSplitting),
    "type" => Some(Token::Type),
    " " => Some(Token::Whitespace),
    "write_structure_pdb" => Some(Token::WriteStructurePDB),
    "indices" => Some(Token::Indices), 
    "elements" => Some(Token::Elements), 
    "serials" => Some(Token::Serials), 
    "residues" => Some(Token::Residues),
    "residue_sequence_numbers" => Some(Token::ResSeqNums),
    "bonded_indices" => Some(Token::BondedIndices),
    "bonded_elements" => Some(Token::BondedElements),
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
        Ok(Token::Int(a.pow(b.try_into().unwrap() ) )) 
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
                |b| a.pow((*b).try_into().unwrap()) ).collect()  ))
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
                |a| a.pow(b.try_into().unwrap()) ).collect()  ))
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
          w.push(u[ii].pow(v[ii].try_into().unwrap()) ); 
        }
        Ok(Token::VectorI32(w))
      },
      (_,_) => Err(CluEError::CannotPowTokens),

    }
  }
}
//------------------------------------------------------------------------------
pub fn is_relational_operator(token: &Token) -> bool{
  match token{
    Token::Equals => true,
    Token::GreaterThan => true,
    Token::GreaterThanEqualTo => true,
    Token::In => true,
    Token::LessThan => true,
    Token::LessThanEqualTo => true,
    Token::NotEqual => true,
    Token::NotIn => true,
    _ => false,
  }
}
//------------------------------------------------------------------------------
pub fn is_opening_delimiter(token: &Token) -> bool{
  match token{
    Token::ParenthesisOpen => true,
    Token::SquareBracketOpen => true,
    Token::CurlyBracketOpen => true,
    _ => false,
  }
}
//------------------------------------------------------------------------------
pub fn is_closing_delimiter(token: &Token) -> bool{
  match token{
    Token::ParenthesisClose => true,
    Token::SquareBracketClose => true,
    Token::CurlyBracketClose => true,
    _ => false,
  }
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
  pub label: Option<String>,
  pub path: Option<String>,
  //pub mode_type: Option<String>,
}

impl Default for ModeAttribute{
  fn default() -> Self{
    ModeAttribute{
      mode:  ConfigMode::Config,
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

    Ok(ModeAttribute{ mode,label, path})
  }
  //----------------------------------------------------------------------------
  pub fn to_string(&self) -> String{
  
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

    format!("#[{}{}{}]",mode_str,label_str,path_str)

  }
}

#[derive(PartialEq, Debug, Clone)]
pub enum ConfigMode{
  Clusters,
  Config,
  Filter,
  Spins,
  Structures,
  Tensors,
}

impl ConfigMode{
  pub fn from(token: Token) -> Result<Self,CluEError>{

    match token{
      Token::Clusters => Ok(ConfigMode::Clusters),
      Token::Config => Ok(ConfigMode::Config),
      Token::Filter => Ok(ConfigMode::Filter),
      Token::Spins => Ok(ConfigMode::Spins),
      Token::Structures => Ok(ConfigMode::Structures),
      Token::Tensors => Ok(ConfigMode::Tensors),
      _ => Err(CluEError::ConfigModeNotRecognized(token.to_string())),
    }
  }
  //----------------------------------------------------------------------------
  pub fn to_string(&self) -> String{
    match self{
      ConfigMode::Clusters => "clusters".to_string(),
      ConfigMode::Config => "config".to_string(),
      ConfigMode::Filter => "filter".to_string(),
      ConfigMode::Spins => "spins".to_string(),
      ConfigMode::Structures => "structures".to_string(),
      ConfigMode::Tensors => "tensors".to_string(),
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
    assert_eq!(identify_token("!").unwrap(), Token::Bang);
    assert_eq!(identify_token("*/").unwrap(), Token::BlockCommentEnd);
    assert_eq!(identify_token("/*").unwrap(), Token::BlockCommentStart);
    assert_eq!(identify_token("clusters").unwrap(), Token::Clusters);
    assert_eq!(identify_token(",").unwrap(), Token::Comma);
    assert_eq!(identify_token("config").unwrap(), Token::Config);
    assert_eq!(identify_token("}").unwrap(), Token::CurlyBracketClose);
    assert_eq!(identify_token("{").unwrap(), Token::CurlyBracketOpen);
    assert_eq!(identify_token("element").unwrap(), Token::Element);
    assert_eq!(identify_token("\n").unwrap(), Token::EOL);
    assert_eq!(identify_token("=").unwrap(), Token::Equals);
    assert_eq!(identify_token("filter").unwrap(), Token::Filter);
    assert_eq!(identify_token(">").unwrap(), Token::GreaterThan);
    assert_eq!(identify_token(">=").unwrap(), Token::GreaterThanEqualTo);
    assert_eq!(identify_token("^").unwrap(), Token::Hat);
    assert_eq!(identify_token("in").unwrap(), Token::In);
    assert_eq!(identify_token("label").unwrap(), Token::Label);
    assert_eq!(identify_token("<").unwrap(), Token::LessThan);
    assert_eq!(identify_token("<=").unwrap(), Token::LessThanEqualTo);
    assert_eq!(identify_token("//").unwrap(), Token::LineComment);
    assert_eq!(identify_token("magnetic_field").unwrap(), Token::MagneticField);
    assert_eq!(identify_token("-").unwrap(), Token::Minus);
    assert_eq!(identify_token("not").unwrap(), Token::Not);
    assert_eq!(identify_token("!=").unwrap(), Token::NotEqual);
    assert_eq!(identify_token("not in").unwrap(), Token::NotIn);
    assert_eq!(identify_token("path").unwrap(), Token::Path);
    assert_eq!(identify_token(")").unwrap(), Token::ParenthesisClose);
    assert_eq!(identify_token("(").unwrap(), Token::ParenthesisOpen);
    assert_eq!(identify_token("+").unwrap(), Token::Plus);
    assert_eq!(identify_token("radius").unwrap(), Token::Radius);
    assert_eq!(identify_token("residue").unwrap(), Token::Residue);
    assert_eq!(identify_token("rng_seed").unwrap(), Token::RNGSeed);
    assert_eq!(identify_token(";").unwrap(), Token::Semicolon);
    assert_eq!(identify_token("#").unwrap(), Token::Sharp);
    assert_eq!(identify_token("/").unwrap(), Token::Slash);
    assert_eq!(identify_token("spins").unwrap(), Token::Spins);
    assert_eq!(identify_token("]").unwrap(), Token::SquareBracketClose);
    assert_eq!(identify_token("[").unwrap(), Token::SquareBracketOpen);
    assert_eq!(identify_token("structures").unwrap(), Token::Structures);
    assert_eq!(identify_token("tensors").unwrap(), Token::Tensors);
    assert_eq!(identify_token("*").unwrap(), Token::Times);
    assert_eq!(identify_token("tunnel_splitting").unwrap(), 
        Token::TunnelSplitting);
    assert_eq!(identify_token("type").unwrap(), Token::Type);
    assert_eq!(identify_token(" ").unwrap(), Token::Whitespace);
    assert_eq!(identify_token("write_structure_pdb").unwrap(), 
        Token::WriteStructurePDB);
    assert_eq!(identify_token("indices").unwrap(), Token::Indices); 
    assert_eq!(identify_token("elements").unwrap(), Token::Elements); 
    assert_eq!(identify_token("serials").unwrap(), Token::Serials); 
    assert_eq!(identify_token("residues").unwrap(), Token::Residues);
    assert_eq!(identify_token("residue_sequence_numbers").unwrap(), 
        Token::ResSeqNums);
    assert_eq!(identify_token("bonded_indices").unwrap(), Token::BondedIndices);
    assert_eq!(identify_token("bonded_elements").unwrap(),
        Token::BondedElements);

  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


