
use crate::clue_errors::CluEError;
use crate::config::token::Token;

#[derive(Debug,Clone,PartialEq)]
pub enum UnitOfClustering{
  Set,
  Spin,
}

impl UnitOfClustering{
  //----------------------------------------------------------------------------
  pub fn from(unit_of_clustering: &[Token], line_number: usize) 
      -> Result<Self,CluEError>
  {
    match unit_of_clustering[0]{
      Token::Set => Ok(UnitOfClustering::Set),
      Token::Spin => Ok(UnitOfClustering::Spin),
      _ => Err(CluEError::CannotParseUnitOfClustering(line_number,
            unit_of_clustering[0].to_string())),
    }
  }
  //----------------------------------------------------------------------------
}
