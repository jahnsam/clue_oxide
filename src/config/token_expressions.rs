use crate::config::token::*;
use crate::CluEError;
use crate::config::ModeAttribute;
use crate::config::token_stream;
use crate::config::token_algebra::*;
use crate::space_3d::Vector3D;
use crate::config::to_i32_token;
use crate::physical_constants::Element;

// TokenExpression holds a line of input, organized for easy parsing.
#[derive(Debug,Clone,Default,PartialEq)]
pub struct TokenExpression{
  pub lhs: Vec::<Token>,
  pub rhs: Option<Vec::<Token>>,
  pub relationship: Option<Token>,
  pub line_number: usize,
}

impl TokenExpression{
  //----------------------------------------------------------------------------
  // This function processes a vector of tokens into a TokenExpression.
  pub fn from(tokens: Vec::<Token>, line_number:usize) -> Result<Self,CluEError>{

    if tokens[0]==Token::Sharp{
      let mode = ModeAttribute::from(tokens)?;
      let mut lhs = Vec::<Token>::with_capacity(1);
      lhs.push(Token::Mode(mode));
      return Ok(TokenExpression{lhs, rhs: None, relationship: None,line_number});
    }

    let idx: usize; 
    let idx_option = find_lhs_rhs_delimiter_index(&tokens, line_number)?;
    
    if let Some(index) = idx_option{
      idx = index;
    }else{
      idx = tokens.len();
    }

    let mut lhs = Vec::<Token>::with_capacity(idx);
    let mut rhs = Vec::<Token>::with_capacity(tokens.len() - idx - 1);

    for (ii,token) in tokens.iter().enumerate(){

      if ii < idx{
        lhs.push(token.clone());
      }else if ii > idx{
        rhs.push(token.clone());
      }
    }

    if rhs.is_empty(){
      return Ok(TokenExpression{lhs, rhs: None,relationship: None,line_number})
    }

    let rhs = Some(rhs);
    let relationship = Some(tokens[idx].clone());
    Ok(TokenExpression{lhs, rhs, relationship, line_number})

  }
}
//------------------------------------------------------------------------------
pub fn find_lhs_rhs_delimiter_index(tokens: &[Token], line_number: usize) 
  -> Result<Option<usize>,CluEError>{
  let mut found_token = false;
  let mut index: usize = 0;

  for ii in 0..tokens.len(){
    if is_relational_operator(&tokens[ii]){
      if !found_token{
        found_token = true;
        index = ii;
      }else{
        return Err(CluEError::TooManyRelationalOperators(line_number));
      }
    }
  }
  if found_token{
    return Ok(Some(index));
  }else{
    //return Err(CluEError::NoRelationalOperators(line_number));
    return Ok(None);
  }
}

//------------------------------------------------------------------------------
pub fn check_target<T>(target: &Option<T>, expression: &TokenExpression)
  -> Result<(),CluEError>
{
  let already_set = ||{
    CluEError::OptionAlreadySet(
        expression.line_number, expression.lhs[0].to_string()) };

  if let Some(_value) = target{
    return Err(already_set());
  }

  Ok(())
}

//------------------------------------------------------------------------------
pub fn set_to_some_f64(target: &mut Option<f64>, expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_f64_token(tokens, expression.line_number)?;
  if let Token::VectorF64(vec) = value_token {
    if vec.len() != 1{
      return Err(CluEError::ExpectedFloatRHS(expression.line_number));
    }
    *target = Some(vec[0]);
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_some_i32(target: &mut Option<i32>, expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_i32_token(tokens, expression.line_number)?;
  if let Token::VectorI32(vec) = value_token {
    if vec.len() != 1{
      return Err(CluEError::ExpectedIntRHS(expression.line_number));
    }
    *target = Some(vec[0]);
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_some_vector3d(
    target: &mut Option<Vector3D>, expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_f64_token(tokens, expression.line_number)?;
  if let Token::VectorF64(vec) = value_token {
    if vec.len() != 3{
      return Err(CluEError::ExpectedVecOfNFloatsRHS(expression.line_number,3));
    }
    *target = Some(Vector3D::from([vec[0],vec[1],vec[2]]));
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn vec_tokens_to_vec_strings(tokens: Vec::<Token>)
  -> Result<Vec::<String>,CluEError>
{
  /*
  let str_token: Vec::<String> = tokens.into_iter()
    .map(|tok| tok.to_string())
    .collect();
  */
  let mut value_token = Vec::<String>::with_capacity(tokens.len());
  for tok in tokens{
    if let Token::UserInputValue(tok_str) = tok{
      value_token.push(tok_str);
    }
  }

  Ok(value_token)
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
pub fn vec_tokens_to_vec_elements(tokens: Vec::<Token>)
  -> Result<Vec::<Element>,CluEError>
{
  let mut value_token = Vec::<Element>::with_capacity(tokens.len());
  for tok in tokens.iter(){
    if let Token::UserInputValue(el_str) = tok{
      let el = Element::from(el_str)?;
      value_token.push(el);
    }
  }

  Ok(value_token)
}
//------------------------------------------------------------------------------


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;  
    //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_TokenExpression(){
    let tokens = vec![Token::Float(1.0), Token::Plus, Token::Float(2.0),
    Token::Equals, Token::Float(3.0)];

    let expression = TokenExpression::from(tokens.clone(),0).unwrap();

    assert_eq!(expression.relationship,Some(tokens[3].clone()));
    for ii in 0..3{
      assert_eq!(expression.lhs[ii],tokens[ii]);
    }
    if let Some(rhs) = expression.rhs{
      assert_eq!(rhs[0],tokens[4]);
    }
    else{
      panic!("Expected Some(rhs), but found None.");
    }
  }

  //----------------------------------------------------------------------------

  #[test]
  fn test_find_lhs_rhs_delimiter_index(){

    let mut tokens = vec![Token::Float(1.0), Token::Plus, Token::Float(2.0),
    Token::Comma, Token::Float(3.0)];

    let result = find_lhs_rhs_delimiter_index(&tokens,0).unwrap();
    assert_eq!(result,None);

    tokens[3] = Token::Equals;
    let result = find_lhs_rhs_delimiter_index(&tokens,0);
    assert_eq!(result,Ok(Some(3)));

    tokens[0] = Token::In;
    let result = find_lhs_rhs_delimiter_index(&tokens,0);
    assert_eq!(result,Err(CluEError::TooManyRelationalOperators(0)));
  }
}
