use crate::config::token::*;
use crate::CluEError;
use crate::config::ModeAttribute;
use crate::config::token_stream;
use crate::config::token_stream::split_on_token;
use crate::config::token_algebra::*;
use crate::config::to_i32_token;
use crate::config::token_stream::{find_brackets,find_outermost_parentheses};
use crate::physical_constants::Element;
use crate::space_3d::Vector3D;
use crate::structure::particle_filter::{SecondaryParticleFilter,
  VectorSpecifier};

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
      let lhs = vec![Token::Mode(mode)];
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
    let mut rhs = Vec::<Token>::with_capacity(tokens.len() - idx );

    for (ii,token) in tokens.iter().enumerate(){

      match ii.cmp(&idx){
        std::cmp::Ordering::Less => lhs.push(token.clone()),
        std::cmp::Ordering::Equal => (),
        std::cmp::Ordering::Greater => rhs.push(token.clone()),
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

  let mut index_opt: Option<usize> = None;

  for (ii, token) in tokens.iter().enumerate(){
    if is_relational_operator(token){
      if index_opt.is_none(){
        index_opt = Some(ii); 
      }else{
        return Err(CluEError::TooManyRelationalOperators(line_number));
      }
    }
  }

  Ok(index_opt)
}

//------------------------------------------------------------------------------
pub fn check_target_option<T>(target: &Option<T>, expression: &TokenExpression)
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
pub fn check_target_vector<T>(target: &Vec<T>, expression: &TokenExpression)
  -> Result<(),CluEError>
{
  let already_set = ||{
    CluEError::OptionAlreadySet(
        expression.line_number, expression.lhs[0].to_string()) };

  if !target.is_empty(){
    return Err(already_set());
  }

  Ok(())
}

//------------------------------------------------------------------------------
pub fn set_to_some_bool(target: &mut Option<bool>, expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target_option(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  if tokens.is_empty(){
    return Err(CluEError::NoRHS(expression.line_number));
  }else if tokens.len() > 1{
    return Err(CluEError::TooManyRHSArguments(expression.line_number));
  }

  match tokens[0]{
    Token::False => *target = Some(false),
    Token::True => *target = Some(true),
    _ => return Err(CluEError::ExpectedBoolRHS(expression.line_number)),
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_some_f64(target: &mut Option<f64>, expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target_option(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_f64_token(tokens, expression.line_number)?;
  if let Token::Float(x) = value_token {
    *target = Some(x);
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_vec_f64(target: &mut Vec<f64>, 
    expression: &TokenExpression)-> Result<(),CluEError>
{

  check_target_vector(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;
  let value_token = to_f64_token(tokens, expression.line_number)?;
  match value_token{
    Token::VectorF64(vec) => *target = vec,
    Token::Float(x) => *target = vec![x],
    _ => return Err(CluEError::NoRHS(expression.line_number)),  
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_some_i32(target: &mut Option<i32>, expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target_option(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_i32_token(tokens, expression.line_number)?;
  if let Token::Int(a) = value_token {
    *target = Some(a);
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_some_usize(target: &mut Option<usize>, 
    expression: &TokenExpression)-> Result<(),CluEError>
{

  check_target_option(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  let value_token = to_i32_token(tokens, expression.line_number)?;
  if let Token::Int(a) = value_token {
    if a < 0{
      return Err(CluEError::ExpectedNonNegativeIntRHS(expression.line_number));
    }
    *target = Some(a as usize);
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_vec_usize(target: &mut Vec<usize>, 
    expression: &TokenExpression)-> Result<(),CluEError>
{

  check_target_vector(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;
  let value_token = to_i32_token(tokens, expression.line_number)?;
  match value_token{
    Token::VectorI32(vec) => {
      *target = Vec::<usize>::with_capacity(vec.len());
      for &a in vec.iter(){
        if a < 0{
          return Err(
              CluEError::ExpectedNonNegativeIntRHS(expression.line_number));
        }
        target.push(a as usize);
      }
    },
    Token::Int(a) => {
      if a < 0{
        return Err(
            CluEError::ExpectedNonNegativeIntRHS(expression.line_number));
      }
      *target = vec![a as usize];
    },
    _ => return Err(CluEError::NoRHS(expression.line_number)),  
  }

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_some_string(target: &mut Option<String>, 
    expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target_option(target,expression)?;

  let Some(rhs) = &expression.rhs else {
    return Err(CluEError::NoRHS(expression.line_number));
  };

  if rhs.is_empty(){
    return Err(CluEError::NoRHS(expression.line_number));
  }else if rhs.len()>1{
    return Err(CluEError::TooManyRHSArguments(expression.line_number));
  }

  *target = Some(rhs[0].to_string());
  /*
  if let Token::UserInputValue(value) = &rhs[0]{
    *target = Some(value.clone());
  }else{
    return Err(CluEError::CannotParseRHS(expression.line_number));
  }
  */

  Ok(())
}
//------------------------------------------------------------------------------
pub fn set_to_some_vector3d(
    target: &mut Option<Vector3D>, expression: &TokenExpression)
  -> Result<(),CluEError>
{

  check_target_option(target,expression)?;

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
// TODO: add more descriptive errors.
pub fn set_to_some_vector_specifier(
    target: &mut Option<VectorSpecifier>, expression: &TokenExpression,
    label: &str)
  -> Result<(),CluEError>
{
  
  check_target_option(target,expression)?;

  let tokens = token_stream::extract_rhs(expression)?;

  if tokens.is_empty(){
    return Err(CluEError::NoRHS(expression.line_number));
  }

  let vector_specifier: VectorSpecifier;

  match tokens[0]{
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Token::Diff => {
      let Some( (idx0,idx1) ) = find_outermost_parentheses(&tokens, 
          expression.line_number)? else{  
         return Err(CluEError::CannotConvertToVector(
              expression.line_number));
      };

      let args = split_on_token(tokens[idx0+1..idx1].to_vec(), Token::Comma);

      if args.len() != 2{
        return Err(CluEError::CannotConvertToVector(
            expression.line_number));
      }

      let mut secondary_particle_filters 
        = Vec::<SecondaryParticleFilter>::with_capacity(2);

      let mut labels = Vec::<String>::with_capacity(2);

      for arg in args.iter(){
        if arg.is_empty(){
          return Err(CluEError::CannotConvertToVector(
            expression.line_number));
        }

        let sec_fltr = SecondaryParticleFilter::from(&arg[0].to_string())?;

        match sec_fltr{
          SecondaryParticleFilter::Particle => labels.push(label.to_string()),
          _ => {
            if arg.len() != 4{
              return Err(CluEError::MissingFilterArgument(
                    expression.line_number, sec_fltr.to_string() ));
            }
            labels.push(arg[2].to_string());
          },
        }
        secondary_particle_filters.push(sec_fltr);
      }

      vector_specifier = VectorSpecifier::Diff(
         secondary_particle_filters[0].clone(),labels[0].clone(),
         secondary_particle_filters[1].clone(),labels[1].clone());
    },
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Vector => {
      let Some( (idx0,idx1) ) = find_brackets(&tokens, expression.line_number)?
      else{  
        return Err(CluEError::CannotConvertToVector(
              expression.line_number));
      };

      let value_token = to_f64_token((tokens[idx0..=idx1]).to_vec(), 
          expression.line_number)?;
      
      if let Token::VectorF64(vec) = value_token{ 
        if vec.len() != 3 {
          return Err(CluEError::WrongVectorLength(
                expression.line_number,3,vec.len()));
        }

        vector_specifier = VectorSpecifier::Vector(
            Vector3D::from([vec[0],vec[1],vec[2]])
            );

      } else{
        return Err(CluEError::CannotConvertToVector(
              expression.line_number));
      }
    },
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    _ => {
      let mut vec3d_opt: Option<Vector3D> =  None;
      set_to_some_vector3d(&mut vec3d_opt, expression)?;
      if let Some(r) = vec3d_opt{
        vector_specifier = VectorSpecifier::Vector(r);
      }else{
        return 
          Err(CluEError::UnrecognizedVectorSpecifier(tokens[0].to_string() ) );
      }
    },
  }
  
  *target = Some(vector_specifier);

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
pub fn get_function_arguments(tokens: &[Token],line_number: usize) 
  -> Result<Vec::<Token>,CluEError>
{
  let parentheses = find_outermost_parentheses(tokens,line_number)?;
  let Some((arg_start,arg_end)) = parentheses else{
    return Err(CluEError::NoArgument(line_number));
  };

  let arg_length = arg_end - arg_start - 1;
  let mut arg_tokens = Vec::<Token>::with_capacity(arg_length);
  for token in tokens.iter().take(arg_end).skip(arg_start+1){
    arg_tokens.push(token.clone());
  }

  Ok(arg_tokens)
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

    let tokens = vec![Token::TunnelSplitting, Token::SquareBracketOpen, 
        Token::UserInputValue("1H".to_string()), Token::SquareBracketClose,
        Token::Equals, Token::Float(80.0e3)];
    let expression = TokenExpression::from(tokens.clone(),0).unwrap();

    for ii in 0..4{
      assert_eq!(expression.lhs[ii],tokens[ii]);
    }
    if let Some(rhs) = expression.rhs{
      assert_eq!(rhs[0],tokens[5]);
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
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_to_some_vector_specifier(){
    
    let mut target: Option<VectorSpecifier> = None;

    let expression = TokenExpression{ 
      lhs: vec![Token::HyperfineX], 
      rhs: Some(vec![Token::Vector, 
          Token::ParenthesisOpen, Token::SquareBracketOpen, Token::Minus, 
          Token::UserInputValue("1".to_string()), Token::Comma, 
          Token::UserInputValue("0".to_string()), Token::Comma, 
          Token::UserInputValue("1".to_string()), Token::SquareBracketClose, 
          Token::ParenthesisClose]), 
      relationship: Some(Token::Equals), 
      line_number: 3 };

    let label = "label";

    set_to_some_vector_specifier(&mut target, &expression,&label).unwrap();

    assert_eq!(target,Some(
          VectorSpecifier::Vector(
            Vector3D::from([-1.0, 0.0, 1.0])
            )
          ));
  }
}
