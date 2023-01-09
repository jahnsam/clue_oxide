use crate::clue_errors::*;
use crate::config::token::*;
use crate::config::token_stream::*;
//use std::ops::{Add,Sub,Mul,Div};

//------------------------------------------------------------------------------
// The function answers the question 
// "does this token represent a number/set of numbers?".
fn is_numeric(token: &Token) -> bool{
  match token{
    Token::Float(_x) => true,
    Token::Int(_n) => true,
    Token::VectorF64(_v) => true,
    Token::VectorI32(_v) => true,
    _ => false

  }
}
//------------------------------------------------------------------------------
// This function takes three tokens, two numeric tokens separated by a basic
// operator token and returns the numeric token that results from the 
// operation.
pub fn basic_token_algebraic_operations(tokens: [Token;3],line_number: usize) 
 ->Result<Token,CluEError>{

  if tokens.len() != 3{
    return Err(CluEError::WrongVectorLength(line_number, 3, tokens.len()));
  }

  match tokens[1] {
    Token::Plus  => tokens[0].clone()  + tokens[2].clone(),
    Token::Minus  => tokens[0].clone() - tokens[2].clone(),
    Token::Times  => tokens[0].clone() * tokens[2].clone(),
    Token::Slash  => tokens[0].clone() / tokens[2].clone(),
    Token::Hat  => tokens[0].clone().pow(tokens[2].clone()),
    _ => Err(CluEError::NotAnOperator(line_number, tokens[1].to_string())),
  }

}
//------------------------------------------------------------------------------
// This function takes a vector of tokens and two tokens with the same
// order of operations, such as + and - or * and /.
// The return vector of tokens is the same as the input, except the operations
// are evaluated.
pub fn contract_operator_inverse_operator_tokens(
    tokens: Vec::<Token>,line_number: usize,
    op: Token, inv_op: Token)
 -> Result<Vec::<Token>, CluEError>{

  let mut out = Vec::<Token>::with_capacity(tokens.len());

  if tokens.len() < 3{
    return Ok(tokens);
  }

  out.push(tokens[0].clone());

  let mut skip_next = false;

  for ii in 1..tokens.len()-1 {

    if skip_next{
      skip_next = false;
      continue;
    }

    if (tokens[ii] == op || tokens[ii] == inv_op) 
      && is_numeric(&out[out.len()-1]) && is_numeric(&tokens[ii+1]){

      let new_token = basic_token_algebraic_operations(
       [out[out.len()-1].clone(), tokens[ii].clone(), tokens[ii+1].clone()],
       line_number)?;

      let idx = out.len() - 1;
      out[idx] = new_token;
      skip_next = true;

    }else{
      out.push(tokens[ii].clone());
    }

  }
  if !skip_next{
    out.push(tokens[tokens.len()-1].clone());
  }
  
  Ok(out)
}
//------------------------------------------------------------------------------
// This function return the input vector of tokens, but with all 
// exponentiations evaluated.
pub fn contract_exponentiation_tokens(tokens: Vec::<Token>,line_number: usize)
 -> Result<Vec::<Token>, CluEError>{

   contract_operator_inverse_operator_tokens(tokens, line_number,
       Token::Hat, Token::Hat)
 }
//------------------------------------------------------------------------------
// This function return the input vector of tokens, but with all 
// multiplications and divisions evaluated.
pub fn contract_multiply_divide_tokens(tokens: Vec::<Token>,line_number: usize)
 -> Result<Vec::<Token>, CluEError>{

   contract_operator_inverse_operator_tokens(tokens, line_number,
       Token::Times, Token::Slash)
 }
//------------------------------------------------------------------------------
// This function return the input vector of tokens, but with all 
// additions and subtractions evaluated.
pub fn contract_add_subtract_tokens(tokens: Vec::<Token>,line_number: usize)
 -> Result<Vec::<Token>, CluEError>{

   contract_operator_inverse_operator_tokens(tokens, line_number,
       Token::Plus, Token::Minus)
 }
 
//------------------------------------------------------------------------------
// This function takes a vector of tokens representing and algebreic sequence
// that does not contain any parentheses, and returns the evaluations.
fn contract_emdas(tokens: Vec::<Token>, line_number: usize) 
 -> Result<Token, CluEError>
{

  // gate keeping section
  if tokens.is_empty(){
    return Err(CluEError::IndexOutOfBounds(line_number,0, 0));
  }

  let n_pairs = count_delimiter_pairs(&tokens,
      &Token::ParenthesisOpen, &Token::ParenthesisClose,
      line_number)?;

  if n_pairs != 0 {
    return Err(CluEError::InvalidToken(line_number,"parenthaesis".to_string()));
  }



  // Evaluate exponentials.
  let tokens = contract_exponentiation_tokens(tokens,line_number)?;

  if count_token(&Token::Hat,&tokens) != 0{
    return Err(CluEError::InvalidToken(line_number,Token::Hat.to_string()));
  }


  // Evaluate mutiplications and divisions.
  let mut tokens = contract_multiply_divide_tokens(tokens,line_number)?;

  if count_token(&Token::Times,&tokens) != 0{
    return Err(CluEError::InvalidToken(line_number,Token::Times.to_string()));
  }
  if count_token(&Token::Slash,&tokens) != 0{
    return Err(CluEError::InvalidToken(line_number,Token::Slash.to_string()));
  }


  // Evaluate additions and subtractions.
  if tokens[0] == Token::Plus || tokens[0] == Token::Minus{
    tokens = add_token_at_index(Token::Float(0.0),0,tokens, line_number)?;
  }

  let tokens = contract_add_subtract_tokens(tokens,line_number)?;
  

  // Check for errors.
  if tokens.len() != 1 {
    return Err(CluEError::CannotCombineTokens(line_number));
  }

  // Get return value
  Ok(tokens[0].clone())
}  

//------------------------------------------------------------------------------
// This function takes a vector of tokens representing and algebreic sequence
// and returns the evaluations.
fn contract_pemdas(mut tokens: Vec::<Token>, line_number: usize) 
 -> Result<Token, CluEError>
{

  let mut is_fully_contracted = false;

  while !is_fully_contracted {
    let idx_option = find_deepest_parentheses(&tokens, line_number)?;

    let index0: usize;
    let index1: usize;

    let found_parentheses: bool;
    if let Some((idx0,idx1)) = idx_option{

      index0 = idx0+1;
      index1 = idx1-1;
      found_parentheses = true;

    }else{ // no parentheses case
       index0 = 0;
       index1 = tokens.len() - 1;
       found_parentheses = false;
       is_fully_contracted = true;
    }

    let mut new_tokens = Vec::<Token>::with_capacity( tokens.len() );

    // Include everthing up to the opening parenthesis.
    if found_parentheses {
      new_tokens.append(&mut tokens[0..index0-1].to_vec());
    }

    // Contract the parentheses.
    if index1 >= index0 {
      let toks = contract_emdas(tokens[index0..=index1].to_vec(),line_number)?;
      new_tokens.push(toks);
    }

    // Add everthing after the closing parenthesis.
    if found_parentheses && index1+2 < tokens.len()-1{
      new_tokens.append(
          &mut tokens[index1+2..tokens.len()].to_vec()
          );
    }

    tokens = new_tokens;
  }

  // Check for errors.
  if tokens.len() != 1 {
    return Err(CluEError::CannotCombineTokens(line_number));
  }

  // Get return value
  Ok(tokens[0].clone())
}

//------------------------------------------------------------------------------
// This function takes a sequence of tokens and combines them to a single
// token representing a vector of floats.
fn to_float_vector(tokens: Vec::<Token>, line_number: usize) 
 -> Result<Token, CluEError>
{

  let vector_elements = get_vector_elements(tokens,line_number)?;

  // Initialize the output.
  let mut out = Vec::<f64>::with_capacity(vector_elements.len());

  // Assign output. 
  for token_element in vector_elements{

    // Contract the parentheses.
    let token = contract_emdas(token_element,line_number)?;
    
    match token{
      Token::Float(x) => out.push(x),
      Token::Int(a) => out.push(a as f64),
      _ => return Err(CluEError::CannotConvertToFloat(line_number, 
            token.to_string() ) ),
    }
  }
  
  Ok(Token::VectorF64(out))

}
//------------------------------------------------------------------------------
// This function takes a sequence of tokens and combines them to a single
// token representing a vector of ints.
fn to_integer_vector(tokens: Vec::<Token>, line_number: usize) 
 -> Result<Token, CluEError>
{

  let vector_elements = get_vector_elements(tokens,line_number)?;

  // Initialize the output.
  let mut out = Vec::<i32>::with_capacity(vector_elements.len());

  // Assign output. 
  for token_element in vector_elements{

    // Contract the parentheses.
    let token = contract_emdas(token_element,line_number)?;
    
    match token{
      Token::Int(a) => out.push(a),
      _ => return Err(CluEError::CannotConvertToFloat(line_number, 
            token.to_string() ) ),
    }
  }
  
  Ok(Token::VectorI32(out))

}
//------------------------------------------------------------------------------
// This function takes a sequence of tokens and combines them to a single
// token representing a vector of either floats or ints dependending on the
// input bool.
fn contract_numeric_vectors(tokens: Vec::<Token>, as_f64: bool,
    line_number: usize)
-> Result<Token, CluEError>
{
  // Find "[" and "]".
  let indices_open = find_token(&Token::SquareBracketOpen, &tokens);
  let indices_close = find_token(&Token::SquareBracketClose, &tokens);

  // Ensure there every "[" is paired with a "]".
  if indices_open.len() != indices_close.len(){
    return Err(CluEError::UnmatchedDelimiter(line_number));
  }

  // Exclude nested bracket such as "[ [] ]".
  for (ii, idx_o) in indices_open.iter().enumerate(){
    let idx_c = indices_close[ii];
    
    if *idx_o >= idx_c 
      || (ii+1<indices_open.len() 
          && indices_open[ii+1] <= idx_c){
      return Err(CluEError::UnmatchedDelimiter(line_number));
      } 
  }

  // Initialized output.
  let mut out = Vec::<Token>::with_capacity(tokens.len());
  let mut read_from = 0;
  let mut read_to;

  // Contract vectors.
  for (ii, idx_o) in indices_open.iter().enumerate(){
    let idx_c = indices_close[ii];

    read_to = *idx_o;
    for idx in read_from .. read_to{
      out.push(tokens[idx].clone());
    }
    read_from = idx_c + 1;


    let array: Token;
    if as_f64{
      array = to_float_vector(tokens[*idx_o..=idx_c].to_vec(),line_number)?;
    }else{
      array = to_integer_vector(tokens[*idx_o..=idx_c].to_vec(),line_number)?;
    }
    out.push(array);
  }

  read_to = tokens.len();
  for idx in read_from .. read_to{
    out.push(tokens[idx].clone());
  }

  contract_pemdas(out, line_number)
}

//------------------------------------------------------------------------------
// This function reads the input vector of tokens and tries to interpret all
// unidentified floats and returns a Token::VectorF64.
pub fn to_f64_token(tokens: Vec::<Token>, line_number: usize)
-> Result<Token, CluEError>
{
  let tokens = read_strings_as_floats(tokens, line_number)?; 
  contract_numeric_vectors(tokens, true, line_number)
}
//------------------------------------------------------------------------------
// This function reads the input vector of tokens and tries to interpret all
// unidentified floats and returns a Token::VectorI32.
pub fn to_i32_token(tokens: Vec::<Token>, line_number: usize)
-> Result<Token, CluEError>
{
  let tokens = read_strings_as_integers(tokens, line_number)?; 
  contract_numeric_vectors(tokens, false, line_number)
}
//------------------------------------------------------------------------------


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
 
  //----------------------------------------------------------------------------
  #[test]
  fn test_basic_token_algebraic_operations(){
  
    let mut tokens = vec![Token::Float(2.0), Token::Plus, Token::Float(3.0),
     Token::Equals, Token::Float(3.0)];
    
    let answer = basic_token_algebraic_operations(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 + 3.0));

    tokens[1] = Token::Minus;
    let answer = basic_token_algebraic_operations(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 - 3.0));

    tokens[1] = Token::Times;
    let answer = basic_token_algebraic_operations(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 * 3.0));

    tokens[1] = Token::Slash;
    let answer = basic_token_algebraic_operations(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 / 3.0));

    tokens[1] = Token::Hat;
    let answer = basic_token_algebraic_operations(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(8.0));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_contract_multiply_divide_tokens(){

    let tokens = vec![Token::Float(2.0), Token::Times, Token::Float(3.0),
     Token::Comma, Token::Float(3.0)];

    let tokens = contract_multiply_divide_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 3);
    assert_eq!(tokens, 
        vec![Token::Float(6.0), Token::Comma, Token::Float(3.0)]);

    let tokens = vec![Token::Float(2.0), Token::Slash, Token::Float(3.0),
     Token::Comma, Token::Float(3.0)];

    let tokens = contract_multiply_divide_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 3);
    assert_eq!(tokens, 
        vec![Token::Float(2.0/3.0), Token::Comma, Token::Float(3.0)]);

    let tokens = vec![Token::Float(6.0), Token::Times, Token::Float(2.0), 
        Token::Slash, Token::Float(3.0), Token::Times, Token::Float(5.0)];

    let tokens = contract_multiply_divide_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 1);
    assert_eq!(tokens, vec![Token::Float(6.0*2.0/3.0*5.0)]);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_contract_add_subtract_tokens(){

    let tokens = vec![Token::Float(2.0), Token::Plus, Token::Float(3.0),
     Token::Comma, Token::Float(3.0)];

    let tokens = contract_add_subtract_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 3);
    assert_eq!(tokens, 
        vec![Token::Float(5.0), Token::Comma, Token::Float(3.0)]);

    let tokens = vec![Token::Float(2.0), Token::Minus, Token::Float(3.0),
     Token::Comma, Token::Float(3.0)];

    let tokens = contract_add_subtract_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 3);
    assert_eq!(tokens, 
        vec![Token::Float(2.0 - 3.0), Token::Comma, Token::Float(3.0)]);

    let tokens = vec![Token::Float(6.0), Token::Plus, Token::Float(2.0), 
        Token::Minus, Token::Float(3.0), Token::Plus, Token::Float(5.0)];

    let tokens = contract_add_subtract_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 1);
    assert_eq!(tokens, vec![Token::Float(6.0 + 2.0 - 3.0 + 5.0)]);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_contract_exponentiation_tokens(){

    let tokens = vec![Token::Float(2.0), Token::Hat, Token::Float(3.0),
     Token::Comma, Token::Float(3.0)];

    let tokens = contract_exponentiation_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 3);

    assert_eq!(tokens, 
        vec![Token::Float(8.0), Token::Comma, Token::Float(3.0)]);

    let tokens = vec![Token::Float(-2.0), Token::Hat, Token::Float(3.0),
     Token::Comma, Token::Float(3.0)];

    let tokens = contract_exponentiation_tokens(tokens,0).unwrap();
    assert_eq!(tokens.len(), 3);

    assert_eq!(tokens, 
        vec![Token::Float(-8.0), Token::Comma, Token::Float(3.0)]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_contract_emdas(){

    let tokens = vec![ Token::Float(1.0) ];
    let result = contract_emdas(tokens,0).unwrap();
    assert_eq!(result, Token::Float(1.0));


    let tokens = vec![
        Token::Float(2.0), Token::Hat, Token::Float(5.0),
        Token::Plus, Token::Float(3.0), Token::Times, Token::Float(12.0),
        Token::Slash, Token::Float(9.0), Token::Minus, Token::Float(-4.0)];
    
    let result = contract_emdas(tokens,0).unwrap();

    assert_eq!(result, Token::Float(40.0));

    let tokens = vec![Token::Minus,
        Token::Float(2.0), Token::Hat, Token::Float(5.0),
        Token::Plus, Token::Float(3.0), Token::Times, Token::Float(12.0),
        Token::Slash, Token::Float(9.0), Token::Minus, Token::Float(-4.0)];
    
    let result = contract_emdas(tokens,0).unwrap();

    assert_eq!(result, Token::Float(-24.0));
  
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_contract_pemdas(){

    let tokens = vec![ Token::ParenthesisOpen,
        Token::Float(2.0), Token::Hat, Token::Float(5.0),
        Token::Plus, Token::Float(3.0), Token::Times, Token::Float(12.0),
        Token::Slash, Token::Float(9.0), Token::Minus, Token::Float(-4.0),
        Token::ParenthesisClose];
    
    let result = contract_pemdas(tokens,0).unwrap();

    assert_eq!(result, Token::Float(40.0));

    let tokens = vec![
      Token::ParenthesisOpen,
        Token::Minus,Token::Float(2.0), 
      Token::ParenthesisClose, Token::Hat, Token::Float(2.0),
      Token::Plus, Token::Float(3.0), 
      Token::Times, Token::Float(12.0),
      Token::Slash, Token::Float(9.0), Token::Minus, Token::Float(-4.0)];

    let result = contract_pemdas(tokens,0).unwrap();
    assert_eq!(result, Token::Float(12.0));


    let tokens = vec![
      Token::ParenthesisOpen,
        Token::ParenthesisOpen,
          Token::Minus,Token::Float(2.0), 
        Token::ParenthesisClose, Token::Hat, Token::Float(2.0),
        Token::Plus, Token::Float(3.0), 
      Token::ParenthesisClose,
      Token::Times, Token::Float(12.0),
      Token::Slash, Token::Float(7.0), Token::Minus, 
      Token::ParenthesisOpen,
        Token::Float(-2.0), Token::Plus, Token::Float(-2.0),
      Token::ParenthesisClose,
    ];
    
    let result = contract_pemdas(tokens,0).unwrap();
    assert_eq!(result, Token::Float(16.0));
  
  }
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_to_float_vector(){
  
    let tokens = vec![
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];

    let vec64 = to_float_vector(tokens,0).unwrap();
    assert_eq!(vec64, Token::VectorF64(vec![1.0,-1.0]));    


    let inv_sqrt2 = f64::sqrt(1.0/2.0);

    let tokens = vec![ 
      Token::SquareBracketOpen,
        Token::Float(1.0*inv_sqrt2), Token::Comma, Token::Float(-1.0*inv_sqrt2),
      Token::SquareBracketClose];

    let vec64 = to_float_vector(tokens,0).unwrap();
    assert_eq!(vec64, Token::VectorF64(vec![inv_sqrt2,-inv_sqrt2]));    



  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_contract_numeric_vector(){
  
    let tokens = vec![
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];

    let vec64 = contract_numeric_vectors(tokens,true,0).unwrap();
    assert_eq!(vec64, Token::VectorF64(vec![1.0,-1.0]));    


    let inv_sqrt2 = f64::sqrt(1.0/2.0);

    let tokens = vec![ 
      Token::Float(inv_sqrt2), Token::Times,
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];

    let vec64 = contract_numeric_vectors(tokens,true, 0).unwrap();
    assert_eq!(vec64, Token::VectorF64(vec![inv_sqrt2,-inv_sqrt2]));    


    let sqrt2 = f64::sqrt(2.0);
    let tokens = vec![ 
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose,
      Token::Slash, Token::Float(sqrt2),
    ];

    let vec64 = contract_numeric_vectors(tokens,true,0).unwrap();
    assert_eq!(vec64, Token::VectorF64(vec![1.0/sqrt2, -1.0/sqrt2]));    


    let tokens = vec![
      Token::Float(1.0), Token::Plus,
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];
    
    let vec64 = contract_numeric_vectors(tokens,true, 0).unwrap();
    assert_eq!(vec64, Token::VectorF64(vec![2.0, 0.0]));    


    let tokens = vec![
      Token::SquareBracketOpen,
        Token::Float(-1.0), Token::Comma, Token::Float(1.0) ,
      Token::SquareBracketClose,
      Token::Plus,
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];
    
    let vec64 = contract_numeric_vectors(tokens,true, 0).unwrap();
    assert_eq!(vec64, Token::VectorF64(vec![0.0, 0.0]));    

    let tokens = vec![Token::SquareBracketOpen, 
        Token::Int(1), Token::Comma, Token::Int(-1), 
      Token::SquareBracketClose];
    let veci32 = contract_numeric_vectors(tokens,false, 0).unwrap();
    assert_eq!(veci32, Token::VectorI32(vec![1, -1]));    
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_to_f64_token(){
    let tokens = vec![Token::SquareBracketOpen, 
        Token::UserInputValue("1".to_string()),Token::Comma, 
        Token::UserInputValue("-1".to_string()),
      Token::SquareBracketClose];

    let result = to_f64_token(tokens,0).unwrap();
    assert_eq!(result, Token::VectorF64(vec![1.0,-1.0]));    
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_to_i32_token(){
    let tokens = vec![Token::SquareBracketOpen, 
        Token::UserInputValue("1".to_string()),Token::Comma, 
        Token::UserInputValue("-1".to_string()),
      Token::SquareBracketClose];

    let result = to_i32_token(tokens,0).unwrap();
    assert_eq!(result, Token::VectorI32(vec![1,-1]));    
  }
  //----------------------------------------------------------------------------
 
}  
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
