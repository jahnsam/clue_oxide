use crate::clue_errors::*;
use crate::config::token::*;

//------------------------------------------------------------------------------
fn are_delimiters_paired(open: &Token, close: &Token, line_number: usize) 
 -> Result<bool,CluEError>{

  if !is_opening_delimiter(open){
    return Err(CluEError::InvalidToken(line_number, open.to_string()));
  }
  if !is_closing_delimiter(close){
    return Err(CluEError::InvalidToken(line_number, close.to_string()));
  } 

  if (*open == Token::ParenthesisOpen && *close == Token::ParenthesisClose)  
    || (*open == Token::SquareBracketOpen 
        && *close == Token::SquareBracketClose)
    || (*open == Token::CurlyBracketOpen && *close == Token::CurlyBracketClose){
      return Ok(true);
    }
  Ok(false)
}
 
//------------------------------------------------------------------------------
fn get_delimiter_depths(tokens: &[Token], line_number: usize)
  -> Result<(Vec::<i32>, i32), CluEError>{

  let mut depths = Vec::<i32>::with_capacity(tokens.len()); 
  let mut depth: i32 = 0; 
  let mut max_depth: i32 = 0;

  for token in tokens.iter(){
    if is_opening_delimiter(token){
      depth += 1;
    } else if is_closing_delimiter(token){
      depth -= 1;
    }
     
    if depth < 0{
      return Err(CluEError::UnmatchedDelimiter(line_number) ) 
    }

    max_depth = std::cmp::max(max_depth, depth);

    depths.push(depth);
  }

  if depth != 0{
    return Err(CluEError::UnmatchedDelimiter(line_number) ) 
  }

  Ok((depths, max_depth))
} 
//------------------------------------------------------------------------------
pub fn split_on_token(tokens: Vec::<Token>, split_on: Token) 
 -> Vec::<Vec::<Token>>
{
    // Count statements.
    let mut num_statements = count_token(&split_on,&tokens);
    if tokens[tokens.len()-1] != split_on{
      num_statements += 1;
    }

    // Get Statements.
    let mut statements = Vec::<Vec::<Token>>::with_capacity(num_statements);

    let mut read_from = 0;
    for _istatement in 0..num_statements{
    
      // Count tokens in statement.
      let mut num_tokens_in_statement = 0;
      
      let mut read_to = tokens.len();
      for ii in read_from..tokens.len(){
        if tokens[ii] == split_on { 
          read_to = ii;
          break;
        }
        num_tokens_in_statement += 1;
      }


      // Get statement.
      let mut statement = Vec::<Token>::with_capacity(num_tokens_in_statement);

      for ii in read_from..read_to{
        statement.push(tokens[ii].clone());
      }

      read_from = read_to + 1;
      statements.push(statement);
    } 

    statements
}

//------------------------------------------------------------------------------
fn find_brackets(tokens: &[Token], line_number: usize) 
 -> Result<Option<(usize,usize)>, CluEError>
{
  // Find brackets.

  let n_pairs = count_delimiter_pairs(&tokens, 
      &Token::SquareBracketOpen, &Token::SquareBracketClose,
      line_number)?;
  
  if n_pairs > 1{
    return Err(CluEError::CannotConvertToVector(line_number));
  }

  if n_pairs == 0{
    return Ok(None); 
  }


  let index0 = find_token(&Token::SquareBracketOpen, &tokens);
  let index0 = index0[0];

  let index1 = find_token(&Token::SquareBracketClose, &tokens);
  let index1 = index1[0];


  if index0 + 1 >= index1{
    return Err(CluEError::CannotConvertToFloat(
          line_number,"[]".to_string()) );
  }
 
  Ok(Some( (index0,index1) ))
}
//------------------------------------------------------------------------------
pub fn read_strings_as_floats(tokens: Vec::<Token>, line_number: usize)
-> Result<Vec::<Token>, CluEError>
{

  // Initialized output.
  let mut out = Vec::<Token>::with_capacity(tokens.len());

  for token in tokens.iter(){
    
    match token{
    
      Token::UserInputValue(val) => {
    
        match val.parse(){
          Ok(x) => out.push(Token::Float(x)),
          Err(_) => return Err(CluEError::CannotConvertToFloat(
                line_number,val.to_string() )),
        }
      }, 
    
      _ => out.push((*token).clone() )
    }
  }

  Ok(out)
  
}
//------------------------------------------------------------------------------
pub fn read_strings_as_integers(tokens: Vec::<Token>, line_number: usize)
-> Result<Vec::<Token>, CluEError>
{

  // Initialized output.
  let mut out = Vec::<Token>::with_capacity(tokens.len());

  for token in tokens.iter(){
    
    match token{
    
      Token::UserInputValue(val) => {
    
        match val.parse(){
          Ok(a) => out.push(Token::Int(a)),
          Err(_) => return Err(CluEError::CannotConvertToFloat(
                line_number,val.to_string() )),
        }
      }, 
    
      _ => out.push((*token).clone() )
    }
  }

  Ok(out)
  
}
//------------------------------------------------------------------------------
pub fn get_vector_elements(tokens: Vec::<Token>, line_number: usize)
 -> Result<Vec::<Vec::<Token>>, CluEError>
{
  let indices = find_brackets(&tokens,line_number)?;


  let vector_elements: Vec::<Vec::<Token>>;

  if let Some((idx0,idx1)) = indices{
    vector_elements =  split_on_token(
      tokens[idx0+1..idx1].to_vec(), Token::Comma);
  }else{
    let mut temp =  Vec::<Vec::<Token>>::with_capacity(1);
    temp.push(tokens);
    vector_elements = temp;
  }

  Ok(vector_elements)
}
//------------------------------------------------------------------------------
pub fn add_token_at_index(add_this_token: Token, idx: usize, tokens: Vec::<Token>,
    line_number: usize)
  -> Result<Vec::<Token>,CluEError>
{
  if idx > tokens.len(){
    return Err(CluEError::IndexOutOfBounds(line_number, idx,tokens.len() +1));
  }

  let mut out = Vec::<Token>::with_capacity(tokens.len() +1);
  
  for (ii,token) in tokens.iter().enumerate(){
    if ii == idx{
      out.push(add_this_token.clone());
    }
    out.push((*token).clone());
  }

  Ok(out)
}  
//------------------------------------------------------------------------------
pub fn count_delimiter_pairs(tokens: &[Token], 
    opening_delimiter: &Token, closing_delimiter: &Token, 
    line_number: usize)
  -> Result<usize, CluEError>{
  
  let n_open = count_token(&opening_delimiter, tokens);  
  let n_close = count_token(&closing_delimiter, tokens);  
  


  if n_open != n_close
  {
    return Err(CluEError::UnmatchedDelimiter(line_number) ) 
  }

  Ok(n_open)
}

//------------------------------------------------------------------------------
pub fn find_deepest_parentheses(tokens: &[Token], line_number: usize) 
  -> Result<Option<(usize,usize)>, CluEError>{

  
  let n_pairs = count_delimiter_pairs(tokens,
      &Token::ParenthesisOpen, &Token::ParenthesisClose,
      line_number)?;
  
  if n_pairs == 0{
    return Ok(None);
  }


  let (depths, max_depth) = get_delimiter_depths(tokens,line_number)?;


  let mut found_open_of_max_depth = false; 
  let mut found_close_of_max_depth = false; 
  let mut idx_open = 0;
  let mut idx_close = tokens.len() - 1;

  for (ii, depth) in depths.iter().enumerate(){
    if !found_open_of_max_depth{
      if *depth == max_depth{
        idx_open = ii;
        found_open_of_max_depth = true;
      }
    } else if !found_close_of_max_depth && *depth < max_depth{
      idx_close = ii;
      found_close_of_max_depth = true;
    }  
  }

  if !found_open_of_max_depth || !found_close_of_max_depth{
    return Err(CluEError::UnmatchedDelimiter(line_number) ) 
  }

  let are_matched = are_delimiters_paired(
      &tokens[idx_open], &tokens[idx_close],line_number)?;
  
  if !are_matched{
    return Err(CluEError::UnmatchedDelimiter(line_number) ) 
  }

  Ok( Some( (idx_open,idx_close) ) )

}





//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  
  #[test]
  fn test_read_strings_as_floats(){
    let tokens = vec![Token::SquareBracketOpen, 
        Token::UserInputValue("1".to_string()),Token::Comma, 
        Token::UserInputValue("-1".to_string()),
      Token::SquareBracketClose];

    let result = read_strings_as_floats(tokens,0).unwrap();
    assert_eq!(result, vec![Token::SquareBracketOpen, 
        Token::Float(1.0),Token::Comma, Token::Float(-1.0),
        Token::SquareBracketClose]
        );
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_read_strings_as_integers(){
    let tokens = vec![Token::SquareBracketOpen, 
        Token::UserInputValue("1".to_string()),Token::Comma, 
        Token::UserInputValue("-1".to_string()),
      Token::SquareBracketClose];

    let result = read_strings_as_integers(tokens,0).unwrap();
    assert_eq!(result, vec![Token::SquareBracketOpen, 
        Token::Int(1),Token::Comma, Token::Int(-1),
        Token::SquareBracketClose]
        );
  }
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_split_on_token(){
  
    let tokens = vec![
      Token::SquareBracketOpen,
        Token::UserInputValue("SOL".to_string()), Token::Comma, 
        Token::UserInputValue("WAT".to_string()), Token::Comma, 
        Token::UserInputValue("GLY".to_string()) ,
      Token::SquareBracketClose];

    let result = split_on_token(tokens,Token::Comma);


    assert_eq!(result, vec![
      vec![Token::SquareBracketOpen, Token::UserInputValue("SOL".to_string())],
      vec![Token::UserInputValue("WAT".to_string())],
      vec![Token::UserInputValue("GLY".to_string()),Token::SquareBracketClose]  
    ]);

    let tokens = vec![
        Token::Float(2.0), Token::Comma, Token::Float(-2.0) 
    ];

    let result = split_on_token(tokens,Token::Comma);

    assert_eq!(result, vec![
        vec![Token::Float(2.0)],
        vec![Token::Float(-2.0)]  ]);

  }
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_find_deepest_parentheses(){
    
    let tokens = vec![Token::ParenthesisClose, Token::ParenthesisOpen];
    let result = find_deepest_parentheses(&tokens,0);
    assert_eq!(result,Err(CluEError::UnmatchedDelimiter(0)));

    let tokens = vec![Token::ParenthesisOpen, Token::SquareBracketClose];
    let result = find_deepest_parentheses(&tokens,0);
    assert_eq!(result,Err(CluEError::UnmatchedDelimiter(0)));

    let tokens = vec![Token::Comma, Token::Equals];
    let option = find_deepest_parentheses(&tokens,0).unwrap();
    assert_eq!(option,None);

    let tokens = vec![Token::ParenthesisOpen, Token::ParenthesisClose];
    let (idx_open, idx_close) = find_deepest_parentheses(&tokens,0)
      .unwrap().unwrap();
    assert_eq!(idx_open,0);
    assert_eq!(idx_close,1);

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
    
    let option = find_deepest_parentheses(&tokens,0).unwrap();
    assert_eq!(option,Some((1,4)));
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
