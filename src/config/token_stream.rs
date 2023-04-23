use crate::clue_errors::*;
use crate::config::token::*;
use crate::config::token_expressions::*;
use substring::Substring;

//------------------------------------------------------------------------------
// The function takes a Vec::<Token> and inserts a token at the specified index.
pub fn add_token_at_index(add_this_token: Token, idx: usize, 
    tokens: Vec::<Token>, line_number: usize)
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
// This function unwraps TokenExpression right hand side or returns an
// appropriate error message.
pub fn extract_rhs(expression: &TokenExpression)
  -> Result<Vec::<Token>,CluEError>
{

  let tokens: Vec::<Token>;
  if let Some(toks) = &expression.rhs{
    tokens = toks.clone();
  }else{
    return Err(CluEError::NoRHS(expression.line_number));
  }

  Ok(tokens)
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
// This function counts the number of delimiter pairs.  
// It will err if the opening and closing delimiters to not have the same
// count, but no type checking is done.
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
pub fn find_outermost_parentheses(tokens: &[Token], line_number: usize) 
  -> Result<Option<(usize,usize)>, CluEError>
{

  let n_pairs = count_delimiter_pairs(tokens,
      &Token::ParenthesisOpen, &Token::ParenthesisClose,
      line_number)?;
  
  if n_pairs == 0{
    return Ok(None);
  }

  let (depths, _max_depth) = get_delimiter_depths(tokens,line_number)?;

  let mut found_open_of_depth_1 = false; 
  let mut found_close_of_depth_1 = false; 
  let mut idx_open = 0;
  let mut idx_close = tokens.len() - 1;

  for (ii, depth) in depths.iter().enumerate(){
    if !found_open_of_depth_1{
      if *depth == 1{
        idx_open = ii;
        found_open_of_depth_1 = true;
      }
    } else if !found_close_of_depth_1 && *depth < 1{
      idx_close = ii;
      found_close_of_depth_1 = true;
    }  
  }

  if !found_open_of_depth_1 || !found_close_of_depth_1{
    return Err(CluEError::UnmatchedDelimiter(line_number) ) 
  }

  let are_matched = are_delimiters_paired(
      &tokens[idx_open], &tokens[idx_close],line_number)?;
  
  if !are_matched{
    return Err(CluEError::UnmatchedDelimiter(line_number) ) 
  }

  Ok( Some( (idx_open,idx_close) ) )
}
//------------------------------------------------------------------------------
// This function finds the first pair of delimiters that do not contain
// other delimiters.
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
//------------------------------------------------------------------------------
// This function takes a Vec::<Token> and replaces each 
// Token::UserInputValue(String) with a Token::Float(f64).
// The function will err if the conversion fails.
pub fn read_strings_as_floats(tokens: Vec::<Token>, line_number: usize)
-> Result<Vec::<Token>, CluEError>
{

  // Initialized output.
  let mut out = Vec::<Token>::with_capacity(tokens.len());

  for token in tokens.iter(){
    
    match token{
    
      Token::UserInputValue(val0) => {
    
        let mut val = val0.to_string();

        let push_time_ten_hat: bool;
        if val.len() >= 2 && val.substring( val.len()-1,val.len())=="e"{
          val = val.substring(0, val.len()-1).to_string();
          push_time_ten_hat = true;
        }else{
          push_time_ten_hat = false;
        }
        match val.parse(){
          Ok(x) => out.push(Token::Float(x)),
          Err(_) => return Err(CluEError::CannotConvertToFloat(
                line_number,val.to_string() )),
        }
        if push_time_ten_hat{
          out.push(Token::Times);
          out.push(Token::Float(10.0));
          out.push(Token::Hat);
        }
      }, 
    
      _ => out.push((*token).clone() )
    }
  }

  Ok(out)
  
}
//------------------------------------------------------------------------------
// This function takes a Vec::<Token> and replaces each 
// Token::UserInputValue(String) with a Token::Int(i32).
// The function will err if the conversion fails.
pub fn read_strings_as_integers(tokens: Vec::<Token>, line_number: usize)
-> Result<Vec::<Token>, CluEError>
{

  // Initialized output.
  let mut out = Vec::<Token>::with_capacity(tokens.len());

  for token in tokens.iter(){
    
    match token{
    
      Token::UserInputValue(val0) => {
        let mut val = val0.to_string();

        let push_time_ten_hat: bool;
        if val.len() >= 2 && val.substring( val.len()-1,val.len())=="e"{
          val = val.substring(0, val.len()-1).to_string();
          push_time_ten_hat = true;
        }else{
          push_time_ten_hat = false;
        }
    
        match val.parse(){
          Ok(a) => out.push(Token::Int(a)),
          Err(_) => return Err(CluEError::CannotConvertToFloat(
                line_number,val.to_string() )),
        }
        if push_time_ten_hat{
          out.push(Token::Times);
          out.push(Token::Int(10));
          out.push(Token::Hat);
        }
      }, 
    
      _ => out.push((*token).clone() )
    }
  }

  Ok(out)
  
}
//------------------------------------------------------------------------------
// This function takes a vector of tokens and return a vector of vector of
// tokens, split on the specified token.  
// The specified token in removed from the output.
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


//------------------------------------------------------------------------------
// This function compares two tokens and checks that they form are an
// opening delimier and closing delimiter of the same kind,
// such as ( and ), [ and ], or { and }, 
// but not { and ] or ) and ).
fn are_delimiters_paired(open: &Token, close: &Token, line_number: usize) 
 -> Result<bool,CluEError>{

  // Ckeck if token is (, [, or {.
  if !is_opening_delimiter(open){
    return Err(CluEError::InvalidToken(line_number, open.to_string()));
  }

  // Ckeck if token is ), ], or }.
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
/// The functions take a token slice and return the index of the first "[" 
/// and the following  "]".  The finction will err if there are nested brackets.
pub fn find_brackets(tokens: &[Token], line_number: usize) 
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
// This function assigns a delimiter depth to each token.
// For example, let t be an unspecified token, then the depths of
// t t ( t t ( t ) t ( ) t ) t are
// 0 0 1 1 1 2 2 1 1 2 2 1 0 0.
// The function will err if the dept is ever negative.
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
  
