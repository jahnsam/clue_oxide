use crate::clue_errors::*;
use crate::config::token::*;
use crate::config::token_stream::*;
use crate::config::token_expressions::*;
use substring::Substring;


//------------------------------------------------------------------------------
// This function reads in a file and returns a Vec::<TokenExpression>; 
// each element contains a line of input.
pub fn get_tokens_from_file(filename: &str) 
  -> Result<Vec::<TokenExpression>,CluEError>
{

    let lexer = Lexer::new(filename)?;
    let tokens = parse_tokens(lexer)?;
    let tokens = find_comments(tokens);
    let (tokens,line_numbers) = prune_tokens(tokens)?;
    let tokens =  simplify_tokens(tokens);
    let expressions = get_token_statements(tokens, &line_numbers)?;

    Ok(expressions)
}
//------------------------------------------------------------------------------
// This function reads in a file and returns a Vec::<TokenExpression>; 
// each element contains a line of input.
pub fn get_tokens_from_line(input: &str) 
  -> Result<Vec::<TokenExpression>,CluEError>
{

    let input = str::replace(input, ";", ";\n");
    let lexer = Lexer{input, position: 0, line_number: 0, 
      quoting: false};
    let tokens = parse_tokens(lexer)?;
    let tokens = find_comments(tokens);
    let (tokens,line_numbers) = prune_tokens(tokens)?;
    let tokens =  simplify_tokens(tokens);
    let expressions = get_token_statements(tokens, &line_numbers)?;

    Ok(expressions)
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
// TokenExpression holds a line of input, organized for easy parsing.
pub struct TokenExpression{
  pub lhs: Vec::<Token>,
  pub rhs: Option<Vec::<Token>>,
  pub relationship: Option<Token>,
  pub line_number: usize,
}

impl TokenExpression{
  //----------------------------------------------------------------------------
  // This function processes a vector of tokens into a TokenExpression.
  fn from(tokens: Vec::<Token>, line_number:usize) -> Result<Self,CluEError>{

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
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Lexer holds the string data the current reading position and line number.
struct Lexer{
  input: String,
  position: usize,
  line_number: usize,
  quoting: bool,
}

impl Lexer{
  //----------------------------------------------------------------------------
  // This function instantiates a ner lexer from a file.
  fn new(filename: &str) -> Result<Self,CluEError>{
    if let Ok(input) = std::fs::read_to_string(filename){
      return Ok(Lexer{input, position: 0, line_number: 1, 
          quoting: false,});
    }
    Err(CluEError::InvalidConfigFile(filename.to_string()))
  }
  //----------------------------------------------------------------------------
  // This function read the next several characters and returns the next
  // token.
  fn next_token(&mut self) -> Token {

    let idx = self.position;
    
    // Check for single character tokens.
    if let Some(token) = identify_token(self.input.substring(idx,idx+1)){
      if !self.quoting || token == Token::DoubleQuote{
        self.position += 1;
        return token;
      }
    }

    // Find multi-character tokens.
    for read_to in idx+1..self.input.len() {

      if !is_token_over(self.input.substring(read_to,read_to+1), self.quoting){ 
        continue; 
      }
        
      self.position = read_to;

      if let Some(token) 
        = identify_token(self.input.substring(idx,read_to) ){
        return token;
      }

      return Token::UserInputValue( 
          self.input.substring(idx, read_to).to_string() );
      
    }

    self.position = self.input.len();

    Token::UserInputValue(
        self.input.substring(idx,self.input.len() ).to_string())

  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// This function reports on whether the inpur character ends a token.
fn is_token_over(word: &str,is_within_quote: bool) -> bool{
  if is_within_quote{
    return word == "\"";
  }
  match word{
    "\n" | " " | "," |"[" | "]" | "{" | "}" |"(" | ")" | "<" | ">" | ";" 
      | "=" | "+" | "-" | "*" | "/" | "^" | "!" | "\"" => true,
      _ => false,
  }
}
//------------------------------------------------------------------------------
// This function consumes a lexer and retutns a vector of raw tokens.
fn parse_tokens(mut lexer: Lexer) -> Result<Vec::<Token>, CluEError>{

  //let mut lexer = Lexer::new(filename)?;

  let mut tokens = Vec::<Token>::with_capacity(lexer.input.len());

  while lexer.position < lexer.input.len(){
    
    let token = lexer.next_token();
    
    if token == Token::DoubleQuote{
      lexer.quoting = !lexer.quoting;
      continue;
    }else if token == Token::EOL {
      lexer.line_number += 1;
    }

    tokens.push(token)
    
  }
  Ok(tokens)
}

//------------------------------------------------------------------------------
// This function replaces the two token comments codes with a single
// token comment token.
fn find_comments(in_tokens: Vec::<Token>) -> Vec::<Token> {

  let mut out_tokens = Vec::<Token>::with_capacity(in_tokens.len());
  
  let mut was_last_token_slash = false;
  let mut was_last_token_times = false;

  for token in in_tokens{

    if !was_last_token_slash && !was_last_token_times {
      if token == Token::Slash{
        was_last_token_slash = true;
        continue;
      }
      else if token == Token::Times{
        was_last_token_times = true;
        continue
      }
    }

    if !was_last_token_slash && !was_last_token_times{ // R'R
      out_tokens.push(token);
    }
    else if was_last_token_slash && token == Token::Slash{ // //
      was_last_token_slash = false;
      out_tokens.push(Token::LineComment);
    }
    else if was_last_token_slash && token == Token::Times { // /*
      was_last_token_slash = false;
      out_tokens.push(Token::BlockCommentStart);
    }
    else if was_last_token_times && token == Token::Slash{ // */
      was_last_token_times = false;
      out_tokens.push(Token::BlockCommentEnd);
    }
    else if was_last_token_slash{ // /R
      was_last_token_slash = false;
      out_tokens.push(Token::Slash);
      out_tokens.push(token);
    }
    else if was_last_token_times{ // *R
      was_last_token_times = false;
      out_tokens.push(Token::Times);
      out_tokens.push(token);
    }


  }

  out_tokens
}
//------------------------------------------------------------------------------
// This function strips the comments from the input.
fn prune_tokens(in_tokens: Vec::<Token>) 
  -> Result<( Vec::<Token>, Vec::<usize>),CluEError> {

  let mut out_tokens = Vec::<Token>::with_capacity(in_tokens.len());
  let mut line_numbers = Vec::<usize>::with_capacity(in_tokens.len());

  let mut block_commenting: i32 = 0;
  let mut is_line_commenting = false;
  let mut is_newline = true;

  let mut is_attribute =  false;

  let mut line_number = 1;
  for token in in_tokens{

    if token == Token::EOL{
      is_newline = true;
      is_line_commenting = false;
      line_number += 1;
      if is_attribute{
        is_attribute = false;
        out_tokens.push(Token::Semicolon);
      }
      continue;
    }

    if token == Token::BlockCommentStart {  
      block_commenting += 1;
    }else if token == Token::BlockCommentEnd{
      block_commenting -= 1;
      continue
    }

    if block_commenting > 0 {
      continue;
    }else if block_commenting < 0 {
      return Err(CluEError::UnmatchedBlockComment(line_number))
    }


    if token == Token::LineComment{
      is_line_commenting = true; 
    }
    if is_line_commenting {continue;}

    if token == Token::Whitespace{ continue; }
      
    if token == Token::Sharp{is_attribute = true;}
    out_tokens.push(token);
    
    if is_newline {
      is_newline = false;
      line_numbers.push(line_number);
    }
    
  }

  Ok((out_tokens,line_numbers))
}
//------------------------------------------------------------------------------
// This function calls build_composit_tokens() to replace pre-defined
// token sequences with single tokens.  It repeats this until all composit
// tokens are found.
fn simplify_tokens(mut tokens: Vec::<Token>) -> Vec::<Token>
{
  let mut n = tokens.len();
  loop{
    tokens = build_composit_tokens(tokens);
    if tokens.len() == n{
      return tokens;
    }
    n = tokens.len();
  }  
}
//------------------------------------------------------------------------------
// This function takes a vector of tokens and identifies token sequences
// that form a new token.  The function returns the vector of tokens
// with the multi-token sequences replaced by a single token representing the
// composit token.
fn build_composit_tokens(mut tokens: Vec::<Token>) -> Vec::<Token>
{
  if tokens.len() < 2{
    return tokens;
  }

  let mut out = Vec::<Token>::with_capacity(tokens.len());

  // Add an extra token at the end to make the loop cleaner.
  tokens.push(Token::Whitespace);

  let mut skip_next = false;

  for ii in 1..tokens.len(){
  
    if skip_next{
      skip_next = false;
      continue;
    }
    skip_next = true;

    // !=
    if tokens[ii-1] == Token::Bang && tokens[ii] == Token::Equals{
      out.push(Token::NotEqual);
      continue;}

    // not in
    if tokens[ii-1] == Token::Not && tokens[ii] == Token::In{
      out.push(Token::NotIn);
      continue;}
  
    // /*
    if tokens[ii-1] == Token::Times && tokens[ii] == Token::Slash{
      out.push(Token::BlockCommentEnd);
      continue;}
  
    // */
    if tokens[ii-1] == Token::Slash && tokens[ii] == Token::Times{
      out.push(Token::BlockCommentStart);
      continue;}
  
    // >=
    if tokens[ii-1] == Token::GreaterThan && tokens[ii] == Token::Equals{
      out.push(Token::GreaterThanEqualTo);
      continue;}
  
    // -+
    if tokens[ii-1] == Token::Minus && tokens[ii] == Token::Plus{
      out.push(Token::Minus);
      continue;}
  
    // --
    if tokens[ii-1] == Token::Minus && tokens[ii] == Token::Minus{
      out.push(Token::Plus);
      continue;}
  
    // <=
    if tokens[ii-1] == Token::LessThan && tokens[ii] == Token::Equals{
      out.push(Token::LessThanEqualTo);
      continue;}
  
    // +-
    if tokens[ii-1] == Token::Plus && tokens[ii] == Token::Minus{
      out.push(Token::Minus);
      continue;}
  
    // ++
    if tokens[ii-1] == Token::Plus && tokens[ii] == Token::Plus{
      out.push(Token::Plus);
      continue;}
  
    // //
    if tokens[ii-1] == Token::Slash && tokens[ii] == Token::Slash{
      out.push(Token::LineComment);
      continue;}

    skip_next = false;
    out.push(tokens[ii-1].clone());
  }

  out
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
fn find_lhs_rhs_delimiter_index(tokens: &[Token], line_number: usize) 
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
*/

//------------------------------------------------------------------------------
// TODO: use or delete
/*
fn combine_statements(statements: Vec::<Vec::<Token>>) -> Vec::<Token>{

  let mut n_tokens = 0;
  for tokens in statements.iter(){
    n_tokens += tokens.len();
  }

  let mut out = Vec::<Token>::with_capacity(n_tokens);


  for tokens in statements.iter(){
    for token in tokens{
      out.push((*token).clone());
    }
  }

  out
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn get_token_statements(tokens: Vec::<Token>, line_numbers: &[usize]) 
 -> Result<Vec::<TokenExpression>,CluEError>
{
  let splits = split_on_token(tokens, Token::Semicolon);

  let mut statements = Vec::<TokenExpression>::with_capacity(splits.len());
  
  for (idx, tokens) in splits.into_iter().enumerate(){
    let line_number = line_numbers[idx];
    let tok_exp = TokenExpression::from(tokens, line_number)?;
    statements.push(tok_exp);
  } 

  Ok(statements)

}


//------------------------------------------------------------------------------
pub fn to_string_vector(tokens: Vec::<Token>, line_number: usize) 
 -> Result<Token, CluEError>
{

  let vector_elements = get_vector_elements(tokens,line_number)?;

  let mut out = Vec::<String>::with_capacity(vector_elements.len());

  for vec_el in vector_elements.iter(){
    let mut str_el = String::from("");
    for el in vec_el.iter(){
      str_el.push_str( &el.to_string() );
    }

    out.push(str_el);
  }

  Ok(Token::VectorString(out))
}



//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

 
  #[test]
  fn test_parse_tokens(){
  
    let lexer = Lexer{
      input: String::from(
"element in [H];
tunnel_splitting = 80e3; // Hz"),
      position: 0, 
      line_number: 0,
      quoting: false
    };

    let ref_tokens = vec![Token::Element, Token::Whitespace, Token::In,
     Token::Whitespace, Token::SquareBracketOpen, 
     Token::UserInputValue("H".to_string()), 
     Token::SquareBracketClose, Token::Semicolon, Token::EOL, 
     Token::TunnelSplitting, Token::Whitespace,
     Token::Equals, Token::Whitespace, 
     Token::UserInputValue("80e3".to_string()), Token::Semicolon,
     Token::Whitespace, Token::Slash, Token::Slash, Token::Whitespace,
     Token::UserInputValue("Hz".to_string()) ];


    let tokens = parse_tokens(lexer).unwrap();

    assert_eq!(tokens.len(), ref_tokens.len());
   
    for ii in 0..tokens.len(){
      assert_eq!(tokens[ii], ref_tokens[ii]);
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_find_comments(){
    let lexer = Lexer{
      input: String::from(
"/*element in [H];*/
tunnel_splitting = 80e3; // Hz"),
      position: 0, 
      line_number: 0,
      quoting: false
    };
    let ref_tokens = vec![ Token::BlockCommentStart,
     Token::Element, Token::Whitespace, Token::In,
     Token::Whitespace, Token::SquareBracketOpen, 
     Token::UserInputValue("H".to_string()), 
     Token::SquareBracketClose, Token::Semicolon, 
     Token::BlockCommentEnd, Token::EOL, 
     Token::TunnelSplitting, Token::Whitespace,
     Token::Equals, Token::Whitespace, 
     Token::UserInputValue("80e3".to_string()), Token::Semicolon,
     Token::Whitespace, Token::LineComment, Token::Whitespace,
     Token::UserInputValue("Hz".to_string()) ];

    let tokens = parse_tokens(lexer).unwrap();
    let tokens = find_comments(tokens);


    assert_eq!(tokens.len(), ref_tokens.len());
   
    for ii in 0..tokens.len(){
      assert_eq!(tokens[ii], ref_tokens[ii]);
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_prune_tokens(){
    let lexer = Lexer{
      input: String::from(
"/*element in [H];*/
tunnel_splitting = 80e3; // Hz
magnetic_field = 1.2; // T"),
      position: 0, 
      line_number: 0,
      quoting: false
    };
    let ref_tokens = vec![ 
     Token::TunnelSplitting, 
     Token::Equals, 
     Token::UserInputValue("80e3".to_string()), Token::Semicolon,
     Token::MagneticField, 
     Token::Equals, 
     Token::UserInputValue("1.2".to_string()), Token::Semicolon,
    ];

    let tokens = parse_tokens(lexer).unwrap();
    let tokens = find_comments(tokens);
    let (tokens,_line_numbers) = prune_tokens(tokens).unwrap();


    assert_eq!(tokens.len(), ref_tokens.len());
   
    for ii in 0..tokens.len(){
      assert_eq!(tokens[ii], ref_tokens[ii]);
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_token_statements(){
    let lexer = Lexer{
      input: String::from(
"/*element in [H];*/
tunnel_splitting = 80e3; // Hz
magnetic_field = 1.2; // T"),
      position: 0, 
      line_number: 0,
      quoting: false
    };
    let ref_statements = vec![ 
     vec![Token::TunnelSplitting, 
     Token::Equals, 
     Token::UserInputValue("80e3".to_string())],
     vec![Token::MagneticField, 
     Token::Equals, 
     Token::UserInputValue("1.2".to_string())],
    ];

    let ref_line_numbers = vec![2,3];
    let tokens = parse_tokens(lexer).unwrap();
    let tokens = find_comments(tokens);
    let (tokens,line_numbers) = prune_tokens(tokens).unwrap();
    let statements = get_token_statements(tokens,&line_numbers).unwrap();

    //let statements = &token_statements.statements;
    //let line_numbers = &token_statements.line_numbers;

    
    assert_eq!(statements.len(), ref_statements.len());
    assert_eq!(line_numbers.len(), ref_line_numbers.len());
   

    for jj in 0..statements.len(){
      assert_eq!(line_numbers[jj], ref_line_numbers[jj] );
      assert_eq!(statements[jj].lhs.len(), 1);
      
      assert_eq!(statements[jj].lhs[0], ref_statements[jj][0]);
      if let Some(rhs) = &statements[jj].rhs{
        assert_eq!(rhs[0], ref_statements[jj][2]);
      }else{
      panic!("Expected Some(rhs), but found None.");
    }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_simplify_tokens(){
  
    let tokens = vec![Token::Float(1.0), Token::Plus, Token::Minus,
     Token::Float(2.0), Token::Minus, Token::Minus, Token::Minus, 
     Token::Float(5.0), Token::Bang, Token::Equals, Token::Float(3.0)];

    let ref_tokens = vec![Token::Float(1.0), Token::Minus, 
        Token::Float(2.0), Token::Minus, Token::Float(5.0), Token::NotEqual,
       Token::Float(3.0)];

    let tokens = simplify_tokens(tokens);

    assert_eq!(tokens.len(), ref_tokens.len());

    for ii in 0..tokens.len(){
      assert_eq!(tokens[ii], ref_tokens[ii]);
    }
  }
  //----------------------------------------------------------------------------
  /*
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
  */
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_to_string_vector(){

    let tokens = vec![Token::UserInputValue("SOL".to_string())]; 
    
    let token = to_string_vector(tokens,0).unwrap();
    
    assert_eq!(token, 
        Token::VectorString(vec!["SOL".to_string()])
        );



    let tokens = vec![
      Token::SquareBracketOpen,
        Token::UserInputValue("SOL".to_string()), Token::Comma, 
        Token::UserInputValue("WAT".to_string()) ,
      Token::SquareBracketClose];
    
    let token = to_string_vector(tokens,0).unwrap();

    assert_eq!(token, 
        Token::VectorString(vec!["SOL".to_string(),"WAT".to_string()])
        );



    let tokens = vec![
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];

    let token = to_string_vector(tokens,0).unwrap();

    assert_eq!(token, 
        Token::VectorString(vec![1.0.to_string(),(-1.0).to_string()])
        );

  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
