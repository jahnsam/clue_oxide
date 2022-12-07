use crate::clue_errors::*;
use substring::Substring;

struct Lexer{
  input: String,
  position: usize,
  line_number: usize,
}

impl Lexer{
  fn new(filename: &str) -> Result<Self,CluEError>{
    if let Ok(input) = std::fs::read_to_string(filename){
      return Ok(Lexer{input, position: 0, line_number: 1});
    }
    Err(CluEError::InvalidConfigFile(filename.to_string()))
  }
  //----------------------------------------------------------------------------
  fn next_token(&mut self) -> Token {

    let idx = self.position;
    
    if let Some(token) = identify_token(self.input.substring(idx,idx+1)){
      self.position += 1;
      return token;
    }

    for read_to in idx+1..self.input.len() {

      if !is_token_over(self.input.substring(read_to,read_to+1)){ continue; }
        
      self.position = read_to;

      if let Some(token) 
        = identify_token(self.input.substring(idx,read_to) ){
        return token;
      }

      return Token::UserInputValue( 
          Box::new(self.input.substring(idx, read_to).to_string()) );
      
    }

    self.position = self.input.len();

    Token::UserInputValue(Box::new(
          self.input.substring(idx,self.input.len() ).to_string()))

  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn is_token_over(word: &str) -> bool{
  match word{
    "\n" | " " | "," | "]" | "}" | ")" | ";" 
      | "=" | "+" | "-" | "*" | "/" => true,
      _ => false,
  }
}
//------------------------------------------------------------------------------
fn parse_tokens(mut lexer: Lexer) -> Result<Vec::<Token>, CluEError>{

  //let mut lexer = Lexer::new(filename)?;

  let mut tokens = Vec::<Token>::with_capacity(lexer.input.len());

  while lexer.position < lexer.input.len(){
    
    let token = lexer.next_token();
    
    if token == Token::EOL {
      lexer.line_number += 1;
    }

    tokens.push(token)
    
  }
  Ok(tokens)
}

//------------------------------------------------------------------------------
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
      was_last_token_times = true;
      out_tokens.push(Token::BlockCommentEnd);
    }
    else if was_last_token_slash{ // /R
      was_last_token_slash = false;
      out_tokens.push(Token::Slash);
      out_tokens.push(token);
    }
    else if was_last_token_times{ // *R
      was_last_token_times = true;
      out_tokens.push(Token::Times);
      out_tokens.push(token);
    }


  }

  out_tokens
}
//------------------------------------------------------------------------------
fn prune_tokens(in_tokens: Vec::<Token>) -> Result<Vec::<Token>,CluEError> {

  let mut out_tokens = Vec::<Token>::with_capacity(in_tokens.len());

  let mut block_commenting: i32 = 0;
  let mut is_line_commenting = false;

  let mut is_attribute =  false;

  let mut line_number = 1;
  for token in in_tokens{

    if token == Token::EOL{
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
    }

    if block_commenting > 0 {
      continue;
    }else if block_commenting < 0 {
      return Err(CluEError::UnmatchedBlockComment(line_number))
    }


    if token == Token::LineComment{
      is_line_commenting = true;
    
    }

    if token == Token::Whitespace{ continue; }
      
    if token == Token::Sharp{is_attribute = true;}
    out_tokens.push(token);
    
  }

  Ok(out_tokens)
}
//------------------------------------------------------------------------------
#[derive(PartialEq, Debug, Clone)]
enum Token{
 BlockCommentEnd, 
 BlockCommentStart, 
 Coma,
 CurlyBracketClose,
 CurlyBracketOpen,
 Element,
 EOL,
 Equals,
 In,
 LineComment, 
 ParenthesisClose,
 ParenthesisOpen,
 Residue,
 Semicolon,
 Sharp,
 Slash,
 SquareBracketClose,
 SquareBracketOpen,
 Times,
 TunnelSplitting,
 UserInputValue(Box<String>),
 Whitespace,
}

//------------------------------------------------------------------------------
fn identify_token(word: &str) -> Option<Token>{
  match word{
    //"*/" => Some(Token::BlockCommetEnd),
    //"/*" => Some(Token::BlockCommetStart),
    "," => Some(Token::Coma),
    "}" => Some(Token::CurlyBracketClose),
    "{" => Some(Token::CurlyBracketOpen),
    "element" => Some(Token::Element),
    "\n" => Some(Token::EOL),
    "=" => Some(Token::Equals),
    "in" => Some(Token::In),
    //"//" => Some(Token::LineComent),
    ")" => Some(Token::ParenthesisClose),
    "(" => Some(Token::ParenthesisOpen),
    "residue" => Some(Token::Residue),
    ";" => Some(Token::Semicolon),
    "#" => Some(Token::Sharp),
    "/" => Some(Token::Slash),
    "]" => Some(Token::SquareBracketClose),
    "[" => Some(Token::SquareBracketOpen),
    "*" => Some(Token::Times),
    "tunnel_splitting" => Some(Token::TunnelSplitting),
    " " => Some(Token::Whitespace),
    _ => None
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
struct TokenStatements{
  statements: Vec::<Vec::<Token>>,
  line_numbers: Vec::<usize>,
}

fn get_token_statements(tokens: Vec::<Token>) -> TokenStatements {

    // Count statements.
    let mut num_statements = 0;

    for token in tokens.iter(){
      if *token == Token::Semicolon { 
        num_statements += 1;
      }
    }


    // Get Statements.
    let mut statements = Vec::<Vec::<Token>>::with_capacity(num_statements);
    let mut line_numbers = Vec::<usize>::with_capacity(num_statements);

    for _istatement in 0..num_statements{
    
      // Count tokens in statement.
      let mut num_tokens_in_statement = 0;
      
      for token in tokens.iter(){
        if *token == Token::Semicolon { break;}
        num_tokens_in_statement += 1;
      }


      // Get statement.
      let mut statement = Vec::<Token>::with_capacity(num_tokens_in_statement);

      for token in tokens.iter(){
        if *token == Token::Semicolon { break;}
        statement.push(token.clone());
      }

      statements.push(statement);
    } 

    TokenStatements{ statements, line_numbers}
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_identify_token(){
  
    assert_eq!(identify_token("element").unwrap(), Token::Element);
    assert_eq!(identify_token(" ").unwrap(), Token::Whitespace);
    assert_eq!(identify_token("in").unwrap(), Token::In);
    assert_eq!(identify_token("[").unwrap(), Token::SquareBracketOpen);
    assert_eq!(identify_token("]").unwrap(), Token::SquareBracketClose);
    assert_eq!(identify_token(";").unwrap(), Token::Semicolon);
    assert_eq!(identify_token("\n").unwrap(), Token::EOL);
  }
  //----------------------------------------------------------------------------

  #[test]
  fn test_parse_tokens(){
  
    let lexer = Lexer{
      input: String::from(
"element in [H];
tunnel_splitting = 80e3; // Hz"),
      position: 0, 
      line_number: 0,
    };

    let ref_tokens = vec![Token::Element, Token::Whitespace, Token::In,
     Token::Whitespace, Token::SquareBracketOpen, 
     Token::UserInputValue(Box::new("H".to_string())), 
     Token::SquareBracketClose, Token::Semicolon, Token::EOL, 
     Token::TunnelSplitting, Token::Whitespace,
     Token::Equals, Token::Whitespace, 
     Token::UserInputValue(Box::new("80e3".to_string())), Token::Semicolon,
     Token::Whitespace, Token::Slash, Token::Slash, Token::Whitespace,
     Token::UserInputValue(Box::new("Hz".to_string())) ];


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
    };
    let ref_tokens = vec![ Token::BlockCommentStart,
     Token::Element, Token::Whitespace, Token::In,
     Token::Whitespace, Token::SquareBracketOpen, 
     Token::UserInputValue(Box::new("H".to_string())), 
     Token::SquareBracketClose, Token::Semicolon, 
     Token::BlockCommentEnd, Token::EOL, 
     Token::TunnelSplitting, Token::Whitespace,
     Token::Equals, Token::Whitespace, 
     Token::UserInputValue(Box::new("80e3".to_string())), Token::Semicolon,
     Token::Whitespace, Token::LineComment, Token::Whitespace,
     Token::UserInputValue(Box::new("Hz".to_string())) ];

    let tokens = parse_tokens(lexer).unwrap();
    let tokens = find_comments(tokens);


    assert_eq!(tokens.len(), ref_tokens.len());
   
    for ii in 0..tokens.len(){
      assert_eq!(tokens[ii], ref_tokens[ii]);
    }


  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
