use crate::clue_errors::*;
use substring::Substring;

pub struct TokenStream{
  statements: Vec::<Vec::<Token>>,
  line_numbers: Vec::<usize>,
}


pub fn get_tokens_from_file(filename: &str) -> Result<TokenStream,CluEError>{

    let lexer = Lexer::new(filename)?;
    let tokens = parse_tokens(lexer)?;
    let tokens = find_comments(tokens);
    let (tokens,line_numbers) = prune_tokens(tokens)?;
    let statements = get_token_statements(tokens);

    Ok(TokenStream{
        statements,
        line_numbers,
        })
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
          self.input.substring(idx, read_to).to_string() );
      
    }

    self.position = self.input.len();

    Token::UserInputValue(
        self.input.substring(idx,self.input.len() ).to_string())

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
#[derive(PartialEq, Debug, Clone)]
pub enum Token{
 BlockCommentEnd, 
 BlockCommentStart, 
 Comma,
 CurlyBracketClose,
 CurlyBracketOpen,
 Element,
 EOL,
 Equals,
 Float(f64),
 In,
 Int(i32),
 LineComment, 
 MagneticField,
 Minus,
 ParenthesisClose,
 ParenthesisOpen,
 Plus,
 Residue,
 Semicolon,
 Sharp,
 Slash,
 SquareBracketClose,
 SquareBracketOpen,
 Times,
 TunnelSplitting,
 UserInputValue(String),
 Whitespace,
}
impl Token{
  pub fn to_string(&self) -> String{
    match self{
      Token::BlockCommentEnd => "*/".to_string(), 
      Token::BlockCommentStart => "/*".to_string(), 
      Token::Comma => ",".to_string(),
      Token::CurlyBracketClose => "}".to_string(),
      Token::CurlyBracketOpen => "{".to_string(),
      Token::Element => "element".to_string(),
      Token::EOL => "\n".to_string(),
      Token::Equals => "=".to_string(),
      Token::Float(x) => x.to_string(),
      Token::In => "in".to_string(),
      Token::Int(n) => n.to_string(),
      Token::LineComment => "//".to_string(), 
      Token::MagneticField => "magnetic_field".to_string(),
      Token::Minus => "-".to_string(),
      Token::ParenthesisClose => ")".to_string(),
      Token::ParenthesisOpen => "(".to_string(),
      Token::Plus => "+".to_string(),
      Token::Residue => "residue".to_string(),
      Token::Semicolon => ";".to_string(),
      Token::Sharp => "#".to_string(),
      Token::Slash => "/".to_string(),
      Token::SquareBracketClose => "]".to_string(),
      Token::SquareBracketOpen => "[".to_string(),
      Token::Times => "*".to_string(),
      Token::TunnelSplitting => "tunnel_splitting".to_string(),
      Token::UserInputValue(string) => (*string).clone(),
      Token::Whitespace => " ".to_string(),

    }
  }
}
//------------------------------------------------------------------------------
fn identify_token(word: &str) -> Option<Token>{
  match word{
    "*/" => Some(Token::BlockCommentEnd),
    "/*" => Some(Token::BlockCommentStart),
    "," => Some(Token::Comma),
    "}" => Some(Token::CurlyBracketClose),
    "{" => Some(Token::CurlyBracketOpen),
    "element" => Some(Token::Element),
    "\n" => Some(Token::EOL),
    "=" => Some(Token::Equals),
    "in" => Some(Token::In),
    "//" => Some(Token::LineComment),
    "magnetic_field" => Some(Token::MagneticField),
    "-" => Some(Token::Minus),
    ")" => Some(Token::ParenthesisClose),
    "(" => Some(Token::ParenthesisOpen),
    "+" => Some(Token::Plus),
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
fn find_lhs_rhs_delimiter_index(tokens: &[Token], line_number: usize) 
  -> Result<usize,CluEError>{
  let mut found_token = false;
  let mut index: usize = 0;

  for ii in 0..tokens.len(){
    if tokens[ii] == Token::Equals || tokens[ii] == Token::In{
      if !found_token{
        found_token = true;
        index = ii;
      }else{
        return Err(CluEError::TooManyRelationalOperators(line_number));
      }
    }
  }
  if found_token{
    return Ok(index);
  }else{
    return Err(CluEError::NoRelationalOperators(line_number));
  }
}
//------------------------------------------------------------------------------
struct TokenExpression{
  lhs: Vec::<Token>,
  rhs: Vec::<Token>,
  relationship: Token
}
impl TokenExpression{
  fn from(tokens: Vec::<Token>, line_number:usize) -> Result<Self,CluEError>{

    let idx = find_lhs_rhs_delimiter_index(&tokens, line_number)?;

    let mut lhs = Vec::<Token>::with_capacity(idx);
    let mut rhs = Vec::<Token>::with_capacity(tokens.len() - idx - 1);
    let relationship = tokens[idx].clone();

    for (ii,token) in tokens.iter().enumerate(){

      if ii < idx{
        lhs.push(token.clone());
      }else if ii > idx{
        rhs.push(token.clone());
      }
    }
    Ok(TokenExpression{lhs, rhs, relationship})

  }
}
//------------------------------------------------------------------------------

fn is_numeric(token: &Token) -> bool{
  match token{
    Token::Float(_x) => true,
    Token::Int(_n) => true,
    _ => false

  }
}
//------------------------------------------------------------------------------
fn basic_token_algebraic_operations_f64(tokens: &[Token;3],line_number: usize) 
 ->Result<Token,CluEError>{

  let a: f64;
  match tokens[0]{
    Token::Float(x) => a = x,
    Token::Int(n) => a = n as f64,
    _ => return 
      Err(CluEError::CannotConvertToFloat(line_number, tokens[0].to_string()))
  }

  let b: f64;
  match tokens[2]{
    Token::Float(x) => b = x,
    Token::Int(n) => b = n as f64,
    _ => return 
      Err(CluEError::CannotConvertToFloat(line_number, tokens[2].to_string()))
  }

  match tokens[1] {
    Token::Plus  => Ok(Token::Float(a+b)),
    Token::Minus => Ok(Token::Float(a-b)),
    Token::Times => Ok(Token::Float(a*b)),
    Token::Slash => Ok(Token::Float(a/b)),
    _ => Err(CluEError::NotAnOperator(line_number, tokens[1].to_string())),
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn get_token_statements(tokens: Vec::<Token>) -> Vec::<Vec::<Token>> {

    // Count statements.
    let mut num_statements = 0;

    for token in tokens.iter(){
      if *token == Token::Semicolon { 
        num_statements += 1;
      }
    }


    // Get Statements.
    let mut statements = Vec::<Vec::<Token>>::with_capacity(num_statements);
    //let mut line_numbers = Vec::<usize>::with_capacity(num_statements);

    let mut read_from = 0;
    for _istatement in 0..num_statements{
    
      // Count tokens in statement.
      let mut num_tokens_in_statement = 0;
      
      let mut read_to = tokens.len();
      for ii in read_from..tokens.len(){
        if tokens[ii] == Token::Semicolon { 
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

    //TokenStatements{ statements, line_numbers}
    statements
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_identify_token(){
    assert_eq!(identify_token("*/").unwrap(), Token::BlockCommentEnd);
    assert_eq!(identify_token("/*").unwrap(), Token::BlockCommentStart);
    assert_eq!(identify_token(",").unwrap(), Token::Comma);
    assert_eq!(identify_token("}").unwrap(), Token::CurlyBracketClose);
    assert_eq!(identify_token("{").unwrap(), Token::CurlyBracketOpen);
    assert_eq!(identify_token("element").unwrap(), Token::Element);
    assert_eq!(identify_token("\n").unwrap(), Token::EOL);
    assert_eq!(identify_token("=").unwrap(), Token::Equals);
    assert_eq!(identify_token("in").unwrap(), Token::In);
    assert_eq!(identify_token("//").unwrap(), Token::LineComment);
    assert_eq!(identify_token("magnetic_field").unwrap(), Token::MagneticField);
    assert_eq!(identify_token("-").unwrap(), Token::Minus);
    assert_eq!(identify_token(")").unwrap(), Token::ParenthesisClose);
    assert_eq!(identify_token("(").unwrap(), Token::ParenthesisOpen);
    assert_eq!(identify_token("+").unwrap(), Token::Plus);
    assert_eq!(identify_token("residue").unwrap(), Token::Residue);
    assert_eq!(identify_token(";").unwrap(), Token::Semicolon);
    assert_eq!(identify_token("#").unwrap(), Token::Sharp);
    assert_eq!(identify_token("/").unwrap(), Token::Slash);
    assert_eq!(identify_token("]").unwrap(), Token::SquareBracketClose);
    assert_eq!(identify_token("[").unwrap(), Token::SquareBracketOpen);
    assert_eq!(identify_token("*").unwrap(), Token::Times);
    assert_eq!(identify_token("tunnel_splitting").unwrap(), 
        Token::TunnelSplitting);
    assert_eq!(identify_token(" ").unwrap(), Token::Whitespace);

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
    let (tokens,line_numbers) = prune_tokens(tokens).unwrap();


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
    let statements = get_token_statements(tokens);

    //let statements = &token_statements.statements;
    //let line_numbers = &token_statements.line_numbers;

    
    assert_eq!(statements.len(), ref_statements.len());
    assert_eq!(line_numbers.len(), ref_line_numbers.len());
   
    println!("\n\n{:?}\n\n",statements);

    for jj in 0..statements.len(){
      assert_eq!(line_numbers[jj], ref_line_numbers[jj] );
      assert_eq!(statements[jj].len(), ref_statements[jj].len());
      for ii in 0..statements[jj].len(){
        assert_eq!(statements[jj][ii], ref_statements[jj][ii]);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_find_lhs_rhs_delimiter_index(){
  
    let mut tokens = vec![Token::Float(1.0), Token::Plus, Token::Float(2.0),
    Token::Comma, Token::Float(3.0)];

    let result = find_lhs_rhs_delimiter_index(&tokens,0);
    assert_eq!(result,Err(CluEError::NoRelationalOperators(0)));

    tokens[3] = Token::Equals;
    let result = find_lhs_rhs_delimiter_index(&tokens,0);
    assert_eq!(result,Ok(3));

    tokens[0] = Token::In;
    let result = find_lhs_rhs_delimiter_index(&tokens,0);
    assert_eq!(result,Err(CluEError::TooManyRelationalOperators(0)));
  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_TokenExpression(){
    let tokens = vec![Token::Float(1.0), Token::Plus, Token::Float(2.0),
    Token::Equals, Token::Float(3.0)];

    let expression = TokenExpression::from(tokens.clone(),0).unwrap();

    assert_eq!(expression.relationship,tokens[3]);
    for ii in 0..3{
      assert_eq!(expression.lhs[ii],tokens[ii]);
    }
    assert_eq!(expression.rhs[0],tokens[4]);
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
