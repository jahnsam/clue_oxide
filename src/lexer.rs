struct Lexer{
  input: String,
  position: usize,
  line_number: usize,
}

impl Lexer{
  fn new(filename: &str) -> Result<Self,CluEError>{
    if let Ok(input) = read_file(filename){
      return Ok(Lexer{input, position: 0, line_number: 1});
    }
    Err(CluEError::InvalidConfigFile(filename.clone()))
  };
  //----------------------------------------------------------------------------
  fn next_token(&mut self) -> Token {

    let idx = self.position;
    
    for read_to in idx..self.input.len() {

      if !is_token_over(self.input[read_to]){ continue; }
        
      self.position = read_to;

      if let Some(token) = identify_token(self.input[idx..=read_to]){
        return token;
      }

      return Token::UserInputValue( Box(self.input[idx .. read_to].clone()) );
      
    }

    self.position = self.input.len();

    Token::UserInputValue(Box(self.input[idx ..].clone()))

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
fn parse_token(filename: &str){

  let mut lexer = Lexer::new(filename);

  let mut tokens = Vec::<Token>::with_capacity(lexer.input.len());

  while lexer.position < lexer.input.len(){
    if let Ok(token) = lexer.next_token(){
      if token == Token::EOL {
        lexer.line_number += 1;
      }
      tokens.push(token)
    }
  }
}

//------------------------------------------------------------------------------
fn find_comments(in_tokens: Vec::<Token>) -> Result<Vec::<Token>,CluEError> {

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

    if !was_last_token_slash && !was_last_token_times{
      out_tokens.push(token);
    }
    else if was_last_token_slash && token == Token::Slash{
      was_last_token_slash = false;
      out_tokens.push(Token::LineComment);
    }
    else if was_last_token_slash && token == Token::Times
      was_last_token_slash = false;
      out_tokens.push(Token::BlockCommentStart);
    }
    else if was_last_token_times && token == Token::Slash{
      was_last_token_times = true;
      out_tokens.push(Token::BlockCommentEnd);
    }
    else if was_last_token_slash{
      was_last_token_slash = false;
      out_tokens.push(Token::Slash);
      out_tokens.push(token);
    }
    else if was_last_token_times{
      was_last_token_times = true;
      out_tokens.push(Token::times);
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
        tokens.push(Token::Semicolon);
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

    if token == Token::WhiteSpace{ continue; }
      
    if token == Token::Sharp{is_attribute = true;}
    out_tokens.push(token);
    
  }

  Ok(out_tokens)
}
//------------------------------------------------------------------------------
enum Token{
 BlockCommetEnd, 
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
 WhiteSpace,
}

//------------------------------------------------------------------------------
fn identify_token(word &str) -> Option<Token>{
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
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
struct TokenStatements{
  statements: Vec::<Vec::<Token>>,
  line_number: Vec::<usize>,
}

fn get_token_statements(tokens: Vec::<Token>) -> TokenStatements {
  Vec::<Vec::<Token>>

    // Count statements.
    let mut num_statements = 0;

    for token in tokens.iter(){
      if token == Token::Semicolon { 
        num_statements += 1;
      }
    }


    // Get Statements.
    let mut statements = Vec::<Vec::<Token>>::with_capacity(num_statements);

    for _istatement in 0..num_statements{
    
      // Count tokens in statement.
      let mut num_tokens_in_statement = 0;
      
      for token in tokens.iter(){
        if token == Token::Semicolon { break;}
        num_tokens_in_statement += 1;
      }


      // Get statement.
      let mut statement = Vec::<Token>::with_capacity(num_tokens_in_statement);

      for token in tokens.iter(){
        if token == Token::Semicolon { break;}
        statement.push(token);
      }

      statements.push(statement);
    } 

    statements
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




