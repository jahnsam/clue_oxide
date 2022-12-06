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
  fn next_token(&mut self) -> Result<Token, CluEError>{

    let err_token: String;

    for read_to in self.position..self.input.len() {

      if let Some(token) == identify_token(self.input[self.position..=read_to]){
        self.position = read_to + 1;
        return Ok(token);

      }
      else if is_token_over(self.input[self.read_to]){
        self.position = read_to;

        return Ok(
            Token::UserInputValue(
              Box(self.input[self.position .. read_to].clone()))
            );
      }
    }
    Err(CluEError::InvalidToken(self.line_number, err_token))
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn is_token_over(word: &str) -> bool{
  match word{
    "\n" | " " | "," | "]" | "}" | ")" | ";" => true,
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

fn prune_tokens(in_tokens: Vec::<Token>) -> Result<Vec::<Token>,CluEError> {

  let out_tokens = Vec::<Token>::with_capacity(in_tokens.len());

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
 SquareBracketClose,
 SquareBracketOpen,
 TunnelSplitting,
 UserInputValue(Box<String>),
 WhiteSpace,
}

//------------------------------------------------------------------------------
fn identify_token(word &str) -> Option<Token>{
  match word{
    "*\\" => Some(Token::BlockCommetEnd),
    "\\*" => Some(Token::BlockCommetStart),
    "," => Some(Token::Coma),
    "}" => Some(Token::CurlyBracketClose),
    "{" => Some(Token::CurlyBracketOpen),
    "element" => Some(Token::Element),
    "\n" => Some(Token::EOL),
    "=" => Some(Token::Equals),
    "in" => Some(Token::In),
    "//" => Some(Token::LineComent),
    ")" => Some(Token::ParenthesisClose),
    "(" => Some(Token::ParenthesisOpen),
    "residue" => Some(Token::Residue),
    ";" => Some(Token::Semicolon),
    "#" => Some(Token::Sharp),
    "]" => Some(Token::SquareBracketClose),
    "[" => Some(Token::SquareBracketOpen),
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




