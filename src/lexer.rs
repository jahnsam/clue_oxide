use crate::clue_errors::*;
use substring::Substring;

pub struct TokenStream{
  pub statements: Vec::<TokenExpression>,
  pub line_numbers: Vec::<usize>,
}


pub fn get_tokens_from_file(filename: &str) -> Result<TokenStream,CluEError>{

    let lexer = Lexer::new(filename)?;
    let tokens = parse_tokens(lexer)?;
    let tokens = find_comments(tokens);
    let (tokens,line_numbers) = prune_tokens(tokens)?;
    let statements = get_token_statements(tokens, &line_numbers)?;

    Ok(TokenStream{
        statements,
        line_numbers,
        })
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct TokenExpression{
  pub lhs: Vec::<Token>,
  pub rhs: Option<Vec::<Token>>,
  pub relationship: Option<Token>
}

impl TokenExpression{
  fn from(tokens: Vec::<Token>, line_number:usize) -> Result<Self,CluEError>{

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
      return Ok(TokenExpression{lhs, rhs: None,relationship: None})
    }

    let rhs = Some(rhs);
    let relationship = Some(tokens[idx].clone());
    Ok(TokenExpression{lhs, rhs, relationship})

  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


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
    "\n" | " " | "," |"[" | "]" | "{" | "}" |"(" | ")" | "<" | ">" | ";" 
      | "=" | "+" | "-" | "*" | "/" | "^" | "!" => true,
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

  tokens
}
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
#[derive(PartialEq, Debug, Clone)]
pub enum Token{
 Bang,
 BlockCommentEnd, 
 BlockCommentStart, 
 Comma,
 CurlyBracketClose,
 CurlyBracketOpen,
 Element,
 EOL,
 Equals,
 Float(f64),
 GreaterThan,
 GreaterThanEqualTo,
 Hat,
 In,
 Int(i32),
 LessThan,
 LessThanEqualTo,
 LineComment, 
 MagneticField,
 Minus,
 Not,
 NotEqual,
 NotIn,
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
      Token::Bang => "!".to_string(), 
      Token::BlockCommentEnd => "*/".to_string(), 
      Token::BlockCommentStart => "/*".to_string(), 
      Token::Comma => ",".to_string(),
      Token::CurlyBracketClose => "}".to_string(),
      Token::CurlyBracketOpen => "{".to_string(),
      Token::Element => "element".to_string(),
      Token::EOL => "\n".to_string(),
      Token::Equals => "=".to_string(),
      Token::Float(x) => x.to_string(),
      Token::GreaterThan => ">".to_string(),
      Token::GreaterThanEqualTo => ">=".to_string(),
      Token::Hat => "^".to_string(),
      Token::In => "in".to_string(),
      Token::Int(n) => n.to_string(),
      Token::LessThan => "<".to_string(),
      Token::LessThanEqualTo => "<=".to_string(),
      Token::LineComment => "//".to_string(), 
      Token::MagneticField => "magnetic_field".to_string(),
      Token::Minus => "-".to_string(),
      Token::Not => "not".to_string(),
      Token::NotEqual => "!=".to_string(),
      Token::NotIn => "not in".to_string(),
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
    "!" => Some(Token::Bang),
    "*/" => Some(Token::BlockCommentEnd),
    "/*" => Some(Token::BlockCommentStart),
    "," => Some(Token::Comma),
    "}" => Some(Token::CurlyBracketClose),
    "{" => Some(Token::CurlyBracketOpen),
    "element" => Some(Token::Element),
    "\n" => Some(Token::EOL),
    "=" => Some(Token::Equals),
    ">" => Some(Token::GreaterThan),
    ">=" => Some(Token::GreaterThanEqualTo),
    "in" => Some(Token::In),
    "^" => Some(Token::Hat),
    "<" => Some(Token::LessThan),
    "<=" => Some(Token::LessThanEqualTo),
    "//" => Some(Token::LineComment),
    "magnetic_field" => Some(Token::MagneticField),
    "-" => Some(Token::Minus),
    "not" => Some(Token::Not),
    "!=" => Some(Token::NotEqual),
    "not in" => Some(Token::NotIn),
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
fn is_relational_operator(token: &Token) -> bool{
  match token{
    Token::Equals => true,
    Token::GreaterThan => true,
    Token::GreaterThanEqualTo => true,
    Token::In => true,
    Token::LessThan => true,
    Token::LessThanEqualTo => true,
    Token::NotEqual => true,
    Token::NotIn => true,
    _ => false,
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
fn basic_token_algebraic_operations_f64(tokens: [Token;3],line_number: usize) 
 ->Result<Token,CluEError>{

  if tokens.len() != 3{
    return Err(CluEError::WrongVectorLength(line_number, 3, tokens.len()));
  }
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
    Token::Hat => Ok(Token::Float( a.powf(b) )),
    _ => Err(CluEError::NotAnOperator(line_number, tokens[1].to_string())),
  }

}
//------------------------------------------------------------------------------
fn contract_operator_inverse_operator_tokens(
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

      let new_token = basic_token_algebraic_operations_f64(
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
fn contract_exponentiation_tokens(tokens: Vec::<Token>,line_number: usize)
 -> Result<Vec::<Token>, CluEError>{

   contract_operator_inverse_operator_tokens(tokens, line_number,
       Token::Hat, Token::Hat)
 }
//------------------------------------------------------------------------------
fn contract_multiply_divide_tokens(tokens: Vec::<Token>,line_number: usize)
 -> Result<Vec::<Token>, CluEError>{

   contract_operator_inverse_operator_tokens(tokens, line_number,
       Token::Times, Token::Slash)
 }
//------------------------------------------------------------------------------
fn contract_add_subtract_tokens(tokens: Vec::<Token>,line_number: usize)
 -> Result<Vec::<Token>, CluEError>{

   contract_operator_inverse_operator_tokens(tokens, line_number,
       Token::Plus, Token::Minus)
 }
//------------------------------------------------------------------------------
fn add_token_at_index(add_this_token: Token, idx: usize, tokens: Vec::<Token>,
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
fn combine_staements(statements: Vec::<Vec::<Token>>) -> Vec::<Token>{

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
//------------------------------------------------------------------------------
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
fn count_token(target: &Token, tokens: &[Token]) -> usize{

  let mut counter = 0;

  for token in tokens.iter(){
    if token == target{
      counter += 1;
    }
  }
  counter   
}
//------------------------------------------------------------------------------
fn is_opening_delimiter(token: &Token) -> bool{
  match token{
    Token::ParenthesisOpen => true,
    Token::SquareBracketOpen => true,
    Token::CurlyBracketOpen => true,
    _ => false,
  }
}
//------------------------------------------------------------------------------
fn is_closing_delimiter(token: &Token) -> bool{
  match token{
    Token::ParenthesisClose => true,
    Token::SquareBracketClose => true,
    Token::CurlyBracketClose => true,
    _ => false,
  }
}
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
fn count_delimiter_pairs(tokens: &[Token], 
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
fn find_deepest_parentheses(tokens: &[Token], line_number: usize) 
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
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn get_token_statements(tokens: Vec::<Token>, line_numbers: &[usize]) 
 -> Result<Vec::<TokenExpression>,CluEError>
{
  let splits = split_on_token(tokens, Token::Semicolon);

  let mut statements = Vec::<TokenExpression>::with_capacity(splits.len());
  
  let mut counter = 0;
  for tokens in splits{
    let line_number = line_numbers[counter];
    counter += 1;
    let tok_exp = TokenExpression::from(tokens, line_number)?;
    statements.push(tok_exp);
  } 

  Ok(statements)

}
//------------------------------------------------------------------------------
fn split_on_token(tokens: Vec::<Token>, split_on: Token) 
 -> Vec::<Vec::<Token>>
{
    // Count statements.
    let num_statements = count_token(&split_on,&tokens);

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
fn find_token(target: &Token, tokens: &[Token]) -> Vec::<usize>{

  let mut out = Vec::<usize>::with_capacity(count_token(target,tokens) );

  for (ii, token) in tokens.iter().enumerate(){
    if *token == *target{
      out.push(ii);
    }
  }
  out
}
//------------------------------------------------------------------------------
/*
// TODO
fn to_string_vector(mut tokens: Vec::<Token>) 
 -> Result<Vec::<String>, CluEError>
{
  let vector_elements =  split_on_token(tokens, Token::Comma);
  let out = Vec::<String>::with_capacity(vector_elements.len());
}
*/
//------------------------------------------------------------------------------
// TODO
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
// TODO
fn to_float_vector(tokens: Vec::<Token>, line_number: usize) 
 -> Result<Vec::<f64>, CluEError>
{

  let indices = find_brackets(&tokens,line_number)?;

  let index0;
  let index1;

  if let Some((idx0,idx1)) = indices{
    index0 = idx0;
    index1 = idx1;
  }else{
    let token = contract_pemdas(tokens,line_number)?;
  
    if let Token::Float(x) = token{
  
      let mut out = Vec::<f64>::with_capacity(1);
      out.push(x);
      return Ok(out);
 
    }else{
      return Err(CluEError::CannotConvertToFloat(
            line_number,token.to_string()) );
    }
  }
  
  // Parse tthe scalar prefactor.
  let prefactor: f64;
  
  let token = contract_pemdas(tokens[0..index0-1].to_vec(), line_number)?;
  
  if let Token::Float(x) = token{
    prefactor = x;
  }else{
    return Err(CluEError::CannotConvertToFloat(line_number, token.to_string()));
  }


  // Parse the scalar post-factor.
  let postfactor: f64;
  
  let token = contract_pemdas(tokens[index1+2..tokens.len()].to_vec(),
      line_number)?;
  
  if let Token::Float(x) = token{
    postfactor = x;
  }else{
    return Err(CluEError::CannotConvertToFloat(line_number, token.to_string()));
  }



  // Get elements
  let vector_elements =  split_on_token(
      tokens[index0+1..index1].to_vec(), Token::Comma);


  // Initialize the output.
  let mut out = Vec::<f64>::with_capacity(vector_elements.len());

  // Assign output. 
  for token_element in vector_elements{

    // Contract the parentheses.
    let token = contract_emdas(token_element,line_number)?;
    
    if let Token::Float(x) = token{
      out.push(prefactor * x *postfactor);
    }else{
      return Err(CluEError::CannotConvertToFloat(line_number, 
            token.to_string()));
    }
  }
  
  Ok(out)

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_identify_token(){
    assert_eq!(identify_token("!").unwrap(), Token::Bang);
    assert_eq!(identify_token("*/").unwrap(), Token::BlockCommentEnd);
    assert_eq!(identify_token("/*").unwrap(), Token::BlockCommentStart);
    assert_eq!(identify_token(",").unwrap(), Token::Comma);
    assert_eq!(identify_token("}").unwrap(), Token::CurlyBracketClose);
    assert_eq!(identify_token("{").unwrap(), Token::CurlyBracketOpen);
    assert_eq!(identify_token("element").unwrap(), Token::Element);
    assert_eq!(identify_token("\n").unwrap(), Token::EOL);
    assert_eq!(identify_token("=").unwrap(), Token::Equals);
    assert_eq!(identify_token(">").unwrap(), Token::GreaterThan);
    assert_eq!(identify_token(">=").unwrap(), Token::GreaterThanEqualTo);
    assert_eq!(identify_token("^").unwrap(), Token::Hat);
    assert_eq!(identify_token("in").unwrap(), Token::In);
    assert_eq!(identify_token("<").unwrap(), Token::LessThan);
    assert_eq!(identify_token("<=").unwrap(), Token::LessThanEqualTo);
    assert_eq!(identify_token("//").unwrap(), Token::LineComment);
    assert_eq!(identify_token("magnetic_field").unwrap(), Token::MagneticField);
    assert_eq!(identify_token("-").unwrap(), Token::Minus);
    assert_eq!(identify_token("not").unwrap(), Token::Not);
    assert_eq!(identify_token("!=").unwrap(), Token::NotEqual);
    assert_eq!(identify_token("not in").unwrap(), Token::NotIn);
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
  fn test_basic_token_algebraic_operations_f64(){
  
    let mut tokens = vec![Token::Float(2.0), Token::Plus, Token::Float(3.0),
     Token::Equals, Token::Float(3.0)];
    
    let answer = basic_token_algebraic_operations_f64(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 + 3.0));

    tokens[1] = Token::Minus;
    let answer = basic_token_algebraic_operations_f64(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 - 3.0));

    tokens[1] = Token::Times;
    let answer = basic_token_algebraic_operations_f64(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 * 3.0));

    tokens[1] = Token::Slash;
    let answer = basic_token_algebraic_operations_f64(
        [tokens[0].clone(),tokens[1].clone(),tokens[2].clone()],0).unwrap();
    assert_eq!(answer, Token::Float(2.0 / 3.0));

    tokens[1] = Token::Hat;
    let answer = basic_token_algebraic_operations_f64(
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
    assert_eq!(vec64, vec![1.0,-1.0]);    


    let inv_sqrt2 = f64::sqrt(1.0/2.0);

    let tokens = vec![ 
      Token::Float(inv_sqrt2), Token::Times,
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];

    let vec64 = to_float_vector(tokens,0).unwrap();
    assert_eq!(vec64, vec![inv_sqrt2,-inv_sqrt2]);    


    let sqrt2 = f64::sqrt(2.0);
    let tokens = vec![ 
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose,
      Token::Slash, Token::Float(sqrt2),
    ];

    let vec64 = to_float_vector(tokens,0).unwrap();
    assert_eq!(vec64, vec![1.0/sqrt2, -1.0/sqrt2]);    


    // TODO: What should this be?
    let tokens = vec![
      Token::Float(1.0), Token::Plus,
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];
    
    assert!(false);


    // TODO: What should this be?
    let tokens = vec![
      Token::SquareBracketOpen,
        Token::Float(-1.0), Token::Comma, Token::Float(1.0) ,
      Token::SquareBracketClose,
      Token::Plus,
      Token::SquareBracketOpen,
        Token::Float(1.0), Token::Comma, Token::Float(-1.0) ,
      Token::SquareBracketClose];
    
    assert!(false);

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
