use crate::clue_errors::CluEError;

use std::io::BufRead;
use std::io::BufWriter;
use std::collections::HashMap;
use substring::Substring;
use std::fs::File;
use std::io::prelude::*;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn read_cif_file(filename: &str)
  -> Result< (), CluEError>
{
  let lexer = CIFLexer::new(filename)?;
  let data = form_cif_dictionary(lexer)?;
  Ok(())
}
//------------------------------------------------------------------------------
fn form_cif_dictionary(lexer: CIFLexer) 
  -> Result<HashMap::<String,Vec::<String>>,CluEError>
{
  let mut data = HashMap::<String,Vec::<String>>::new();
  let mut keys = Vec::<String>::new();
  let mut loop_index = 0;
  let mut mode = CIFMode::Open;

  for (ii,token) in lexer.tokens.iter().enumerate(){
    let line_number = lexer.line_numbers[ii];

    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if token == "loop_"{
      if mode != CIFMode::Open && mode != CIFMode::LoopDataEntry{  
        return Err(CluEError::CIFCannotStartLoop(line_number));
      }
      mode = CIFMode::LoopItemDef;
      keys = Vec::<String>::new();
    }
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    else if token.substring(0,1) == "_"{ 
      if data.get(token).is_some(){
        return Err(CluEError::CIFDataItemAlreadyDeclared(
          line_number,token.to_string()));
      }

      match mode{
        CIFMode::Open => {
          keys = vec![token.to_string()];
          data.insert(token.to_string(),Vec::<String>::new());
          mode = CIFMode::DataEntry;
        },
        CIFMode::LoopItemDef => {
          keys.push(token.to_string());
          data.insert(token.to_string(),Vec::<String>::new());
        },
        CIFMode::LoopDataEntry => {
          keys = vec![token.to_string()];
          data.insert(token.to_string(),Vec::<String>::new());
          mode = CIFMode::DataEntry;
        },
        _ => return Err(CluEError::CIFCannotDeclareDataItem(
              line_number,token.to_string())),
      }
    }else if mode == CIFMode::DataEntry{

      if keys.len() != 1 {
        return Err(CluEError::CIFNoDeclaredDataItem(
              line_number,token.to_string()));
      }

      let Some(data_vec) = data.get_mut(&keys[0]) else{
        return Err(CluEError::CIFNoDeclaredDataItem(
              line_number,token.to_string()));
      };

      data_vec.push(token.to_string());
      keys = Vec::<String>::new();
      mode = CIFMode::Open;
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }else if mode == CIFMode::LoopDataEntry 
      || mode == CIFMode::LoopDataEntryLocked
      || mode == CIFMode::LoopItemDef {
      if keys.is_empty() {
        return Err(CluEError::CIFNoDeclaredDataItem(
              line_number,token.to_string()));
      }

      let Some(data_vec) = data.get_mut(&keys[loop_index]) else{
        return Err(CluEError::CIFNoDeclaredDataItem(
              line_number,token.to_string()));
      };

      data_vec.push(token.to_string());

      loop_index = (loop_index + 1) % keys.len();
      
      if loop_index == 0 {
        mode = CIFMode::LoopDataEntry
      }else{
        mode = CIFMode::LoopDataEntryLocked;
      }


    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
  }

  Ok(data)
}
//------------------------------------------------------------------------------
#[derive(PartialEq,Debug)]
enum CIFMode{
  Open,
  DataEntry,
  LoopItemDef,
  LoopDataEntryLocked,
  LoopDataEntry,
}
//------------------------------------------------------------------------------
/*
pub fn parse_cif(filename: &str)
  -> Result< Structure, CluEError>
{

  let bath_particles = parse_cif_atoms()?

  Ok( Structure::new(
    bath_particles,
    connections,
    cell_offsets)
    )
}
//------------------------------------------------------------------------------
struct ParticleMap{
  particles: Vec::<Particle>,
  map: HashMap::<u32,Option<usize>>
}
//------------------------------------------------------------------------------
fn parse_cif_atoms(filename: &str,n_atoms: usize)
  -> Result<ParticleMap,CluEError>
{
  let mut bath_particles = Vec::<Particle>::with_capacity(n_atoms);

  let Ok(file) = std::fs::File::open(filename) else{
    return Err(CluEError::CannotOpenFile(filename.to_string()));
  };

  let lines = std::io::BufReader::new(file).lines();
  for line_result in lines {
    let Ok(line) = line_result else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };
  }
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Lexer holds the string data the current reading position and line number.
struct CIFLexer{
  pub tokens: Vec::<String>,
  pub line_numbers: Vec::<usize>,
}

impl CIFLexer{
  
  //----------------------------------------------------------------------------
  // This function instantiates a ner lexer from a file.
  fn new(filename: &str) -> Result<Self,CluEError>{

    let Ok(file) = std::fs::File::open(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };

    let lines = std::io::BufReader::new(file).lines();


    let mut tokens = Vec::<String>::new();
    let mut line_numbers = Vec::<usize>::new();

    let mut line_number = 0;
    for line_result in lines {
      line_number += 1;

      let Ok(line) = line_result else{
        return Err(CluEError::CannotOpenFile(filename.to_string()));
      };

      for word in line
        .replace("\'"," \' ")
        .replace("\""," \" ")
        .replace("#"," # ")
        .split_whitespace(){
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if word == "#"{
          break;
        }
        tokens.push(word.to_string());
        line_numbers.push(line_number);
      }
    }

    (tokens, line_numbers) = collapse_quotes(tokens,";", &line_numbers);
    (tokens, line_numbers) = collapse_quotes(tokens,"\'", &line_numbers);
    (tokens, line_numbers) = collapse_quotes(tokens,"\"", &line_numbers);
    Ok(CIFLexer{tokens, line_numbers})
  }

}
//------------------------------------------------------------------------------
fn collapse_quotes(tokens: Vec::<String>,quote_char: &str,
    line_numbers: &[usize]) 
  -> (Vec::<String>, Vec::<usize>)
{
  let mut out = Vec::<String>::with_capacity(tokens.len());
  let mut out_ln = Vec::<usize>::with_capacity(tokens.len());

  let mut is_quoting = false;
  let mut next_token = String::new();
  
  let mut idx = -1;
  let mut line_num = 0;
  for token in tokens{
    idx += 1;

    if token == quote_char{
      is_quoting = !is_quoting;
    }else{

    if next_token.is_empty(){
      next_token = token;
      line_num = line_numbers[idx as usize];
    }else{
      next_token = format!("{} {}",next_token,token);
    }
    }
    if is_quoting {continue; }
    
    out.push(next_token.clone());
    out_ln.push(line_num);
    next_token = String::new();  
  }
  (out, out_ln)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_read_cif_file(){
    let filename = "assets/VO(TCPP-Zn2-bpy).cif";
    read_cif_file(&filename).unwrap();
    assert!(false);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_form_cif_dictionary(){
    let filename = "assets/VO(TCPP-Zn2-bpy).cif";
    let lexer = CIFLexer::new(filename).unwrap();
    let data = form_cif_dictionary(lexer).unwrap();
    
    assert_eq!(data["_cell_length_a"],vec![String::from("16.644(2)")]);
    
    let ref_atom_site_label: Vec::<String> = vec![
      "C1", "C2", "H1", "C3", "H2", 
      "C4", "C5", "C6", "C7", "H3", 
      "C8", "H4", "C9", "H5", "C10", 
      "H6", "C11", "C12", "C13", "H13", 
      "C14", "H14", "C15", "C16", "C17", 
      "H17", "C18", "H18", "N1", "N2", 
      "N3", "O1", "O2", "O3", "O4", 
      "Zn1", "Zn2", "V1", "V2"
    ].iter().map(|s| s.to_string()).collect();
    assert_eq!(data["_atom_site_label"],ref_atom_site_label);

    let ref_symmetry_equiv_pos_as_xyz = vec![
      "x, y, z",
      "-x, -y, z",
      "-y, x, z",
      "y, -x, z",
    ].iter().map(|s| s.to_string()).collect::<Vec::<String>>();
    assert_eq!(data["_symmetry_equiv_pos_as_xyz"],
        ref_symmetry_equiv_pos_as_xyz);
  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_CIFLexer(){
    let filename = "assets/VO(TCPP-Zn2-bpy).cif";
    let lexer = CIFLexer::new(filename).unwrap();
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
