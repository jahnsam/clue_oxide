use crate::clue_errors::CluEError;
use crate::Config;
use crate::config::Token;
use crate::config::token_stream;
use crate::config::particle_config::{CellType,ParticleConfig};
use crate::config::token_expressions::*;
use crate::config::to_i32_token;
use crate::physical_constants::ANGSTROM;
use crate::structure::particle_filter::ParticleFilter;

impl Config{
  ///  This function parses lines found under "#[group(...)]"
  pub fn parse_filter_line(&mut self, expression: &TokenExpression,
      label_opt: &Option<String>)
    -> Result<(),CluEError>
  {
    // Get label.
    let Some(label) = label_opt else{
      return Err(CluEError::MissingFilterLabel(expression.line_number));
    };

    // Find the correct particle_config.
    let mut found_particle_config: Option<usize> = None;
    for (idx, p_cfg) in self.particles.iter().enumerate(){
      if *label == p_cfg.label{
        found_particle_config = Some(idx);
        break;
      }
    }

    let particle_config: &mut ParticleConfig = match found_particle_config{
      Some(idx) => &mut self.particles[idx],
      None => {
        self.particles.push(ParticleConfig::new((*label).clone()));
        let idx = self.particles.len() - 1; 
        &mut self.particles[idx]}
    };

    if particle_config.filter.is_none(){
      particle_config.filter = Some(ParticleFilter::new());
    }
    let Some(filter) = &mut particle_config.filter else{
      return Err(CluEError::MissingFilter( (*label).clone() ));
    };


    let already_set = ||{
      CluEError::OptionAlreadySet(
          expression.line_number, expression.lhs[0].to_string()) };
    match expression.lhs[0]{
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Indices => {
        let tokens = token_stream::extract_rhs(expression)?;
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as usize).collect();
          match expression.relationship{
            Some(Token::In) => {
              if !filter.indices.is_empty(){return Err(already_set());}
              filter.indices = vec;
            },
            Some(Token::NotIn) => {
              if !filter.not_indices.is_empty(){return Err(already_set());}
              filter.not_indices = vec;
            },
            _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Elements => {
        let tokens = token_stream::extract_rhs(expression)?;
        let value_token = vec_tokens_to_vec_elements(tokens);
        if let Ok(vec) = value_token{
          match expression.relationship{
            Some(Token::In) => {
              if !filter.elements.is_empty(){return Err(already_set());}
              filter.elements = vec;
            },
            Some(Token::NotIn) => {
              if !filter.not_elements.is_empty(){return Err(already_set());}
              filter.not_elements = vec;
            },
            _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Extracells => {
        if particle_config.cell_type != CellType::AllCells{
          return Err(already_set());
        }
        particle_config.cell_type = CellType::Extracells;
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Distance => {
        match expression.relationship { 
          Some(Token::LessThanEqualTo) => {
            set_to_some_f64(&mut filter.within_distance, expression)?;
            if let Some(r) = &mut filter.within_distance{
              *r *= ANGSTROM;
            }else{ 
              return Err(CluEError::FilterNoMaxDistance(label.to_string()));
            }
          },
          Some(Token::GreaterThanEqualTo) => {
            set_to_some_f64(&mut filter.not_within_distance, expression)?;
            if let Some(r) = &mut filter.not_within_distance{
              *r *= ANGSTROM;
            }else{ return Err(CluEError::FilterNoMinDistance(label.to_string()));}
          },
          _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::PrimaryCell => {
        if particle_config.cell_type != CellType::AllCells{
          return Err(already_set());
        }
        particle_config.cell_type = CellType::PrimaryCell;
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Serials => {
        let tokens = token_stream::extract_rhs(expression)?;
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as u32).collect();
          match expression.relationship{
            Some(Token::In) => {
              if !filter.serials.is_empty(){return Err(already_set());}
              filter.serials = vec;
            },
            Some(Token::NotIn) => {
              if !filter.not_serials.is_empty(){return Err(already_set());}
              filter.not_serials = vec;
            },
            _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Residues => {
        let tokens = token_stream::extract_rhs(expression)?;
        let value_token = vec_tokens_to_vec_strings(tokens);
        if let Ok(vec) = value_token{
          match expression.relationship{
            Some(Token::In) => {
              if !filter.residues.is_empty(){return Err(already_set());}
              filter.residues = vec;
            },
            Some(Token::NotIn) => {
              if !filter.not_residues.is_empty(){return Err(already_set());}
              filter.not_residues = vec;
            },
            _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      Token::ResSeqNums => {
        let tokens = token_stream::extract_rhs(expression)?;
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as u32).collect();
          match expression.relationship{
            Some(Token::In) => {
              if !filter.residue_sequence_numbers.is_empty(){
                return Err(already_set());
              }
              filter.residue_sequence_numbers = vec;
            },
            Some(Token::NotIn) => {
              if !filter.not_residue_sequence_numbers.is_empty(){
                return Err(already_set());
              }
              filter.not_residue_sequence_numbers = vec;
            },
            _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::BondedIndices => {
        let tokens = token_stream::extract_rhs(expression)?;
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as usize).collect();
          match expression.relationship{
            Some(Token::In) => {
              if !filter.bonded_indices.is_empty(){return Err(already_set());}
              filter.bonded_indices = vec;
            },
            Some(Token::NotIn) => {
              if !filter.not_bonded_indices.is_empty(){return Err(already_set());}
              filter.not_bonded_indices = vec;
            },
            _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::BondedElements => {
        let tokens = token_stream::extract_rhs(expression)?;
        let value_token = vec_tokens_to_vec_elements(tokens);
        if let Ok(vec) = value_token{
          match expression.relationship{
            Some(Token::In) => {
              if !filter.bonded_elements.is_empty(){return Err(already_set());}
              filter.bonded_elements = vec;
            },
            Some(Token::NotIn) => {
              if !filter.not_bonded_elements.is_empty(){return Err(already_set());}
              filter.not_bonded_elements = vec;
            },
            _ => return Err(CluEError::NoRelationalOperators(
                expression.line_number)),
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
    Ok(())

  }
}



#[cfg(test)]
mod tests{
  use super::*;
  use crate::config::get_tokens_from_line;
  use crate::physical_constants::Element;
  #[test]
  fn test_parse_filter_line(){
    let expressions = get_tokens_from_line(
        "extracells;
        serials in [17,18]; serials not in [28,29];
        indices in [1,2,3]; indices not in [6];
        bonded_indices in [5,4]; bonded_indices not in [42];
        residue_sequence_numbers in [1501]; 
        residue_sequence_numbers not in [9];
        elements in [H,O]; elements not in [N];
        bonded_elements not in [C]; bonded_elements in [N];
        residues in SOL; residues not in TEM;
        distance <= 4;
        distance >= 1;")
      .unwrap();

    let mut config = Config::new();
 
    for expression in expressions.iter(){
      config.parse_filter_line(expression,&Some("filter_label".to_string()))
        .unwrap();
    }

    assert_eq!(config.particles[0].cell_type,CellType::Extracells);

    let filter = config.particles[0].filter.as_ref().unwrap();
    
    assert_eq!(filter.serials,vec![17,18]);
    assert_eq!(filter.not_serials,vec![28,29]);
    assert_eq!(filter.indices,vec![1,2,3]);
    assert_eq!(filter.not_indices,vec![6]);
    assert_eq!(filter.bonded_indices,vec![5,4]);
    assert_eq!(filter.not_bonded_indices,vec![42]);
    assert_eq!(filter.residue_sequence_numbers,vec![1501]);
    assert_eq!(filter.not_residue_sequence_numbers,vec![9]);
    assert_eq!(filter.elements,vec![Element::Hydrogen,Element::Oxygen]);
    assert_eq!(filter.not_elements,vec![Element::Nitrogen]);
    assert_eq!(filter.bonded_elements,vec![Element::Nitrogen]);
    assert_eq!(filter.not_bonded_elements,vec![Element::Carbon]);
    assert_eq!(filter.residues,vec!["SOL".to_string()]);
    assert_eq!(filter.not_residues,vec!["TEM".to_string()]);
    assert_eq!(filter.within_distance,Some(4e-10));
    assert_eq!(filter.not_within_distance,Some(1e-10));


  }
}
