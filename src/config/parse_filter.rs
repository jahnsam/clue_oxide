use crate::clue_errors::CluEError;
use crate::Config;
use crate::config::Token;
use crate::config::token_stream;
use crate::config::particle_config::ParticleConfig;
use crate::config::token_expressions::*;
use crate::config::to_i32_token;
use crate::physical_constants::ANGSTROM;
use crate::structure::particle_filter::ParticleFilter;

impl Config{
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



    // Get relational symbol.
    let include: bool = match expression.relationship{
      Some(Token::In) | Some(Token::Equals) 
        | Some(Token::LessThanEqualTo) => true,
      Some(Token::NotIn) | Some(Token::NotEqual) 
        | Some(Token::GreaterThanEqualTo) => false,
      _ => return Err(CluEError::NoRelationalOperators(expression.line_number)),
    };
   


    let tokens = token_stream::extract_rhs(expression)?;

    let already_set = ||{
      CluEError::OptionAlreadySet(
          expression.line_number, expression.lhs[0].to_string()) };
    match expression.lhs[0]{
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Indices => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as usize).collect();
          if include{
            if !filter.indices.is_empty(){return Err(already_set());}
            filter.indices = vec;
          }else{
            if !filter.not_indices.is_empty(){return Err(already_set());}
            filter.not_indices = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::CellIDs => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as usize).collect();
          if include{
            if !filter.cell_ids.is_empty(){return Err(already_set());}
            filter.cell_ids = vec;
          }else{
            if !filter.not_cell_ids.is_empty(){return Err(already_set());}
            filter.not_cell_ids = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Elements => {
        let value_token = vec_tokens_to_vec_elements(tokens);
        if let Ok(vec) = value_token{
          if include{
            if !filter.elements.is_empty(){return Err(already_set());}
            filter.elements = vec;
          }else{
            if !filter.not_elements.is_empty(){return Err(already_set());}
            filter.not_elements = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Distance => {
        if include{
          set_to_some_f64(&mut filter.within_distance, &expression)?;
          if let Some(r) = &mut filter.within_distance{
            *r *= ANGSTROM;
          }else{ return Err(CluEError::FilterNoMaxDistance(label.to_string()));}
        }else{
          set_to_some_f64(&mut filter.not_within_distance, &expression)?;
          if let Some(r) = &mut filter.not_within_distance{
            *r *= ANGSTROM;
          }else{ return Err(CluEError::FilterNoMinDistance(label.to_string()));}
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Serials => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as u32).collect();
          if include{
            if !filter.serials.is_empty(){return Err(already_set());}
            filter.serials = vec;
          }else{
            if !filter.not_serials.is_empty(){return Err(already_set());}
            filter.not_serials = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Residues => {
        let value_token = vec_tokens_to_vec_strings(tokens);
        if let Ok(vec) = value_token{
          if include{
            if !filter.residues.is_empty(){return Err(already_set());}
            filter.residues = vec;
          }else{
            if !filter.not_residues.is_empty(){return Err(already_set());}
            filter.not_residues = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      Token::ResSeqNums => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as u32).collect();
          if include{
            if !filter.residue_sequence_numbers.is_empty(){
              return Err(already_set());
            }
            filter.residue_sequence_numbers = vec;
          }else{
            if !filter.not_residue_sequence_numbers.is_empty(){
              return Err(already_set());
            }
            filter.not_residue_sequence_numbers = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::BondedIndices => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          let vec = vec.iter().map(|&x| x as usize).collect();
          if include{
            if !filter.bonded_indices.is_empty(){return Err(already_set());}
            filter.bonded_indices = vec;
          }else{
            if !filter.not_bonded_indices.is_empty(){return Err(already_set());}
            filter.not_bonded_indices = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::BondedElements => {
        let value_token = vec_tokens_to_vec_elements(tokens);
        if let Ok(vec) = value_token{
          if include{
            if !filter.bonded_elements.is_empty(){return Err(already_set());}
            filter.bonded_elements = vec;
          }else{
            if !filter.not_bonded_elements.is_empty(){return Err(already_set());}
            filter.not_bonded_elements = vec;
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
        "serials in [17,18]; serials not in [28,29];
        indices = [1,2,3]; indices not in [6];
        bonded_indices in [5,4]; bonded_indices not in [42];
        residue_sequence_numbers in [1501]; residue_sequence_numbers != [9];
        elements in [H,O]; elements not in [N];
        bonded_elements != [C]; bonded_elements in [N];
        residues in SOL; residues not in TEM;
        distance <= 4;
        distance >= 1;")
      .unwrap();

    let mut config = Config::new();
 
    for expression in expressions.iter(){
      config.parse_filter_line(expression,&Some("filter_label".to_string()))
        .unwrap();
    }

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
