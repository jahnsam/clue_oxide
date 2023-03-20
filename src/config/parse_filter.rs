use crate::Config;
use crate::config::Token;
use crate::config::token_stream;
use crate::config::particle_config::ParticleConfig;
use crate::config::token_expressions::*;
use crate::clue_errors::CluEError;
use crate::structure::particle_filter::ParticleFilter;
use crate::config::to_i32_token;

impl Config{
  pub fn parse_filter_line(&mut self, expression: &TokenExpression,
      label_opt: &Option<String>)
    -> Result<(),CluEError>
  {
    let Some(label) = label_opt else{
      return Err(CluEError::MissingFilterLabel(expression.line_number));
    };

    let mut found_particle_config: Option<usize> = None;
    for (idx, p_cfg) in self.particles.iter().enumerate(){
      if *label == p_cfg.label{
        found_particle_config = Some(idx);
        break;
      }
    }

    let particle_config: &mut ParticleConfig; 
    match found_particle_config{
      Some(idx) => particle_config = &mut self.particles[idx],
      None => {
        self.particles.push(ParticleConfig::new((*label).clone()));
        let idx = self.particles.len() - 1; 
        particle_config = &mut self.particles[idx];}
    }

    if particle_config.filter == None{
      particle_config.filter = Some(ParticleFilter::new());
    }
    let Some(filter) = &mut particle_config.filter else{
      return Err(CluEError::MissingFilter( (*label).clone() ));
    };

    let include: bool;
    match expression.relationship{
      Some(Token::In) | Some(Token::Equals) => include = true,
      Some(Token::NotIn) | Some(Token::NotEqual) => include = false,
      _ => return Err(CluEError::NoRelationalOperators(expression.line_number)),
    }
   
    let tokens = token_stream::extract_rhs(expression)?;

    let already_set = ||{
      CluEError::OptionAlreadySet(
          expression.line_number, expression.lhs[0].to_string()) };
    match expression.lhs[0]{
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Indices => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          if !vec.is_empty(){return Err(already_set());}
          let vec = vec.iter().map(|&x| x as usize).collect();
          if include{
            filter.indices = vec;
          }else{
            filter.not_indices = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Elemments => (),
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Serials => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          if !vec.is_empty(){return Err(already_set());}
          let vec = vec.iter().map(|&x| x as u32).collect();
          if include{
            filter.serials = vec;
          }else{
            filter.not_serials = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Residues => (),
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      /*
      Token::ResSeqNums => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          if !vec.is_empty(){return Err(already_set());}
          let vec = vec.iter().map(|&x| x as u32).collect();
          if include{
            filter.indices = vec;
          }else{
            filter.not_indices = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      */
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::BondedIndices => {
        let value_token = to_i32_token(tokens, expression.line_number)?;
        if let Token::VectorI32(vec) = value_token{
          if !vec.is_empty(){return Err(already_set());}
          let vec = vec.iter().map(|&x| x as usize).collect();
          if include{
            filter.bonded_indices = vec;
          }else{
            filter.not_bonded_indices = vec;
          }
        }else{
          return Err(CluEError::NoRHS(expression.line_number));
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::BondedElements => (),
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
    Ok(())

  }
}
