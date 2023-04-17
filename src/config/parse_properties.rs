use crate::clue_errors::CluEError;
use crate::Config;
use crate::config::Token;
use crate::config::token_stream;
use crate::config::particle_config::ParticleConfig;
use crate::config::token_expressions::*;
use crate::config::particle_config::ParticleProperties;

impl Config{
  pub fn parse_properties_line(&mut self, expression: &TokenExpression,
      label_opt: &Option<String>)
    -> Result<(),CluEError>
  {
    // Get label.
    let Some(label) = label_opt else{
      return Err(CluEError::MissingPropertiesLabel(expression.line_number));
    };

    // Find the correct particle_config.
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

    if particle_config.properties == None{
      particle_config.properties = Some(ParticleProperties::new());
    }
    let Some(properties) = &mut particle_config.properties else{
      return Err(CluEError::MissingProperties( (*label).clone() ));
    };



    // Get relational symbol.
    match expression.relationship{
      Some(Token::Equals) => (),
      _ => return Err(CluEError::NoRelationalOperators(expression.line_number)),
    }
   


    let tokens = token_stream::extract_rhs(expression)?;

    let already_set = ||{
      CluEError::OptionAlreadySet(
          expression.line_number, expression.lhs[0].to_string()) };
    match expression.lhs[0]{
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      //Token::HyperfineCoupling =>,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      /*
      Token::HyperfineX => {
        is isotope specified?
        is axis already set?
        parse SecondaryParticleFilter
      },
      */
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
  //#[test]
  //fn test_parse_properties_line(){
  //}
}
