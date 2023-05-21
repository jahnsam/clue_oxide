use crate::clue_errors::CluEError;
use crate::Config;
use crate::config::Token;
//use crate::config::token_stream;
use crate::config::particle_config::ParticleConfig;
use crate::config::token_expressions::*;
use crate::config::particle_config::{IsotopeProperties,ParticleProperties,
  TensorSpecifier};
use crate::physical_constants::Isotope;

impl Config{
  pub fn parse_properties_line(&mut self, expression: &TokenExpression,
      label_opt: &Option<String>,isotope_opt: Option<Isotope>)
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

    let particle_config: &mut ParticleConfig = match found_particle_config{
      Some(idx) => &mut self.particles[idx],
      None => {
        self.particles.push(ParticleConfig::new((*label).clone()));
        let idx = self.particles.len() - 1; 
        &mut self.particles[idx]}
    };

    if particle_config.properties.is_none(){
      particle_config.properties = Some(ParticleProperties::new());
    }
    let Some(properties) = &mut particle_config.properties else{
      return Err(CluEError::MissingProperties( (*label).clone() ));
    };


    match isotope_opt{
      Some(isotope) => parse_isotope_properties(properties,expression,
          label, &isotope),
      None => Ok(()),//parse_element_properties()
    } 
  }
}
  //----------------------------------------------------------------------------
  fn parse_isotope_properties(properties: &mut ParticleProperties, 
      expression: &TokenExpression, label: &str, isotope: &Isotope) 
  -> Result<(),CluEError>
  {


    let key = isotope.to_string();
    
    
    match properties.isotope_properties.get(&key){
      Some(_) => (),
      None => {
        properties.isotope_properties.insert(key.clone(),IsotopeProperties::new());
      },
    }
    let mut isotope_properties = properties.isotope_properties[&key].clone();

    // Get relational symbol.
    match expression.relationship{
      Some(Token::Equals) => (),
      _ => return Err(CluEError::NoRelationalOperators(expression.line_number)),
    }
   
    //let tokens = token_stream::extract_rhs(expression)?;

    let already_set = ||{
      CluEError::OptionAlreadySet(
          expression.line_number, expression.lhs[0].to_string()) };
    match expression.lhs[0]{
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::HyperfineCoupling | Token::HyperfineX | 
        Token::HyperfineY | Token::HyperfineZ =>{

        let mut hf_values = Vec::<f64>::new();
      
        set_to_vec_f64(&mut hf_values,expression)?;
      

        if isotope_properties.hyperfine_coupling.is_none() {
          isotope_properties.hyperfine_coupling = Some(TensorSpecifier::new());
        }
        let Some(hyperfine) = &mut isotope_properties.hyperfine_coupling else{
          return Err(CluEError::NoHyperfineSpecifier(label.to_string(),
                isotope.to_string() ));
        }; 

        match expression.lhs[0]{
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::HyperfineCoupling => {
            if hyperfine.values.is_some(){
              return Err(already_set());
            }

            if hf_values.len() == 1{
              hyperfine.values = Some([hf_values[0],hf_values[0],hf_values[0]]);
            }else if hf_values.len() == 3{
              hyperfine.values = Some([hf_values[0],hf_values[1],hf_values[2]]);
            }else{
              return Err(CluEError::CannotInferEigenvalues(
                    expression.line_number));
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::HyperfineX => {
            set_to_some_vector_specifier(&mut hyperfine.x_axis, expression,
                label)?;
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::HyperfineY => {
            set_to_some_vector_specifier(&mut hyperfine.y_axis, expression,
                label)?;
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::HyperfineZ => {
            set_to_some_vector_specifier(&mut hyperfine.z_axis, expression,
                label)?;
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ElectricQuadrupoleCoupling | Token::ElectricQuadrupoleX
        | Token::ElectricQuadrupoleY | Token::ElectricQuadrupoleZ =>{

        let mut qc_values = Vec::<f64>::new();
      
        set_to_vec_f64(&mut qc_values, expression)?;
      

        if isotope_properties.electric_quadrupole_coupling.is_none() {
          isotope_properties.electric_quadrupole_coupling =
            Some(TensorSpecifier::new());
        }
        let Some(quadrupole) 
          = &mut isotope_properties.electric_quadrupole_coupling else{
          return Err(CluEError::NoQuadrupoleSpecifier(label.to_string(),
                isotope.to_string() ));
        }; 

        match expression.lhs[0]{
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ElectricQuadrupoleCoupling => {
            if quadrupole.values.is_some(){
              return Err(already_set());
            }
 
            if qc_values.len() == 1{
              quadrupole.values 
                = Some([qc_values[0],qc_values[0],qc_values[0]]);
            }else if qc_values.len() == 3{
              quadrupole.values 
                = Some([qc_values[0],qc_values[1],qc_values[2]]);
            }else{
              return Err(CluEError::CannotInferEigenvalues(
                    expression.line_number));
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ElectricQuadrupoleX => {
            set_to_some_vector_specifier(&mut quadrupole.x_axis, expression,
                label)?;
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ElectricQuadrupoleY => {
            set_to_some_vector_specifier(&mut quadrupole.y_axis, expression,
                label)?;
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ElectricQuadrupoleZ => {
            set_to_some_vector_specifier(&mut quadrupole.z_axis, expression,
                label)?;
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::TunnelSplitting => {
        let mut nut_opt = isotope_properties.exchange_coupling;
        set_to_some_f64(&mut nut_opt,expression)?;
        
        if let Some(nut) = nut_opt{
          isotope_properties.exchange_coupling = Some(-2.0*nut/3.0);
        }else{
          return Err(CluEError::CannotSetExchangeCoupling(
                expression.line_number));
        }
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
    
    properties.isotope_properties.insert(key, isotope_properties);
    Ok(())

  }



#[cfg(test)]
mod tests{
  use super::*;
  use crate::config::get_tokens_from_line;
  use crate::structure::particle_filter::VectorSpecifier;
  use crate::space_3d::Vector3D;
  use crate::structure::particle_filter::SecondaryParticleFilter;

  #[test]
  fn test_parse_isotope_properties(){

    let expressions = get_tokens_from_line("\
        tunnel_splitting = 60e3;\
        hyperfine_coupling = 12.3e2;\
        hyperfine_x = vector([-1,0,1]);\
        hyperfine_y = diff(particle, bonded(test_label1));\
        hyperfine_z = diff(particle, same_molecule(test_label1));\
        electric_quadrupole_coupling = [1.2, 3.4, 5.6];\
        electric_quadrupole_x = vector([-1,0,1]);\
        electric_quadrupole_y = diff(bonded(test_label1),\
          bonded(test_label2));\
        electric_quadrupole_z = diff(same_molecule(test_label1),\
          particle);\
        ").unwrap();

    let isotope = Isotope::Hydrogen2;

    let mut properties = ParticleProperties::new();
    for expression in expressions.iter() {
      parse_isotope_properties(&mut properties, expression, "test_label0", 
          &isotope).unwrap(); 
    }

    let d_properties = &properties.isotope_properties["2H"];
    let hyperfine_coupling = d_properties.hyperfine_coupling.as_ref().unwrap();
    let quadrupole_coupling = d_properties.electric_quadrupole_coupling.as_ref().unwrap();
    
    assert_eq!(d_properties.exchange_coupling, Some(-40e3) );

    assert_eq!(hyperfine_coupling.values , Some([12.3e2, 12.3e2, 12.3e2]));
    assert_eq!(hyperfine_coupling.x_axis , Some(
          VectorSpecifier::Vector(Vector3D::from([-1.0, 0.0, 1.0]))));
    assert_eq!(hyperfine_coupling.y_axis , Some(
          VectorSpecifier::Diff(
            SecondaryParticleFilter::Particle, "test_label0".to_string(),
            SecondaryParticleFilter::Bonded, "test_label1".to_string(),
            ) ) );
    assert_eq!(hyperfine_coupling.z_axis , Some(
          VectorSpecifier::Diff(
            SecondaryParticleFilter::Particle, "test_label0".to_string(),
            SecondaryParticleFilter::SameMolecule, "test_label1".to_string(),
            ) ) );
    
    assert_eq!(quadrupole_coupling.values , Some([1.2, 3.4, 5.6]));
    assert_eq!(quadrupole_coupling.x_axis , Some(
          VectorSpecifier::Vector(Vector3D::from([-1.0, 0.0, 1.0]))));
    assert_eq!(quadrupole_coupling.y_axis , Some(
          VectorSpecifier::Diff(
            SecondaryParticleFilter::Bonded, "test_label1".to_string(),
            SecondaryParticleFilter::Bonded, "test_label2".to_string(),
            ) ) );
    assert_eq!(quadrupole_coupling.z_axis , Some(
          VectorSpecifier::Diff(
            SecondaryParticleFilter::SameMolecule, "test_label1".to_string(),
            SecondaryParticleFilter::Particle, "test_label0".to_string(),
            ) ) );
    }
}
