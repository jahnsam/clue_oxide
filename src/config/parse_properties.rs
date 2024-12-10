use crate::clue_errors::CluEError;
use crate::config::{
  Config,
  particle_config::{
    IsotopeAbundance,IsotopeProperties,
    ParticleConfig, ParticleProperties,
    EigSpecifier, TensorSpecifier},
  to_f64_token,
  Token,
  token_expressions::*,
  token_stream::split_on_token,
};
use crate::structure::particle_filter::SecondaryParticleFilter;
use crate::physical_constants::Isotope;
use crate::space_3d::SymmetricTensor3D;


impl Config{
  ///  This function parses lines found under "#[spin_properties(...)]"
  /// or "#[structure_properties(...)]".
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
      None => parse_structure_properties(properties,expression,label),
    } 
  }
}
//----------------------------------------------------------------------------
///  This function parses lines found under "#[spin_properties(...)]".
fn parse_isotope_properties(properties: &mut ParticleProperties, 
    expression: &TokenExpression, label: &str, isotope: &Isotope) 
-> Result<(),CluEError>
{


  let key = isotope.to_string();
  
  
  match properties.isotope_properties.get(&key){
    Some(_) => (),
    None => {
      properties.isotope_properties.insert(key.clone(),
          IsotopeProperties::new());
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
    Token::Active => {
      if isotope_properties.active.is_some(){
        return Err(already_set());
      }

      set_to_some_bool(&mut isotope_properties.active,expression)?;
    },
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Token::GMatrix | Token::GX | Token::GY | Token::GZ 
      => set_symmetric_tensor_3d(&mut isotope_properties.g_matrix, 
          &expression.lhs[0], expression, label)?,
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Token::HyperfineCoupling | Token::HyperfineX | 
        Token::HyperfineY | Token::HyperfineZ 
      => set_symmetric_tensor_3d(&mut isotope_properties.hyperfine_coupling, 
          &expression.lhs[0], expression, label)?,
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Token::ElectricQuadrupoleCoupling | Token::ElectricQuadrupoleX
        | Token::ElectricQuadrupoleY | Token::ElectricQuadrupoleZ 
      => set_symmetric_tensor_3d(
             &mut isotope_properties.electric_quadrupole_coupling, 
             &expression.lhs[0], expression, label)?,
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
//------------------------------------------------------------------------------
///  This function parses lines found under "#[structure_properties(...)]".
fn parse_structure_properties(properties: &mut ParticleProperties, 
    expression: &TokenExpression, _label: &str) 
-> Result<(),CluEError>
{

  // Get relational symbol.
  match expression.relationship{
    Some(Token::Equals) => (),
    _ => return Err(CluEError::NoRelationalOperators(expression.line_number)),
  }
 

  let already_set = ||{
    CluEError::OptionAlreadySet(
        expression.line_number, expression.lhs[0].to_string()) };
  match expression.lhs[0]{
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Token::Cosubstitute => {
      if properties.cosubstitute.is_some(){
        return Err(already_set());
      }

      let Some(rhs) = expression.rhs.as_ref() else{
        return Err(CluEError::NoRHS(expression.line_number));
      }; 

      let sec_fltr = SecondaryParticleFilter::from(&rhs[0].to_string())?;

      if sec_fltr != SecondaryParticleFilter::SameMolecule{
        return Err(CluEError::InvalidSecondaryFilter(expression.line_number,
              sec_fltr.to_string()));
      }

      properties.cosubstitute = Some(sec_fltr);

    },
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Token::IsotopeAbundances => {
      if !properties.isotopic_distribution.isotope_abundances.is_empty(){
        return Err(already_set());
      }
      properties.isotopic_distribution.isotope_abundances = 
        parse_isotope_abundances(expression)?;

    },
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Token::VoidProbability => {
      if properties.isotopic_distribution.void_probability.is_some(){
        return Err(already_set());
      }
      set_to_some_f64(&mut properties.isotopic_distribution
          .void_probability, expression)?;
    },
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    _ => return Err(CluEError::InvalidToken(expression.line_number,
          expression.lhs[0].to_string())),
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  }
  
  Ok(())
}      
//------------------------------------------------------------------------------
///  This function parses the lines found under "#[structure_properties(...)]"
/// that determine isotope identities and abundances.
fn parse_isotope_abundances(expression: &TokenExpression)
  -> Result< Vec::<IsotopeAbundance>, CluEError>
{

  let Some(rhs) = &expression.rhs else{
    return Err(CluEError::IncorrectFormattingIsotopeAbundances(
          expression.line_number));
  };

  // Check formatting.
  let n_tok = rhs.len();
  if n_tok < 5
    || rhs[0] != Token::CurlyBracketOpen 
    || rhs[n_tok-1] != Token::CurlyBracketClose{

    return Err(CluEError::IncorrectFormattingIsotopeAbundances(
          expression.line_number));
  }

  // Extract key and values.
  let key_values = split_on_token(rhs[1..n_tok-1].to_vec(), Token::Comma);

  // Initialize output.
  let mut isotope_abundances =
    Vec::<IsotopeAbundance>::with_capacity(key_values.len());

  // Initialize nomalization.
  let mut normalization = 0.0;

  // Loop over entries.
  for field in key_values.iter(){

    // Check formatting.
    if field.len() != 3 || field[1] != Token::Colon{
      return Err(CluEError::IncorrectFormattingIsotopeAbundances(
            expression.line_number));
    }

    // Extract isotope.
    let isotope = Isotope::from(&field[0].to_string())?;

    // Extract abundance.
    let value_token = to_f64_token(vec![field[2].clone()], 
        expression.line_number)?;
    let Token::Float(abundance) = value_token else {
      return Err(CluEError::IncorrectFormattingIsotopeAbundances(
            expression.line_number));
    };

    if abundance < 0.0{
      return Err(CluEError::IsotopeAbundancesMustBeNonnegative(
            expression.line_number));
    }

    // Update.

    normalization += abundance;

    isotope_abundances.push(IsotopeAbundance{isotope,abundance});
  }

  if normalization <= 0.0{
      return Err(CluEError::IsotopeAbundancesCannotBeNormalized(
            expression.line_number));
  }

  // Ensure sum(abundances) = 1.
  for isotope_abundance in isotope_abundances.iter_mut(){
    isotope_abundance.abundance /= normalization;
  }

  Ok(isotope_abundances)
}
//------------------------------------------------------------------------------
// TODO: update errors
pub fn set_symmetric_tensor_3d(matrix_opt: &mut Option<TensorSpecifier>, 
    token: &Token, expression: &TokenExpression, label: &str)
  -> Result<(),CluEError>
{

  if matrix_opt.is_none() {
    *matrix_opt = Some(TensorSpecifier::new());
  }

  let already_set = ||{
    CluEError::OptionAlreadySet(
        expression.line_number, expression.lhs[0].to_string()) };

  match token{
    // TODO: There should be a beter method to catch these cases.
    Token::DetectedSpinGMatrix | Token::GMatrix | 
      Token::HyperfineCoupling | Token::ElectricQuadrupoleCoupling 
    => {
    
      // Extract the user input values.
      let mut vals = Vec::<f64>::new();    
      set_to_vec_f64(&mut vals,expression)?;
      

      // The number of values determines how to set the matrix.
      match vals.len(){
        3 => {

          let new_values = Some([vals[0],vals[1],vals[2]]);
      
          let mut matrix = match matrix_opt{
            Some(TensorSpecifier::Unspecified) => EigSpecifier::new(),

            Some(TensorSpecifier::Eig(eig_specifier)) => eig_specifier.clone(),

            _ => return Err(already_set()),
          };
        
          if matrix.values.is_some(){ 
            return Err(already_set()); 
          };

          matrix.values = new_values;
        
          *matrix_opt = Some(TensorSpecifier::Eig(matrix));
          
        },

        1 | 6 | 9 =>{   
          if *matrix_opt != Some(TensorSpecifier::Unspecified){
            return Err(CluEError::CannotInferEigenvalues(
                expression.line_number));
          } 
          let values = if vals.len() == 1{
            [vals[0],     0.0,     0.0,
                      vals[0],     0.0,
                               vals[0]]

          }else if vals.len() == 6{
            [ vals[0], vals[1], vals[2],
                       vals[3], vals[4],
                                vals[5]] 
          }else{
            [ vals[0], vals[1], vals[2],
                       vals[4], vals[5],
                                vals[8]] 
          };

          *matrix_opt = Some(TensorSpecifier::SymmetricTensor3D(
              SymmetricTensor3D::from(values)));
        },

        _ => return Err(CluEError::CannotInferEigenvalues(
                expression.line_number)),
      }
    },
    // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    // TODO: There should be a beter method to catch these cases.
    Token::DetectedSpinGX | Token::DetectedSpinGY | Token::DetectedSpinGZ |
        Token::GX | Token::GY | Token::GZ  |
        Token::HyperfineX | Token::HyperfineY | Token::HyperfineZ |
        Token::ElectricQuadrupoleX | Token::ElectricQuadrupoleY |
        Token::ElectricQuadrupoleZ
    => {

      // Check that axes can be assigned.
      let mut matrix = match matrix_opt{
        Some(TensorSpecifier::Unspecified) => EigSpecifier::new(),

        Some(TensorSpecifier::Eig(eig_specifier)) => eig_specifier.clone(),

        _ => return Err(already_set()),
      };


      // The the appropriate axis.
      match expression.lhs[0]{
        // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        Token::DetectedSpinGX | Token::GX | 
          Token::HyperfineX | Token::ElectricQuadrupoleX
        => {
          set_to_some_vector_specifier(&mut matrix.x_axis, expression,
              label)?;
        },
        // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        Token::DetectedSpinGY | Token::GY | 
          Token::HyperfineY | Token::ElectricQuadrupoleY
        => {
          set_to_some_vector_specifier(&mut matrix.y_axis, expression,
              label)?;
        },
        // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        Token::DetectedSpinGZ | Token::GZ | 
          Token::HyperfineZ | Token::ElectricQuadrupoleZ
        => {

          set_to_some_vector_specifier(&mut matrix.z_axis, expression,
              label)?;

        },
        // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        _ => return Err(CluEError::InvalidToken(expression.line_number,
          expression.lhs[0].to_string())),
      }

      *matrix_opt = Some(TensorSpecifier::Eig(matrix));
    },
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    _ => return Err(CluEError::InvalidToken(expression.line_number,
          expression.lhs[0].to_string())),
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  }
  Ok(())
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        active = false;
        tunnel_splitting = 60e3;\
        hyperfine_coupling = 12.3e2;\
        electric_quadrupole_coupling = [1.2, 3.4, 5.6];\
        electric_quadrupole_x = vector([-1,0,1]);\
        electric_quadrupole_y = diff(bonded(test_label1),\
          bonded(test_label2));\
        electric_quadrupole_z = diff(same_molecule(test_label1),\
          particle);\
        g_matrix = [1, 2, 3,
                    4, 5, 6,
                    7, 8, 9];\
        ").unwrap();

    let isotope = Isotope::Hydrogen2;

    let mut properties = ParticleProperties::new();
    for expression in expressions.iter() {
      parse_isotope_properties(&mut properties, expression, "test_label0", 
          &isotope).unwrap(); 
    }

    let d_properties = &properties.isotope_properties["2H"];

    let hyperfine_coupling = match d_properties.hyperfine_coupling
        .as_ref().unwrap(){
          TensorSpecifier::SymmetricTensor3D(tensor) => tensor.clone(),
          _ => panic!("Expected TensorSpecifier::SymmetricTensor3D(tensor)."),
    };

    let quadrupole_coupling 
      = match d_properties.electric_quadrupole_coupling.as_ref().unwrap(){
          TensorSpecifier::Eig(eig_specifier) => eig_specifier.clone(),
          _ => panic!("Expected TensorSpecifier::Eig(eig_specifier)."),
    };

    let g_matrix = match d_properties.g_matrix.as_ref().unwrap(){
          TensorSpecifier::SymmetricTensor3D(tensor) => tensor.clone(),
          _ => panic!("Expected TensorSpecifier::SymmetricTensor3D(tensor)."),
    };
    

    assert_eq!(d_properties.active, Some(false) );
    
    assert_eq!(d_properties.exchange_coupling, Some(-40e3) );

    assert_eq!(hyperfine_coupling, SymmetricTensor3D::from([12.3e2, 0.0, 0.0,
                                                                 12.3e2, 0.0,
                                                                      12.3e2]));
    
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

    assert_eq!(g_matrix , SymmetricTensor3D::from([ 1.0, 2.0, 3.0,
                                                         5.0, 6.0,
                                                              9.0])
    );
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_structure_properties(){
    let expressions = get_tokens_from_line("\
      isotope_abundances = {1H: 8, 2H: 2};
      void_probability = 0.5;
    ").unwrap();

    let mut properties = ParticleProperties::new();
    for expression in expressions.iter() {
      parse_structure_properties(&mut properties, expression, "label").unwrap(); 
    }

    
    assert_eq!(properties.isotopic_distribution.isotope_abundances, 
      vec![
      IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: 0.8},
      IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 0.2},
      ]);

    assert_eq!(properties.isotopic_distribution.void_probability,
        Some(0.5));
      
      
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
