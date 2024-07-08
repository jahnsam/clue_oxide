
use crate::structure::particle_filter::{ParticleFilter,VectorSpecifier,
  SecondaryParticleFilter};
//use super::particle_specifier::*;
use crate::physical_constants::*;
use std::collections::HashMap;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ParticleConfig` stores user setting for how to treat particles/
/// `label` is an identifier.
/// `filter` selects the set of particles to modify.
/// `properties` specify the custom properties.
/// `cell_type` allows for treating the main cell and PBC copies differently.
#[derive(Debug,Clone,PartialEq)]
pub struct ParticleConfig{
  pub label: String,
  pub filter: Option<ParticleFilter>,
  pub properties: Option<ParticleProperties>,
  pub cell_type: CellType,
}

impl ParticleConfig{
  /// This function generates a default `ParticleConfig` with the input `label`.
  pub fn new(label: String) -> Self{
    ParticleConfig{
      label,
      filter: None,
      properties: None,
      cell_type: CellType::AllCells,
    }
  }
  //----------------------------------------------------------------------------
  /// This function looks at all spins that the `ParticleConfig` might apply to
  /// and returns the largest spin multiplicity, 2S+1.
  pub fn max_possible_spin_multiplicity(&self) -> Option<usize> {
    
    let properties = &self.properties.as_ref()?;

    if properties.isotopic_distribution.isotope_abundances.is_empty(){
      return None;
    }

    let mut max_spin_mult = 0;
    for isotope_abundance in properties.isotopic_distribution
      .isotope_abundances.iter()
    {
      max_spin_mult = usize::max(max_spin_mult,
          isotope_abundance.isotope.spin_multiplicity());
    }

    Some(max_spin_mult)
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This enum specifies different sets of PBC cells.
#[derive(Debug,Clone,PartialEq)]
pub enum CellType{
  AllCells,
  PrimaryCell,
  Extracells,
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ParticleProperties` specifies custom particle properties.
/// 'cosubstitute` selects the set of particles that should always be the same
/// isotope when the isotopic distribution is randomized.
/// `isotopic_distribution` specifies how elements are assigned an isotope.
/// `isotope_properties` defines some physical properties of the spin.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct ParticleProperties{
  pub cosubstitute: Option<SecondaryParticleFilter>,
  pub isotopic_distribution:  IsotopeDistribution,
  pub isotope_properties: HashMap::<String,IsotopeProperties>,
}

impl ParticleProperties{
  /// This function creates a new instance of `ParticleProperties`.
  pub fn new() -> Self{
    Default::default()
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `IsotopeProperties` defines some physical properties of a spin.
/// `active` specifies whether the spin should be included in simulations.
/// `electric_quadrupole_coupling` defines the coupling of the spin to the 
/// local electric field gradient.
/// `exchange_coupling` defines an effective coupling coming from
/// symmetry requirements of the wavefunction upon the exchange of 
/// identical particles.
/// `g_matrix` specifies the effective coupling between the spin and the applied
/// magnetic field.
/// `hyperfine_coupling` specifices the coupling to the detected electron.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct IsotopeProperties{
  pub active: Option<bool>,
  pub electric_quadrupole_coupling: Option<TensorSpecifier>,
  pub exchange_coupling: Option<f64>,
  pub g_matrix: Option<TensorSpecifier>,
  pub hyperfine_coupling: Option<TensorSpecifier>,
}

impl IsotopeProperties{
  /// This function creates a new instance of `IsotopeProperties`.
  pub fn new() -> Self{ IsotopeProperties::default() }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `IsotopeDistribution` determine how isotopes are assigned to elements.
/// `isotope_abundances` lists the possible isotope and propabilities.
/// 'void_probability' is the probability that the particle will be included
/// at all.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct IsotopeDistribution{
  pub isotope_abundances: Vec::<IsotopeAbundance>,
  pub void_probability: Option<f64>,
}
//------------------------------------------------------------------------------
/// `IsotopeAbundances` lists the possible isotope and abundances.
#[derive(Debug,Clone,PartialEq)]
pub struct IsotopeAbundance{
  pub isotope: Isotope,
  pub abundance: f64,         
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `TensorSpecifier` specifies a symmetric 3-by-3 coupling matrix.
/// `values` contains the three eigenvalues.
/// `e_axis` for `e` in {'x','y','z'} specify eigenvectors.
/// Note that only two axes should be specified. 
#[derive(Debug,Clone,PartialEq,Default)]
pub struct TensorSpecifier{
  pub values: Option<[f64; 3]>,
  pub x_axis: Option<VectorSpecifier>,
  pub y_axis: Option<VectorSpecifier>,
  pub z_axis: Option<VectorSpecifier>,
}

impl TensorSpecifier{
  /// This function creates a new instance of `TensorSpecifier`.
  pub fn new() -> Self{
    TensorSpecifier::default()
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



