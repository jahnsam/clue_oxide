use crate::space_3d::Vector3D;
use crate::physical_constants::{Element,Isotope};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct Particle{
  pub element: Element,
  pub isotope: Isotope,
  pub coordinates: Vector3D,
  pub active: bool,

  pub serial: Option<u32>,
  pub residue: Option<String>,
  pub residue_sequence_number: Option<u32>,
  //pub exchange_group: Option<ExchangeGroup>,
  //pub exchange_group_id: Option<usize>,  

  //pub isotope_distribution: Option<IsotopeDistribution>,

  //pub cell_id: usize,
}
impl Particle{
  pub fn new(element: Element, x: f64, y: f64, z: f64) -> Self{
    let isotope = Isotope::most_common_for(&element);
    let coordinates = Vector3D::from([x,y,z]);

    Particle{
      element,
      isotope,
      coordinates,
      active: true,
      serial: None,
      residue: None,
      residue_sequence_number: None,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
#[derive(Debug,Clone)]
pub struct IsotopeDistribution{
  pub isotope_abundances: Vec::<IsotopeAbundance>,
  pub force_no_pbc: bool,
  pub extracell_void_probability: f64,
}

impl Default for IsotopeDistribution{
  fn default() -> Self{
    IsotopeDistribution{
      isotope_abundances: Vec::<IsotopeAbundance>::new(),
      force_no_pbc: false,
      extracell_void_probability: 0.0,
    }
  }
}

#[derive(Debug,Clone)]
pub struct IsotopeAbundance{
  pub isotope: Isotope,
  pub abundance: f64,
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




