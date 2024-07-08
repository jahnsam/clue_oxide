use crate::space_3d::Vector3D;
use crate::physical_constants::{Element,Isotope};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `Particle` specifies possible particles.
/// 'Element' specifies a generalized chemical element.
/// 'Isotope' specifies a generalized isotope.
/// 'Coordinate' specifies the coordinates in 3D space of the particle.
/// 'active' indicates whether the particle is used for simulations.
/// 'serial' is the PDB serial ID.
/// 'residue' is the PDB residue.
/// 'residue_sequence_number' is the PDB residue sequence number.
#[derive(Debug,Clone)]
pub struct Particle{
  pub element: Element,
  pub isotope: Isotope,
  pub coordinates: Vector3D,
  pub active: bool,

  pub serial: Option<u32>,
  pub residue: Option<String>,
  pub residue_sequence_number: Option<u32>,
}

impl Particle{
  /// This function generates a 'Particle' of the specified element at
  /// the specified coordinates.
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






