use super::particle_config::{ParticleConfig,find_particle_configs};
use super::particle_specifier::SpecifiedParticle;
use super::physical_constants::*;
use super::vector3::*;

#[derive(Debug,Clone)]
pub struct Config{
  pub radius: f64,
  pub inner_radius: f64,
  //max_number_of_cells: usize,
  //pbc_style: PBCStyle,
  //error_tolerance: f64,
  pub magnetic_field: Vector3,
  pub central_spin_coordinates: Option<CentralSpinCoordinates>,
  //use_periodic_boundary_conditions: bool,
  particles: Vec::<ParticleConfig>,
}

impl Default for Config{
  fn default() -> Self {
    Config{
      radius: f64::INFINITY,
      inner_radius: 0.0,
      magnetic_field: Vector3::from([0.0, 0.0, 1.2]),
      central_spin_coordinates: None,
      particles: Vec::<ParticleConfig>::new(),
    }
  }
}

impl Config{
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------

  pub fn get_max_spin_multiplicity_for_any_isotope(&self, element: Element)
    -> usize{

      let mut target = SpecifiedParticle::new();
      target.element = Some(element);

      let found_configs = find_particle_configs(&target, &self.particles);

      let mut max_spin_multiplicity = 1;

      for particle_config in found_configs{
        for iso_abu in &particle_config
          .config.isotopic_distribution.isotope_abundances{

          max_spin_multiplicity = max_spin_multiplicity.max(
              iso_abu.isotope.spin_multiplicity() );
        }
      }

      max_spin_multiplicity
    }
}
#[derive(Debug,Clone)]
pub enum CentralSpinCoordinates{
  Atoms (Vec::<u32>),
  XYZ (Vector3),
}


#[derive(Debug,Clone)]
pub enum PBCSyle{
  CRYST1,
  TIGHT,
}
