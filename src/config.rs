use strum::IntoEnumIterator;

use super::particle_config::*;
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
  pub particles: Vec::<ParticleConfig>,
}

impl Default for Config{
  fn default() -> Self {
    let mut particles = Vec::<ParticleConfig>::new();
    for element in Element::iter(){

      let mut particle_config = ParticleConfig::new();
      particle_config.filter.elements.push(element);

      particle_config.config.isotopic_distribution.isotope_abundances.push(
          IsotopeAbundance{
            isotope: Isotope::most_common_for(&element),
            abundance: 1.0,
          });


      particles.push(particle_config);
    }

    Config{
      radius: f64::INFINITY,
      inner_radius: 0.0,
      magnetic_field: Vector3::from([0.0, 0.0, 1.2]),
      central_spin_coordinates: None,
      particles: particles,
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
