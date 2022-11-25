use super::ParticleSpecifier;

pub struct Config{

  radius: f64,
  inner_radius: f64,
  max_number_of_cells: usize,
  pbc_style: PBCStyle,
  error_tolerance: f64,
  magnetic_field: Vec::<f64>,
  central_spin_coordinates: CentralSpinCoordinates,
  use_periodic_boundary_conditions: bool,
  isotope_abundances: IsotopeAbundaces,
  particles: Vec::<SpecifiedParticleProperties>
}

enum CentralSpinCoordinate{
  Atoms (Vec::<u32>),
  XYZ (f64, f64, f64),
}


pub enum PBCSyle{
  CRYST1,
  TIGHT,
}

struct SpecifiedParticleProperties{
  specifiers: ParticleSpecifier,
  properties: ParticleConfigProperties,
}

struct ParticleConfigProperties{
  isotopes: IsotopeInfo,
  tunnel_splitting:  f64,
}

impl ParticleConfigProperties{
  fn new() -> ParticleConfigProperties {
 
    let isotopes = most_common_isotopes();
    let abundaces = vec![vec![1.0]; isotopes.len()];
    let isotope_info = IsotopeInfo::new(isotopes, abundances);

    ParticleConfigProperties{
      isotopes: isotope_info,
      tunnel_splitting:  0.0,
    }
  }
}

struct IsotopeInfo{
particles: Vec::<Vec::<Particles>>,
  abundances: Vec::<Vec::<f64>>,
  force_no_pbc: bool,
  extracell_void_probability: f64,
}

impl IsotopeInfo{
  fn new(partilce_switch_list: Vec::<Particle>, probabilities: Vec::<f64>) 
    -> Result<IsotopeInfo, TODO_ERROR> {

      if partilce_switch_list,len() != probabilities.len() {
        return Err();
      }

      if sum(probabilities) != 1.0 {
        return Err();
      }

      if !all_positive(probabilities) {
        return Err();
      }
      
      Ok(IsotopeInfo{
          particles: particles,
          abundances: probabilities,
          force_no_pbc: false,
          extracell_void_probability: 0.0,
          })

  }
}
