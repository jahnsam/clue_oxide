

pub struct Config{

  radius: f64,
  inner_radius: f64,
  magnetic_field: Vec::<f64>,
  central_spin_coordinates: CentralSpinCoordinates,
  use_periodic_boundary_conditions: bool,
  isotope_abundances: IsotopeAbundaces,
}

enum CentralSpinCoordinate{
  Atoms (Vec::<u32>),
  XYZ (f64, f64, f64),
}

struct IsotopeAbundaces{

}
