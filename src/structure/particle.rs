use crate::space_3d::Vector3D;
use crate::physical_constants::{Element,Isotope};

/* Pseudocode Overview   
  
   let config = read_input(file);

   for pdb in config.pdbs{
     let structure = parse_pdb()
                     .set_spins()
                     .apply_pbc();
    
    let spin_op = get_spin_operators(max_cluster_size, max_spin_mult);                 

    let tensors = get_tensors(&structure,&config);

    let clusters = find_clusters();

    clusters.remove_partial_methyls();

    let clusters.calculate_cluster_signals();

    let clusters.calculate_auxiliary_signals();

    let structure_signal = do_cce(clusters);

    return structure_signal;

   }
*/
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct Particle{
  pub element: Element,
  pub isotope: Isotope,
  pub coordinates: Vector3D,

  pub serial: Option<u32>,
  pub residue: Option<String>,
  pub residue_sequence_number: Option<u32>,
  //pub exchange_group: Option<ExchangeGroup>,
  //pub exchange_group_id: Option<usize>,  

  //pub isotope_distribution: Option<IsotopeDistribution>,

  //pub cell_id: usize,
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




