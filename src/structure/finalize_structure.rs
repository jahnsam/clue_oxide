use crate::clue_errors::CluEError;
use crate::config::Config;
use crate::space_3d::Vector3D;
use crate::structure::Structure;
use crate::math;  

use rand::seq::index::sample;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Binomial, Distribution, Uniform};
use std::collections::HashMap;
use crate::config::particle_config::IsotopeDistribution;

impl Structure{
  /// This function takes a structure from build_extended_structure(), 
  /// trims the system and sets isotope identities.
  pub fn finalize_structure(&mut self,
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    //self.trim_system(config)?;

    // Set isotpic identities after adding voidable particles since the
    // non-void particles can potentially have multiple isotopic options. 
    self.set_isotopologue(rng, config)?;


    Ok(())
  }
  //----------------------------------------------------------------------------
  /*
  fn trim_system(&mut self, configL &Config) -> Result<(),CluEError>{
    let Some(load_geometry) = &config.load_geometry else {
      return Err(CluEError::NoLoadGeometry);
    };
    match load_geometry{
      LoadGeometry::Cube => (),
      LoadGeometry::Sphere =>{  
        let Some(radius) = config.radius else {
          return Err(CluEError::NoRadius);
        };

        for particle in self.bath_particles.iter_mut(){

          if (particle.coordinates - r_electron).norm() > radius{
            particle.active = false;
          }
        }
      }

    Ok(())
  }
  */
  //----------------------------------------------------------------------------
  fn set_isotopologue(&mut self, rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {
    // Initialize rng range.
    let range = Uniform::new(0.0f64, 1.0);

    // Check if there are any configurations.
    let particle_configs = &config.particles;

    for unit_cell_id in 0..self.cell_offsets.len(){

      // Loop over bath particle.
      for cosubstitution_group in self.cosubstitution_groups.iter(){

        // Generate a random number.
        let random_number = range.sample(rng);

        for particle_idx0 in cosubstitution_group.iter(){
          if !self.bath_particles[*particle_idx0].active {continue;}
        
          let Some(particle_idx) 
            = self.cell_indices[*particle_idx0][unit_cell_id] else{continue;};

          let particle = &mut self.bath_particles[particle_idx];
      
          // Check if this particle has a custom config.
          let Some(config_id) = self.particle_config_ids[*particle_idx0] else {
            continue;
          };

          // Check if this particle has any custom properties.
          let Some(properties) = &particle_configs[config_id].properties else{
            continue;
          };


          let mut cdf = 0.0;
     
          // Select a isotope for this particle.
          let isotopic_distribution: &IsotopeDistribution;
          if unit_cell_id == 0{
            isotopic_distribution = &properties.isotopic_distribution;
          }else{
            if let Some(iso) = &properties.extracell_isotopic_distribution{
              isotopic_distribution = iso;
            }else{
              isotopic_distribution = &properties.isotopic_distribution;
            } 
          }

          for iso in isotopic_distribution.isotope_abundances.iter(){
            cdf += iso.abundance;
            if cdf >= random_number{
              particle.isotope = iso.isotope;
              particle.active = particle.isotope.spin_multiplicity() > 1;
              break;
            } 
          }

        }
      }
    }

    Ok(())

  }
  //----------------------------------------------------------------------------

}



#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb;
  use crate::structure::particle_filter::*;
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

  use crate::config::particle_config::{ParticleConfig,
    ParticleProperties,IsotopeDistribution,IsotopeAbundance};
  use crate::physical_constants::{Element,Isotope};

  
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_extended_structure_lone_tempo(){
    let filename = "./assets/TEMPO.pdb";
    let mut structures = pdb::parse_pdb(&filename).unwrap();

    let mut particle_configs =vec![ParticleConfig::new("hydrogen".to_string())];

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    particle_configs[0].filter = Some(filter);
    
    let mut properties = ParticleProperties::new();
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: 0.5});
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 0.5});

    particle_configs[0].properties = Some(properties);
    
    let mut config = Config::new();
    config.load_geometry = Some(LoadGeometry::Cube);
    config.particles = particle_configs;
    config.radius = Some(73.5676e-10);


    structures[0].build_primary_structure(&config).unwrap();
    let num_particles = structures[0].bath_particles.len();

    let mut rng =  ChaCha20Rng::from_entropy();

    let mut structure = structures[0].clone();
    structure.build_extended_structure(&mut rng, &config);

    // Test: set_cell_shifts().
    assert_eq!(structure.cell_offsets.len(),125);

    // Test: extend_structure().
    let expected_num = num_particles + 124*19;
    assert_eq!(structure.bath_particles.len(),expected_num);

    // Test: set_isotopologue().
    let mut n_h1 = 0;
    let mut n_h2 = 0;
    for particle in structure.bath_particles.iter(){
      if (*particle).isotope == Isotope::Hydrogen1{
        n_h1 += 1;
      }else if (*particle).isotope == Isotope::Hydrogen2{
        n_h2 += 1; 
      }
    }
    assert_eq!(n_h1+n_h2,125*18);
    assert!(n_h1 >= 1000);
    assert!(n_h1 <= 1200);
    assert!(n_h2 >= 1000);
    assert!(n_h2 <= 1200);
    
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_add_voidable_particles(){
    // n    : Molecules    : Hydrons 
    // 1    : TEMPO        :    18
    // 1500 : glycerols    : 12000
    // 7469 : waters       : 14938
    //
    // total hydrons =  26956.
    //
    let n_wat = 7469;
    let n_gly = 1500;
    let n_ex = 2*n_wat + 3*n_gly; 
    let n_nx = 5*n_gly + 18;
    let filename = "./assets/TEMPO_wat_gly_70A.pdb";
    let mut structures = pdb::parse_pdb(&filename).unwrap();


    let mut particle_configs =vec![
      ParticleConfig::new("exchangeable_hydrogens".to_string()),
      ParticleConfig::new("nonexchangeable_hydrogens".to_string()),
    ];

    // filter 0 
    let mut filter_ex = ParticleFilter::new();
    filter_ex.elements = vec![Element::Hydrogen];
    filter_ex.bonded_elements = vec![Element::Oxygen];
    particle_configs[0].filter = Some(filter_ex.clone());

    // filter 1
    let mut filter_nx = ParticleFilter::new();
    filter_nx.elements = vec![Element::Hydrogen];
    filter_nx.not_bonded_elements = vec![Element::Oxygen];
    particle_configs[1].filter = Some(filter_nx.clone());
    
    let mut properties = ParticleProperties::new();
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: 0.5});
    
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 0.5});
    
    properties.isotopic_distribution.extracell_void_probability = Some(0.5);

    let mut extracell_isotopic_distribution = IsotopeDistribution::default();
    extracell_isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: 1.0});
    properties.extracell_isotopic_distribution 
      = Some(extracell_isotopic_distribution);

    particle_configs[0].properties = Some(properties.clone());

    properties.cosubstitute = Some(SecondaryParticleFilter::SameMolecule);
    particle_configs[1].properties = Some(properties);

    let mut config = Config::new();
    config.load_geometry = Some(LoadGeometry::Cube);
    config.particles = particle_configs;
    config.radius = Some(73.5676e-10);
    let n_uc = 125;


    structures[0].build_primary_structure(&config).unwrap();


    

    let mut rng =  ChaCha20Rng::from_entropy();

    let structure = &mut structures[0];
    structure.build_extended_structure(&mut rng, &config).unwrap();

    assert_eq!(structure.cell_offsets.len(),n_uc);
    //assert_eq!(structure.bath_particles.len(),n_uc*num_particles);

    // Check that outside the primary cell deuterons are removed.
    for (idx, particle) in structure.bath_particles.iter().enumerate(){
      if structure.cell_id(idx).unwrap() > 0{
        assert!(
            (*particle).isotope==Isotope::Hydrogen1 
            || (*particle).isotope==Isotope::Nitrogen14 
            );

      }else{
        let idx0 = structure.primary_cell_indices[idx];
        assert_eq!(idx,idx0);
      }

      if (*particle).isotope == Isotope::Hydrogen2{
        let idx0 = structure.primary_cell_indices[idx];
        assert_eq!(idx,idx0);
      }
    }


    // Check for for the correct exchanges.
    let [n_h, n_h1, n_h2,n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2,_n_mol] =
      get_conditional_hydron_stats(&structure,&filter_ex);

    let approx_eq = |m: usize,n: usize| -> bool{
      (1.0 - (m as f64)/(n as f64) ).abs() <1e3 
    };
    assert!(n_h1_h1>0);
    assert!(n_h2_h2>0);
    assert!(n_h2_h1>0);
    assert!(n_h1_h2>0);
    assert!( approx_eq(n_h1,n_uc*n_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h1) );
    assert!( approx_eq(n_h1_h1,n_uc*n_h1_h2) );
    assert!( approx_eq(n_h, n_ex*(1 +(n_uc-1)/2)  ));

    let [n_h, n_h1, n_h2, n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2,_n_mol] =
      get_conditional_hydron_stats(&structure,&filter_nx);

    assert!(n_h1_h1>0);
    assert!(n_h2_h2>0);
    assert_eq!(n_h2_h1,0);
    assert_eq!(n_h1_h2,0);
    assert!( approx_eq(n_h1,n_uc*n_h2) );
    assert!( approx_eq(n_h1_h1,n_uc*n_h2_h2) );
    assert!( approx_eq(n_h, n_nx*(1 +(n_uc-1)/2)  ));

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_extended_structure_tempo_wat_gly_7nm(){
    // n    : Molecules    : Hydrons 
    // 1    : TEMPO        :    18
    // 1500 : glycerols    : 12000
    // 7469 : waters       : 14938
    //
    // total hydrons =  26956.
    //
    let n_wat = 7469;
    let n_gly = 1500;
    let n_ex = 2*n_wat + 3*n_gly; 
    let n_nx = 5*n_gly + 18;

    let filename = "./assets/TEMPO_wat_gly_70A.pdb";
    let mut structures = pdb::parse_pdb(&filename).unwrap();


    let mut particle_configs =vec![
      ParticleConfig::new("exchangeable_hydrogens".to_string()),
      ParticleConfig::new("nonexchangeable_hydrogens".to_string()),
    ];

    let mut filter_ex = ParticleFilter::new();
    filter_ex.elements = vec![Element::Hydrogen];
    filter_ex.bonded_elements = vec![Element::Oxygen];
    particle_configs[0].filter = Some(filter_ex.clone());

    let mut filter_nx = ParticleFilter::new();
    filter_nx.elements = vec![Element::Hydrogen];
    filter_nx.not_bonded_elements = vec![Element::Oxygen];
    particle_configs[1].filter = Some(filter_nx.clone());
    
    let mut properties = ParticleProperties::new();
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: 0.5});
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 0.5});


    particle_configs[0].properties = Some(properties.clone());

    properties.cosubstitute = Some(SecondaryParticleFilter::SameMolecule);
    particle_configs[1].properties = Some(properties);
    
    let mut config = Config::new();
    config.load_geometry = Some(LoadGeometry::Cube);
    config.particles = particle_configs;
    config.radius = Some(73.5676e-10);
    let n_uc = 125;


    structures[0].build_primary_structure(&config).unwrap();
    let num_particles = structures[0].bath_particles.len();


    

    let mut rng =  ChaCha20Rng::from_entropy();

    let mut structure = structures[0].clone();
    structure.build_extended_structure(&mut rng, &config);

    assert_eq!(structure.cell_offsets.len(),n_uc);



    // Check for for the correct exchanges.
    let [n_h, n_h1, n_h2,n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2,_n_mol] =
      get_conditional_hydron_stats(&structure,&filter_ex);

    assert_eq!(n_h,n_uc*n_ex);
    let approx_eq = |m: usize,n: usize| -> bool{
      (1.0 - (m as f64)/(n as f64) ).abs() <1e3 
    };
    assert!(n_h1_h1>0);
    assert!(n_h2_h2>0);
    assert!(n_h2_h1>0);
    assert!(n_h1_h2>0);
    assert!( approx_eq(n_h1,n_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h1) );
    assert!( approx_eq(n_h1_h1,n_h1_h2) );

    let [n_h, n_h1, n_h2, n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2,_n_mol] =
      get_conditional_hydron_stats(&structure,&filter_nx);

    assert_eq!(n_h,n_uc*n_nx);
    assert!(n_h1_h1>0);
    assert!(n_h2_h2>0);
    assert_eq!(n_h2_h1,0);
    assert_eq!(n_h1_h2,0);
    assert!( approx_eq(n_h1,n_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h2) );


   
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  fn get_conditional_hydron_stats(structure: &Structure,filter: &ParticleFilter)
    -> [usize;8]
  {

    let mut n_h = 0;
    let mut n_h1 = 0;
    let mut n_h2 = 0;
    let mut n_h1_h1 = 0;
    let mut n_h1_h2 = 0;
    let mut n_h2_h1 = 0;
    let mut n_h2_h2 = 0;
    let mut n_mol = 0;

    let mut is_ref_proton = true;
    let mut ref_mol_id = 1000000;
    let mut ref_cell_id = 1000000;

    let indices = filter.filter(structure);

    for idx in indices{
      let particle = &structure.bath_particles[idx];
      if (*particle).element == Element::Hydrogen{
        n_h += 1;
      }else{
        assert!((*particle).element != Element::Oxygen);
      }
      let mol_id = structure.molecule_id(idx);
      let idx0 = structure.primary_cell_indices[idx];
      let cell_indices = &structure.cell_indices[idx0];
      let cell_id = cell_indices.iter().position(|&x| {
            match x{ Some(id) => id == idx, None => false,}
          }).unwrap();
      let is_proton = (*particle).isotope == Isotope::Hydrogen1;
      let is_deuteron = (*particle).isotope == Isotope::Hydrogen2;
      
      if is_proton{
        n_h1 += 1;
      }else if is_deuteron{
        n_h2 += 1;
      }

      if !is_proton && !is_deuteron{ continue; }

      if mol_id != ref_mol_id || cell_id != ref_cell_id{
        ref_mol_id = mol_id;
        ref_cell_id = cell_id;
        is_ref_proton = is_proton;
        n_mol += 1;
        continue;
      }

      match (is_proton,is_ref_proton){
        (true,true) => n_h1_h1 +=1, 
        (true,false) => n_h1_h2 +=1, 
        (false,true) => n_h2_h1 +=1, 
        (false,false) => n_h2_h2 +=1, 
      }

    }
    [n_h, n_h1, n_h2, n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2, n_mol]
  }

  //----------------------------------------------------------------------------
  // TODO: set non-exchangeable hydrons as cosubstitution groups.
  // TODO: set extra-cell non-protons to voids.
  #[test]
  fn test_build_extended_structure_tempo_3gly_1npr(){
    let filename = "./assets/TEMPO_3gly_1npr_50A.pdb";
    let mut structures = pdb::parse_pdb(&filename).unwrap();
    assert_eq!(structures[0].bath_particles.len(),13891);

    let mut particle_configs =vec![
      ParticleConfig::new("exchangeable_hydrogens".to_string()),
      ParticleConfig::new("nonexchangeable_hydrogens".to_string()),
    ];

    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Hydrogen];
    filter.bonded_elements = vec![Element::Oxygen];
    particle_configs[0].filter = Some(filter.clone());

    filter.bonded_elements = vec![];
    filter.not_bonded_elements = vec![Element::Oxygen];
    particle_configs[1].filter = Some(filter);
    
    let mut properties = ParticleProperties::new();
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: 0.5});
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 0.5});

    //properties.isotopic_distribution.extracell_void_probability = Some(0.5);


    particle_configs[0].properties = Some(properties.clone());
    
    particle_configs[1].properties = Some(properties);
    
    let mut config = Config::new();
    config.load_geometry = Some(LoadGeometry::Cube);
    config.particles = particle_configs;
    config.radius = Some(73.5676e-10);
    let n_uc = 125;


    structures[0].build_primary_structure(&config).unwrap();

    let mut rng =  ChaCha20Rng::from_entropy();

    let mut structure = structures[0].clone();
    structure.build_extended_structure(&mut rng, &config);

    // Test: set_cell_shifts().
    assert_eq!(structure.cell_offsets.len(),n_uc);

    // Test: extend_structure().
    //assert!(structure.bath_particles.len()-n_uc*num_particles);

    // Test: set_isotopologue().
    let mut n_h1 = 0;
    let mut n_h2 = 0;
    for particle in structure.bath_particles.iter(){
      if (*particle).isotope == Isotope::Hydrogen1{
        n_h1 += 1;
        assert_eq!((*particle).element,Element::Hydrogen);

      }else if (*particle).isotope == Isotope::Hydrogen2{
        n_h2 += 1; 
        assert_eq!((*particle).element,Element::Hydrogen);
      }
    }

    let n_hydrogens = 18 + 775*8 + 251*8;
    let n_hydrogens_uc = (n_uc)*n_hydrogens;
    let n_up = 3*n_hydrogens_uc/5;  
    let n_low = 2*n_hydrogens_uc/5;  
    assert!(n_h1 >= n_low);
    assert!(n_h1 <= n_up);
    assert!(n_h2 >= n_low);
    assert!(n_h2 <= n_up);
    assert_eq!(n_h1+n_h2,n_hydrogens_uc);
  }
}
