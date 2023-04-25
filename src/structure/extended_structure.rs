use crate::clue_errors::CluEError;
use crate::config::{Config,LoadGeometry};
use crate::space_3d::Vector3D;
use crate::structure::Structure;
use crate::structure::exchange_groups::*;
use crate::math;  

use rand::seq::index::sample;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Binomial, Distribution, Uniform};
use std::collections::HashMap;
use crate::config::particle_config::IsotopeDistribution;

impl Structure{
  /// This function takes a structure from build_primary_structure(), applies
  /// periodic boundary condition (if applicable), and sets isotope identities.
  pub fn build_extended_structure(&mut self,
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    // Construct offset vectors for each pbc application.
    self.set_cell_shifts(config)?;

    // Copy non-voidable particles over to each PBC. 
    self.extend_structure(config)?;

    // For particles that are to either be kept or dropped entirely,
    // determine which particles to keep.
    self.add_voidable_particles(rng, config)?;

    // Set isotpic identities after adding voidable particles since the
    // non-void particles can potentially have multiple isotopic options. 
    self.set_isotopologue(rng, config)?;

    self.update_exhange_groups(config)?;

    self.trim_system(config)?;

    Ok(())
  }
  //----------------------------------------------------------------------------
  
  fn trim_system(&mut self, config: &Config) -> Result<(),CluEError>{
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

          if particle.coordinates.norm() > radius{
            particle.active = false;
          }
        }
      },
    }

    Ok(())
  }
  
  //----------------------------------------------------------------------------
  fn update_exhange_groups(&mut self, config: &Config) -> Result<(),CluEError>
  {
    let Some(exchange_group_manager0) = &self.exchange_groups else {
      return Ok(());
    };

    let n_ex = exchange_group_manager0.exchange_groups.len()
      *self.cell_offsets.len();
    
    let mut exchange_groups = Vec::<ExchangeGroup>::with_capacity(n_ex);
    let mut exchange_group_ids = self.bath_particles.iter().map(|p| None)
      .collect::<Vec::<Option<usize>>>();
    let mut exchange_couplings = Vec::<f64>::with_capacity(n_ex);

    // Loop over unit cells.
    for icell in 0..self.cell_offsets.len(){ 
      for (ex_id, exchange_group_0) in exchange_group_manager0
        .exchange_groups.iter().enumerate()
      {
        let indices_0 = exchange_group_0.indices();
        let mut indices = Vec::<usize>::with_capacity(indices_0.len());

        for &idx0 in indices_0.iter(){
          let Some(idx) = self.cell_indices[idx0][icell] else{
            break;
          };
          indices.push(idx);
        }
        if indices.len() != indices_0.len(){
          continue;
        }

        let center = exchange_group_0.centroid() + &self.cell_offsets[icell];
        let exchange_group: ExchangeGroup;
        match exchange_group_0{
          ExchangeGroup::Methyl(rotor) 
          => exchange_group = ExchangeGroup::Methyl(
              C3Rotor{ center, normal: rotor.normal.clone(), 
              indices: [indices[0],indices[1],indices[2]]}),
          ExchangeGroup::PrimaryAmonium(rotor)
          => exchange_group = ExchangeGroup::PrimaryAmonium(
              C3Rotor{ center, normal: rotor.normal.clone(), 
              indices: [indices[0],indices[1],indices[2]]}),
        }
        exchange_groups.push(exchange_group);


        for (ii,&idx) in indices.iter().enumerate(){
          exchange_group_ids[idx] = exchange_group_manager0
            .exchange_group_ids[ii].clone();
        }


        let exchange_coupling: f64 
          = self.extract_exchange_coupling(indices[0],config);
        exchange_couplings.push(exchange_coupling);

     
      }
    }
    self.exchange_groups = Some(ExchangeGroupManager{
      exchange_groups,
      exchange_group_ids,
      exchange_couplings,
    });
    
    Ok(())
  }
  //----------------------------------------------------------------------------
  // TODO: TEST add_voidable_particles()
  fn add_voidable_particles(&mut self, 
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    /*
    let n: usize = 20;
    let bin = Binomial::new(n as u64, 0.3).unwrap();
    let v = bin.sample(&mut rand::thread_rng());
    println!("Choosing {} out of {} from a binomial distribution", v,n);

    let x = sample(rng,n, v as usize);
    println!("{:?}",x);
    */
    
    // Check if there are any configurations.
    let particle_configs = &config.particles;

    let mut property_map = HashMap::<u64, Vec::<usize>>::with_capacity(
        self.cosubstitution_groups.len());

    // Find all the cosubstitution_groups with the same substitution chance.
    for (cosub_idx,cosubstitution_group) in self.cosubstitution_groups.iter()
      .enumerate()
    {
      if cosubstitution_group.is_empty(){ continue; }

      let particle_idx0 = cosubstitution_group[0];
      if !self.bath_particles[particle_idx0].active {continue;}
      
      // Check if this particle has a custom config.
      let Some(config_id) = self.particle_config_ids[particle_idx0] else {
        continue;
      };

      // Check if this particle has any custom properties.
      let Some(properties) = &particle_configs[config_id].properties else{
        continue;
      };

      let Some(p_remove) = properties.isotopic_distribution
        .extracell_void_probability else{
          continue;
      };
      let p_keep = 1.0 - p_remove;
      let ppt = (p_keep*1e12) as u64;
      if let Some(indices) = &mut property_map.get_mut(&ppt){
        indices.push(cosub_idx);
      }else{
        let indices = vec![cosub_idx];
        property_map.insert(ppt,indices);
      }

    }

    
    for (cell_id,offset) in self.cell_offsets.iter().enumerate().skip(1){
      for (ppt, indices) in property_map.iter(){ 

        let p_keep = (*ppt as f64)/1e12; 
        let n_groups = indices.len();
        let Ok(bin) = Binomial::new(n_groups as u64, p_keep) else {
          return Err(CluEError::CannotSampleBinomialDistribution(
                n_groups,p_keep));
        };

        let n_indices = bin.sample(rng) as usize;
        let non_void_indices = sample(rng,n_groups,n_indices);

        for cosub_idx in non_void_indices.iter(){
          for particle_idx0 in self.cosubstitution_groups[cosub_idx].iter(){

            let mut new_particle = self.bath_particles[*particle_idx0].clone();
            new_particle.coordinates = &new_particle.coordinates + offset;
            self.bath_particles.push(new_particle);
        
            self.primary_cell_indices.push(*particle_idx0);
            let particle_idx = self.bath_particles.len() - 1;
            self.cell_indices[*particle_idx0][cell_id] = Some(particle_idx);
          }
        }

      }
    } 
    Ok(())

  }
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
  fn extend_structure(&mut self, config: &Config) -> Result<(),CluEError>{

    let n_particles0 = self.number();

    let n_to_reserve = n_particles0*(self.cell_offsets.len()-1);
    self.bath_particles.reserve( n_to_reserve );
    self.primary_cell_indices.reserve( n_to_reserve );

    for cindices in self.cell_indices.iter_mut(){
      cindices.reserve(self.cell_offsets.len()-1);
    }

    for offset in self.cell_offsets.iter().skip(1){

      for particle_idx0 in 0..n_particles0{
        if !self.bath_particles[particle_idx0].active {
          self.cell_indices[particle_idx0].push(None);
          continue;
        }

        // TODO: implement force_no_pbc?
        // No because properties.extracell_void_probability = Some(0.0)
        // does the same thing
        // if properties.force_no_pbc{}

        
        // Check if this particle has a custom config.
        if let Some(config_id) = self.particle_config_ids[particle_idx0]{
          // Check if this particle has any custom properties.
          if let Some(properties) = &config.particles[config_id].properties{
            // TODO: check void probability.
            if let Some(_) = properties.isotopic_distribution
              .extracell_void_probability{
              self.cell_indices[particle_idx0].push(None); 
              continue;
            }
          }
        }
      
        let mut new_particle = self.bath_particles[particle_idx0].clone();
        new_particle.coordinates = &new_particle.coordinates + offset;
        self.bath_particles.push(new_particle);


        self.primary_cell_indices.push(particle_idx0);
        let particle_idx = self.bath_particles.len() - 1;
        self.cell_indices[particle_idx0].push(Some(particle_idx));
      }
    }


    Ok(())
  }
  //----------------------------------------------------------------------------
  fn set_cell_shifts(&mut self, config: &Config) -> Result<(),CluEError>{
    let cell_edges = self.cell_offsets.clone();

    if cell_edges.is_empty(){
      self.cell_offsets = vec![Vector3D::zeros()];
      println!("No unit cell information found.  
Periodic boudary conditions will not be applied.");
      return Ok(());
    }else if cell_edges.len() != 3 {
      return Err(CluEError::IncorrectNumberOfAxes(cell_edges.len(),3));
    }

    self.cell_offsets = build_cell_shifts(cell_edges, config)?;
    Ok(())
  }
}
  //----------------------------------------------------------------------------
  fn build_cell_shifts(cell_edges: Vec::<Vector3D>, config: &Config) 
    -> Result<Vec::<Vector3D>,CluEError>
  {
    
    /*
    let Some(load_geometry) = &config.load_geometry else {
      return Err(CluEError::NoLoadGeometry);
    };
    */
    let mut n_cells_per_dim = [1,1,1];

    let Some(radius) = config.radius else {
      return Err(CluEError::NoRadius);
    };


    let mut n_cells = 1;
    for ix in 0..3{
      n_cells_per_dim[ix] = std::cmp::max(1,
       math::ceil(radius/cell_edges[ix].norm()) as i32
          );

      n_cells *= 1 + 2*(n_cells_per_dim[ix] as usize);
    }

    let mut cell_offsets = Vec::<Vector3D>::with_capacity(n_cells);
  
    for ix in -n_cells_per_dim[0]..=n_cells_per_dim[0]{
      for iy in -n_cells_per_dim[1]..=n_cells_per_dim[1]{
        for iz in -n_cells_per_dim[2]..=n_cells_per_dim[2]{

          /*
          // For spherical geometries, do not use unit cell that cannot
          // have any particles within the load radius.
          if ix !=0 && iy != 0 && iz != 0 
            && *load_geometry == LoadGeometry::Sphere{
              let jx = (ix.abs() - 1)*ix.signum();
              let jy = (iy.abs() - 1)*iy.signum();
              let jz = (iz.abs() - 1)*iz.signum();
              
            let min_possible_r = &(&cell_edges[0].scale(jx as f64)
              + &cell_edges[1].scale(jy as f64)) 
              + &cell_edges[2].scale(jz as f64);
            
            if min_possible_r.norm() > radius{
              continue;
            }

          }
          */

          let offset = &(&cell_edges[0].scale(ix as f64)
            + &cell_edges[1].scale(iy as f64)) 
            + &cell_edges[2].scale(iz as f64);

          cell_offsets.push(offset);
    }}}
    cell_offsets.sort();
    Ok(cell_offsets)
  }




#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb;
  use crate::structure::particle_filter::*;
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
  use crate::config::DetectedSpinCoordinates;

  use crate::config::particle_config::{ParticleConfig,
    ParticleProperties,IsotopeDistribution,IsotopeAbundance};
  use crate::physical_constants::{Element,Isotope};

  use crate::config::LoadGeometry;
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_cell_shifts(){
    let cell_edges = vec![
      Vector3D::from([1.0, 0.0, 0.0]),
      Vector3D::from([0.0, 0.5, 0.0]),
      Vector3D::from([0.0, 0.0, 2.0]),
    ];
    let mut config = Config::new();
    config.radius = Some(1.0);
    let cell_shifts = build_cell_shifts(cell_edges, &config).unwrap(); 
    assert_eq!(cell_shifts.len(), 3*5*3);
    assert_eq!(cell_shifts, vec![
      Vector3D::from([0.0, 0.0, 0.0]),
      Vector3D::from([0.0, -0.5, 0.0]),
      Vector3D::from([0.0, 0.5, 0.0]),
      Vector3D::from([-1.0, 0.0, 0.0]),
      Vector3D::from([0.0, -1.0, 0.0]),
      Vector3D::from([0.0, 1.0, 0.0]),
      Vector3D::from([1.0, 0.0, 0.0]),
      Vector3D::from([-1.0, -0.5, 0.0]),
      Vector3D::from([-1.0, 0.5, 0.0]),
      Vector3D::from([1.0, -0.5, 0.0]),
      Vector3D::from([1.0, 0.5, 0.0]),
      Vector3D::from([-1.0, -1.0, 0.0]),
      Vector3D::from([-1.0, 1.0, 0.0]),
      Vector3D::from([1.0, -1.0, 0.0]),
      Vector3D::from([1.0, 1.0, 0.0]),
      Vector3D::from([0.0, 0.0, -2.0]),
      Vector3D::from([0.0, 0.0, 2.0]),
      Vector3D::from([0.0, -0.5, -2.0]),
      Vector3D::from([0.0, -0.5, 2.0]),
      Vector3D::from([0.0, 0.5, -2.0]),
      Vector3D::from([0.0, 0.5, 2.0]),
      Vector3D::from([-1.0, 0.0, -2.0]),
      Vector3D::from([-1.0, 0.0, 2.0]),
      Vector3D::from([0.0, -1.0, -2.0]),
      Vector3D::from([0.0, -1.0, 2.0]),
      Vector3D::from([0.0, 1.0, -2.0]),
      Vector3D::from([0.0, 1.0, 2.0]),
      Vector3D::from([1.0, 0.0, -2.0]),
      Vector3D::from([1.0, 0.0, 2.0]),
      Vector3D::from([-1.0, -0.5, -2.0]),
      Vector3D::from([-1.0, -0.5, 2.0]),
      Vector3D::from([-1.0, 0.5, -2.0]),
      Vector3D::from([-1.0, 0.5, 2.0]),
      Vector3D::from([1.0, -0.5, -2.0]),
      Vector3D::from([1.0, -0.5, 2.0]),
      Vector3D::from([1.0, 0.5, -2.0]),
      Vector3D::from([1.0, 0.5, 2.0]),
      Vector3D::from([-1.0, -1.0, -2.0]),
      Vector3D::from([-1.0, -1.0, 2.0]),
      Vector3D::from([-1.0, 1.0, -2.0]),
      Vector3D::from([-1.0, 1.0, 2.0]),
      Vector3D::from([1.0, -1.0, -2.0]),
      Vector3D::from([1.0, -1.0, 2.0]),
      Vector3D::from([1.0, 1.0, -2.0]),
      Vector3D::from([1.0, 1.0, 2.0]),
    ]);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_extended_structure_lone_tempo(){
    let filename = "./assets/TEMPO.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();

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
    config.detected_spin_position = Some( 
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );


    config.set_defaults();
    structure.build_primary_structure(&config).unwrap();
    let num_particles = structure.bath_particles.len();

    let mut rng =  ChaCha20Rng::from_entropy();

    structure.build_extended_structure(&mut rng, &config).unwrap();

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
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();


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
    config.detected_spin_position = Some( 
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    let n_uc = 125;


    config.set_defaults();
    structure.build_primary_structure(&config).unwrap();


    

    let mut rng =  ChaCha20Rng::from_entropy();

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
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();


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
    let r_nitrogen=Vector3D::from([36.440*1e-10, 36.900*1e-10,  37.100*1e-10]);
    let r_oxygen = Vector3D::from([35.290*1e-10, 36.430*1e-10, 37.810*1e-10]);
    let r_e = (&r_nitrogen + &r_oxygen).scale(0.5);
    config.detected_spin_position = Some( 
        DetectedSpinCoordinates::XYZ(r_e) );


    config.set_defaults();
    structure.build_primary_structure(&config).unwrap();


    

    let mut rng =  ChaCha20Rng::from_entropy();

    structure.build_extended_structure(&mut rng, &config).unwrap();

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
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    assert_eq!(structure.bath_particles.len(),13891);

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
    config.detected_spin_position = Some( 
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    let n_uc = 125;


    config.set_defaults();
    structure.build_primary_structure(&config).unwrap();

    let mut rng =  ChaCha20Rng::from_entropy();

    structure.build_extended_structure(&mut rng, &config).unwrap();

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
