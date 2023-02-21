use crate::clue_errors::CluEError;
use crate::config::{Config,
  particle_config::{ParticleConfig,ParticleProperties}};
use crate::space_3d::Vector3D;
use crate::structure::{Structure,particle::Particle,
  particle_filter::ParticleFilter};
use crate::math;  

use rand::seq::index::sample;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Binomial, Distribution, Uniform};
use std::collections::HashMap;


impl Structure{
  pub fn build_extended_structure(&mut self, rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    // Construct offset vectors for each pbc application.
    self.set_cell_shifts(config)?;

    // Copy non-voidable particles over to each PBC. 
    self.extend_structure(config)?;

    // TODO: reduce calls to pair_particle_configs().
    self.pair_particle_configs(config);

    self.find_cosubstitution_groups(config); // TODO: where should this line go?

    // For particles that are to either be kept or dropped entirely,
    // determine which particles to keep.
    self.add_voidable_particles(rng, config)?;

    self.pair_particle_configs(config);// TODO: where should this line go?

    // Set isotpic identities after adding voidable particles since the
    // non-void particles can potentially have multiple isotopic options. 
    self.set_isotopologue(rng, config)?;


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

    for (cosub_idx,cosubstitution_group) in self.cosubstitution_groups.iter()
      .enumerate()
    {
      if cosubstitution_group.is_empty(){ continue; }

      // Check if this particle has a custom config.
      let Some(config_id) = self.particle_config_ids[0] else {
        continue;
      };

      // Check if this particle has any custom properties.
      let Some(properties) = &particle_configs[config_id].properties else{
        continue;
      };

      let Some(p) = properties.isotopic_distribution
        .extracell_void_probability else{
          continue;
      };
      let p = 1.0 - p;
      let ppt = (p*1e12) as u64;
      if let Some(indices) = property_map.get_mut(&ppt){
        indices.push(cosub_idx);
      }else{
        let indices = vec![cosub_idx];
        property_map.insert(ppt,indices);
      }

    }

    
    for offset in self.cell_offsets.iter(){
      for (ppt, indices) in property_map.iter(){ 

        let p = (*ppt as f64)/1e12; 
        let n_groups = indices.len();
        let Ok(bin) = Binomial::new(n_groups as u64, p) else {
          return Err(CluEError::CannotSampleBinomialDistribution(n_groups,p));
        };

        let n_indices = bin.sample(rng) as usize;
        let non_void_indices = sample(rng,n_groups,n_indices);

        for cosub_idx in non_void_indices.iter(){
          for particle_idx in self.cosubstitution_groups[cosub_idx].iter(){

            let mut new_particle = self.bath_particles[*particle_idx].clone();
            new_particle.coordinates = &new_particle.coordinates + offset;
            self.bath_particles.push(new_particle);
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

    // TODO: loop over cosubstitution_groups not individual particles.
    // Loop over bath particle.
    for cosubstitution_group in self.cosubstitution_groups.iter(){
      // Generate a random number.
      let random_number = range.sample(rng);

      for particle_idx in cosubstitution_group.iter(){
        
        let particle = &mut self.bath_particles[*particle_idx];
      
        // Check if this particle has a custom config.
        let Some(config_id) = self.particle_config_ids[*particle_idx] else {
          continue;
        };

        // Check if this particle has any custom properties.
        let Some(properties) = &particle_configs[config_id].properties else{
          continue;
        };


        let mut cdf = 0.0;
     
        // Select a isotope for this particle.
        for iso in properties.isotopic_distribution.isotope_abundances.iter(){
          cdf += iso.abundance;
          if cdf >= random_number{
            particle.isotope = iso.isotope;
            break;
          } 
        }

      }
    }

    Ok(())

  }
  //----------------------------------------------------------------------------
  fn extend_structure(&mut self, config: &Config) -> Result<(),CluEError>{

    let mut extra_particles = Vec::<Particle>::with_capacity(
        self.number()*self.cell_offsets.len());


    for offset in self.cell_offsets.iter(){

      for (particle_idx,particle) in self.bath_particles.iter().enumerate(){

        // TODO: implement force_no_pbc?
        // No because properties.extracell_void_probability = Some(0.0)
        // does the same thing
        // if properties.force_no_pbc{}

        
      // Check if this particle has a custom config.
      if let Some(config_id) = self.particle_config_ids[particle_idx]{

      // Check if this particle has any custom properties.
      if let Some(properties) = &config.particles[config_id].properties{
        // TODO: check void probability.
        if let Some(_) = properties.isotopic_distribution
          .extracell_void_probability{
          continue;
        }
      }
      }
      
        let mut new_particle = (*particle).clone();
        new_particle.coordinates = &new_particle.coordinates + offset;
        extra_particles.push(new_particle);
      }
    }

    self.bath_particles = extra_particles;

    Ok(())
  }
  //----------------------------------------------------------------------------
  fn set_cell_shifts(&mut self, config: &Config) -> Result<(),CluEError>{
    let cell_edges = self.cell_offsets.clone();

    // TODO: Decide on better error handeling here.
    assert_eq!(cell_edges.len(),3);

    self.cell_offsets = Structure::build_cell_shifts(cell_edges, config)?;
    Ok(())
  }
  //----------------------------------------------------------------------------
  fn build_cell_shifts(cell_edges: Vec::<Vector3D>, config: &Config) 
    -> Result<Vec::<Vector3D>,CluEError>
  {
    
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

          let offset = &(&cell_edges[0].scale(ix as f64)
            + &cell_edges[1].scale(iy as f64)) 
            + &cell_edges[2].scale(iz as f64);

          cell_offsets.push(offset);
    }}}
    cell_offsets.sort();
    Ok(cell_offsets)
  }

}



#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::{pdb, primary_structure,particle_filter::ParticleFilter};
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

  use crate::config::particle_config::{ParticleConfig,
    ParticleProperties,IsotopeDistribution,IsotopeAbundance};
  use crate::physical_constants::{Element,Isotope};

  #[test]
  fn test_build_extended_structure(){
    let filename = "./assets/TEMPO.pdb";
    let file = std::fs::read_to_string(filename).unwrap();
    let mut structures = pdb::parse_pdb(&file).unwrap();

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
    config.particles = particle_configs;
    config.radius = Some(73.5676e-10);


    structures[0].build_primary_structure(&config);
    let num_particles = structures[0].bath_particles.len();

    let mut rng =  ChaCha20Rng::from_entropy();

    let mut structure = structures[0].clone();
    structure.build_extended_structure(&mut rng, &config);

    // Test: set_cell_shifts().
    assert_eq!(structure.cell_offsets.len(),125);

    // Test: extend_structure().
    assert_eq!(structure.bath_particles.len(),125*num_particles);

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
    assert!(n_h1 >= 1000);
    assert!(n_h1 <= 1200);
    assert!(n_h2 >= 1000);
    assert!(n_h2 <= 1200);
    assert_eq!(n_h1+n_h2,125*18);
    
  }

}
