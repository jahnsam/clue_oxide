use crate::clue_errors::CluEError;
use crate::config::{Config,
  particle_config::{ParticleConfig,ParticleProperties}};
use crate::space_3d::Vector3D;
use crate::structure::{Structure,particle::Particle,
  particle_filter::ParticleFilter};

use rand::seq::index::sample;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Binomial, Distribution, Uniform};


impl Structure{
  fn build_extended_structure(&mut self, rng: &mut ChaCha20Rng, config: &Config) 
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


    Ok(())
  }
  //----------------------------------------------------------------------------
  // TODO: implement.
  fn add_voidable_particles(&mut self, 
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    // TODO: TEST CODE, DOES NOT DO ANYTHING YET.
    let n: usize = 20;
    let bin = Binomial::new(n as u64, 0.3).unwrap();
    let v = bin.sample(&mut rand::thread_rng());
    println!("Choosing {} out of {} from a binomial distribution", v,n);

    let x = sample(rng,n, v as usize);
    println!("{:?}",x);

    Ok(())

  }
  //----------------------------------------------------------------------------
  fn set_isotopologue(&mut self, rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {
    
    let range = Uniform::new(0.0f64, 1.0);
    let particle_configs: &Vec::<ParticleConfig>;
    match &config.particles{
      Some(part_confs) => particle_configs = part_confs,
      None => return Ok(()),
    }

    for (particle_idx,particle) in self.bath_particles.iter_mut().enumerate(){
      
      let config_id: usize;
      match self.particle_config_ids[particle_idx]{
        Some(id) => config_id = id,
        None => continue,
      }

      let properties: &ParticleProperties;
      match &particle_configs[config_id].properties{
        Some(props) => properties = props,
        None => continue,
      }


      let random_number = range.sample(rng);
      let mut cdf = 0.0;
      for iso in properties.isotopic_distribution.isotope_abundances.iter(){
        cdf += iso.abundance;
        if cdf >= random_number{
          particle.isotope = iso.isotope;
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

      for particle in self.bath_particles.iter(){

        // TODO: implement force_no_pbc?
        // if properties.force_no_pbc{}

        // TODO: check void probability.
        // if properties.extracell_void_probability{}
        let mut new_particle = (*particle).clone();
        new_particle.coordinates = &new_particle.coordinates + offset;
        extra_particles.push(new_particle);
      }
    }

    self.pair_particle_configs(config);

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

    let radius: f64;
    match config.radius{
      Some(r) => radius = r,
      None => return Err(CluEError::NoRadius),
    }


    let mut n_cells = 1;
    for ix in 0..3{
      n_cells_per_dim[ix] = std::cmp::max(1,
       ceil(radius/cell_edges[ix].norm()) as i32
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

// TODO: move ceil() to a different module.
pub fn ceil(x: f64) -> f64{

  let mut a = x as i32;

  let err = (x-a as f64).abs();
  if err > 1e-12 && x >= 0.0{
    a += 1;
  }

  a as f64
}

#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_ceil(){
    assert_eq!(ceil(2.0),2.0);
    assert_eq!(ceil(2.1),3.0);
    assert_eq!(ceil(-2.1),-2.0);
    assert_eq!(ceil(-2.0),-2.0);
  }
}
