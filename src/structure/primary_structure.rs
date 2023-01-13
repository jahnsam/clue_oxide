use crate::structure::Structure; 
use crate::config::Config;

impl Structure{
  /// This method uses an input `Config` to set the structure's
  /// spins and exchange groups.  The number of bath particle is unchanged.
  pub fn build_primary_structure(&mut self, config: &Config){

    // TODO: one PBC on each side should be used to ensure reconect_bonds()
    // does not take a spin out of range.
    // self.reconect_bonds();

    self.set_spins(config);

    // Set methyls.
  }  
  //----------------------------------------------------------------------------
  // This method self.bath_spins_indices so that each element corresponds
  // to an element of self_bath_particles that has the potential to have a spin.
  fn set_spins(&mut self,config: &Config){
    self.pair_particle_configs(config);

    // Find spins.
    let n_spins = self.count_spins(config);
    self.bath_spins_indices = Vec::<usize>::with_capacity(n_spins);

    for (idx, particle) in self.bath_particles.iter().enumerate(){

      if particle.isotope.spin_multiplicity() <= 1 { continue; }

      if let Some(id) = self.particle_config_ids[idx]{
        if config.max_spin_multiplicity_for_particle_config(id) <= 1{
          continue;
        }
      }
      self.bath_spins_indices.push(idx);
    }

   
  }
  //----------------------------------------------------------------------------
  // This method counts the to an elements of self_bath_particles that have
  // the potential to have a spin.
  fn count_spins(&self,config: &Config) -> usize{
    let mut n_spins = 0;
    for (idx, particle) in self.bath_particles.iter().enumerate(){

      if particle.isotope.spin_multiplicity() <= 1 { continue; }

      if let Some(id) = self.particle_config_ids[idx]{
        if config.max_spin_multiplicity_for_particle_config(id) <= 1{
          continue;
        }
      }
      n_spins += 1;
    }
    n_spins
  }
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  /*
  // The method uses the PBCs to move atoms near atoms they are bonded to.
  fn reconnect_bonds(&mut self){

    const MAXBOND: f64 = 3.0;
    for connections in self.connections.iter() {

      let mut it = connections.indices.iter();
      let idx0: &usize = it.next().unwrap();

      let r0 = self.bath_particles[*idx0].coordinates.clone();

      for idx in it{

        let r = self.bath_particles[*idx].coordinates.clone();

        if (&r-&r0).magnitude() < MAXBOND { continue; }

        self.bath_particles[*idx].coordinates.x()
         = minimize_absolute_difference_for_step(r.x, r0.x, self.crystal.a);

        self.bath_particles[*idx].coordinates.y()
         = minimize_absolute_difference_for_step(r.y, r0.y, self.crystal.b);

        self.bath_particles[*idx].coordinates.z()
         = minimize_absolute_difference_for_step(r.z, r0.z, self.crystal.c);
      }





    }
 }
  */
 //-----------------------------------------------------------------------------

}

//------------------------------------------------------------------------------
// TODO: find better module for this function.
// Let n be an integer and x,x0,x1, and step be real numbers, 
// This function returns the x1 := x + n*step founds by
//
//   x1 = argmin(|x0-x1|) = argmin( |x0 - (x + n*step)| ).
//
fn minimize_absolute_difference_for_step( x: f64, x0: f64, step: f64 )
 -> f64 {

   assert!(step > 0.0);

   let mut x1 = x;

   loop{

     if (x1 - x0).abs() <= (x1 + step - x0).abs()
       && (x1 - x0).abs() <= (x1 - step - x0).abs(){
         break;
       }
     else if (x1 - x0).abs() > (x1 + step - x0).abs(){
      x1 += step;
     }
     else if (x1 - x0).abs() > (x1 - step - x0).abs(){
      x1 -= step;
     }
   }

  x1
}
