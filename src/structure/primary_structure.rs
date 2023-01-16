use crate::structure::Structure; 
use crate::config::Config;
use crate::space_3d;

impl Structure{
  /// This method uses an input `Config` to set the structure's
  /// spins and exchange groups.  The number of bath particle is unchanged.
  fn build_primary_structure(&mut self, config: &Config){

    // TODO: one PBC on each side should be used to ensure reconect_bonds()
    // does not take a spin out of range.
    self.reconnect_bonds();

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
  // The method uses the PBCs to move atoms near atoms they are bonded to.
  fn reconnect_bonds(&mut self){

    // TODO: Decide on better error handeling here.
    match self.cell_offsets.len(){
      0 => return,
      3 => (),
      _ => panic!("There should be 3 cell offsets."),  
    }

    const MAXBOND: f64 = 3.0; // Angstroms.

    //for connections in self.connections.iter() {
    for idx0 in 0..self.bath_particles.len() {


      // Get the list of bonded indices, if there are any.
      let connections: &Vec::<usize>;
      if let Some(cnc)= self.connections.get_neighbors(idx0){
        connections = cnc;
      }else{
        continue;
      }

      // Get the coordinates of the atom of interest.
      let r0 = self.bath_particles[idx0].coordinates.clone();

      for idx in connections.iter(){

        // Get the coordinates of the bonded atom.
        let r = self.bath_particles[*idx].coordinates.clone();

        // Check if the atoms are near each other.
        let delta_r = &r -&r0;
        if delta_r.norm() < MAXBOND { continue; }

        // Loop through all 3 spatial dimensions.
        for ix in 0..3 {

          // Find the PBC copy of the neghbor that is closest.
          let r1 = space_3d::minimize_absolute_difference_for_vector3d_step(
            &r,&r0,&self.cell_offsets[ix]);

          self.bath_particles[*idx].coordinates = r1;
        }
      }

    }
 }
 //-----------------------------------------------------------------------------

}

