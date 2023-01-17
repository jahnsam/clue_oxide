use crate::structure::{Structure, exchange_groups::*}; 
use crate::config::Config;
use crate::space_3d;
use crate::physical_constants::Element;

impl Structure{
  /// This method uses an input `Config` to set the structure's
  /// spins and exchange groups.  The number of bath particle is unchanged.
  fn build_primary_structure(&mut self, config: &Config){

    // TODO: one PBC on each side should be used to ensure reconect_bonds()
    // does not take a spin out of range.
    self.reconnect_bonds();

    self.set_spins(config);

    // Set methyl and primary amonium groups..
    self.set_exchange_groups();
  }  
  //----------------------------------------------------------------------------
  // This method sets self.bath_spins_indices so that each element corresponds
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
 fn set_exchange_groups(&mut self){

   let exchange_groups = self.find_exchange_groups();
   
   let mut exchange_group_ids 
     = Vec::<Option<usize>>::with_capacity(self.number());
   for ii in 0..self.number(){
     exchange_group_ids.push(None);
   }

   let mut exchange_coupling = Vec::<f64>::with_capacity(exchange_groups.len());
   for ii in 0..exchange_groups.len(){
     exchange_coupling.push(0.0);
     for h in exchange_groups[ii].indices(){
       exchange_group_ids[h] = Some(ii);
     }
   }

   self.exchange_groups = Some( ExchangeGroupManager{
     exchange_groups,
     exchange_group_ids,
     exchange_coupling,
   });
 }  
 //-----------------------------------------------------------------------------
 fn find_exchange_groups(&self) -> Vec::<ExchangeGroup> {
 
   // Five atoms are required to form a methyl group.
   let n_max: usize = self.number()/5;

   let mut exchange_groups = Vec::<ExchangeGroup>::with_capacity(n_max);

   for ii in 0..self.number(){

     let hydrogens = self.get_methyl_hydrogen_indices(ii);
      if let Some([h0,h1,h2]) = hydrogens {
        let r_carbon = self.bath_particles[ii].coordinates.clone();
        let r_h0 = self.bath_particles[h0].coordinates.clone();
        let r_h1 = self.bath_particles[h2].coordinates.clone();
        let r_h2 = self.bath_particles[h2].coordinates.clone();

        
        if self.bath_particles[ii].element == Element::Carbon{
          exchange_groups.push(ExchangeGroup::Methyl(
                C3Rotor::from(r_carbon, r_h0, r_h1, r_h2, [h0,h1,h2]) ));

        }else if self.bath_particles[ii].element == Element::Nitrogen{
          exchange_groups.push(ExchangeGroup::PrimaryAmonium(
                C3Rotor::from(r_carbon, r_h0, r_h1, r_h2, [h0,h1,h2]) ));
        }
      } 
   }

   exchange_groups
 }  
 //-----------------------------------------------------------------------------
 fn get_methyl_hydrogen_indices(&self, index: usize) ->Option<[usize;3] >{
   
   let particle = &self.bath_particles[index];

   if particle.element != Element::Carbon 
    && particle.element != Element::Nitrogen {
     return None;
   }

   // Get the list of bonded indices, if there are any.
   let connections: &Vec::<usize>;
   if let Some(cnc)= self.connections.get_neighbors(index){
     connections = cnc;
   }else{
     return None;
   }
   
   // Methyls and primary amoniums have four atoms bonded to the central
   // carbon and nitrogen respectively.
   if self.connections.len() != 4 { 
     return None;
   }

   let mut hydrons = Vec::<usize>::with_capacity(3);

   for idx in connections.iter() {

     if self.bath_particles[*idx].element == Element::Hydrogen{
       
       // Skip methane and amonium, which are not implemented. 
       if hydrons.len()==3{ return None;}

       hydrons.push(*idx);
     }
   }


   if hydrons.len()!=3{ return None;}

   Some([hydrons[0], hydrons[1], hydrons[2] ])
 }  
 //-----------------------------------------------------------------------------

}

