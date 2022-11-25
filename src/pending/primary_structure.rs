use super::pdb;
use super::Structure;
use matrix_oxide as mox;



pub fn build_primary_structure(pdb: &PDB) -> Structure {

  let pdb = center_pdb_central_spin(pdb);

  let mut primary_structure = PrimaryStructure::with_capacity(pdb.number());
  let mut methyl_counter = 0;

  
  let mut kept_atoms = vec![false; pdb.number()];
  let mut exchange_groups = vec![0 as usize; pdb.number()];
  let mut tunnel_splittings = Vec::<f64>::with_capacity(pdb.number());
  let mut exchange_group_sizes = Vec::<usize>::with_capacity(pdb.number());

  for ipdb in 0..pdb.number(){
    
    let element = Element::from(&pdb.element(ii));

    if element == Element::Carbon{
      // Flag the hydrogens.
      match get_associated_hydrogens() {
        Some(hydrogens) => {
          
          tunnel_splittings.push( get_tunnel_splitting(ipdb, &pdb, &config) );
          
          exchange_group_sizes.push(3);

          methyl_counter += 1;

          for ih in hydrogens{
            let idx_h = find_in_pdb_serial(ih, &pdb);
            exchange_groups[ih] = methyl_counter;
          }

        },
        None => (),
      }
    }

    if get_max_spin_multiplicity_for_any_isotope(element,config) == 1 {
      continue;
    }


    primary_structure.push_element(element, &pdb);
    kept_atoms[ipdb] = true;

  }

  primary_structure.set_tunnel_splittings(tunnel_splittings);
  primary_structure.set_number_exchange_groups(methyl_counter);
  primary_structure.set_exchange_groups(exchange_groups);
  primary_structure.set_exchange_group_ids(exchange_groups,kept_atoms);
  primary_structure.trim();
  primary_structure.set_pdb(pdb);

  primary_structure
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn get_tunnel_splitting(idx_pdb: usize, pdb: &PDB, config: &Config) -> f64{

  let specified_particle = specify_particle(idx_pdb,&pdb);

  let specified_index = find_specifying_index(&specified_particle, &config);

  config.particles.properties[specified_index].tunnel_splitting

}

//------------------------------------------------------------------------------
fn get_methyl_hydrogens(idx_pdb: usize, pdb: &PDB) 
  -> Option<[usize; 3]> {

    let mut hydrogen_indices: [usize; 3] = [0; 3];
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fn center_pdb_on_central_spin(pdb: &PDB, config: &Config) -> PDB{
  let origin = get_electron_coordinates(&pdb, &config);

  for ii in 0..pdb.number() {
    pdb.x[ii] -= origin[0];
    pdb.y[ii] -= origin[0];
    pdb.z[ii] -= origin[0];
  
  }

  pdb
}

//------------------------------------------------------------------------------
fn get_electron_coordinates(pdb: &PDB, config: &Config) -> Vec::<f64>{
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


