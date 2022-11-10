use super::pdb;
use matrix_oxide as mox;

struct Structure{

  // Core Infomation.
  number: usize,         
  coordinates: mox::Mat,
  elements: Vec::<Elements>
  exchange_group_ids: Vec::<usize>,
  tunnel_splittings: Vec::<f64>,

  // Identifying informtion
  residues: Vec::<String>,
  pdb_ids: Vec::<usize>,

}


pub fn build_primary_structure(pdb: &PDB) -> Structure {

  let pdb = center_pdb(pdb);

  let mut primary_structure = PrimaryStructure::with_capacity(pdb.number());
  let mut methyl_counter = 0;

  let exchange_groups = Vec::<usize>::with_capacity(pdb.number());
  let tunnel_splittings = Vec::<f64>::with_capacity(pdb.number());

  for ipdb in 0..pdb.number(){
    
    let element = get_element(pdb.element(ii));

    if element == Element::Carbon{
      // Flag the hydrogens.
      match get_associated_hydrogens() {
        Some(hydrogens) => {
          methyl_counter += 1;
          for ih in hydrogens{
            exchange_groups[ih] = methyl_counter;
          }
        },
        None => continue,
      }
    }

    if get_max_spin_multiplicity_for_any_isotope(element,config) == 0 {
      continue;
    }


    primary_structure.push_element(element, &pdb)

  }

  primary_structure.set_tunnel_splittings(tunnel_splittings);
  primary_structure.set_number_exchange_groups(methyl_counter);
  primary_structure.set_exchange_groups(exchange_groups);
  primary_structure.trim();
  primary_structure.set_pdb(pdb);

  primary_structure
}



pub fn build_extended_structure(){

  // Determine how many unit cells are needed.

  // Find void particles.

  // Copy non-void particles
}


fn get_element(element_str: &str) -> Option<Element>{

  match element_str.trim() {
  
    "e" => return Some(Element::Electron),
    "H" => return Some(Element::Hydrogen),
    "D" => return Some(Element::Deuterium),
    "He" => return Some(Element::Helium),
    "Li" => return Some(Element::Lithium),
    "B" => return Some(Element::Boron),
    "C" => return Some(Element::Carbon),
    "N" => return Some(Element::Nitrogen),
    "O" => return Some(Element::Oxygen),
    "F" => return Some(Element::Fluorine),
    "Ne" => return Some(Element::Neon),

  }


}
