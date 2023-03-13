use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::physical_constants::{ANGSTROM,PI,Element,Isotope};
use crate::space_3d::Vector3D;
use crate::structure::{Structure, particle::Particle};

use pdbtbx;

pub fn parse_pdb(filename: &str) -> Result< Vec::<Structure>, CluEError>
{
  // Read PDB file.
  let Ok( (pdb, _errors) ) = pdbtbx::open(filename, 
      pdbtbx::StrictnessLevel::Medium) 
    else {
    return Err(CluEError::CannotOpenFile((*filename).to_owned()));
  };

  // Connections
  let connections = build_connections(&pdb)?;

  // Get cell offsets.
  let cell_offsets: Vec::<Vector3D>; 
  if let Some(unit_cell) = &pdb.unit_cell{
     cell_offsets = read_unit_cell(unit_cell)?;
  }else{
    cell_offsets = Vec::<Vector3D>::with_capacity(1);
  }

  let mut structures = Vec::<Structure>::with_capacity(pdb.model_count());  
  for model in pdb.models(){
    
    let bath_particles = parse_atoms(model)?;

      structures.push(Structure::new(
        bath_particles,
        connections.clone(),
        cell_offsets.clone(),
      ));
    }


  Ok(structures)
}
//------------------------------------------------------------------------------
fn read_unit_cell(unit_cell: &pdbtbx::UnitCell)
  -> Result<Vec::<Vector3D>,CluEError>
{
    let alpha = unit_cell.alpha()*PI/180.0;
    let beta = unit_cell.beta()*PI/180.0;
    let gamma = unit_cell.gamma()*PI/180.0;

    let a = unit_cell.a()*ANGSTROM;
    let b = unit_cell.b()*ANGSTROM;
    let c = unit_cell.c()*ANGSTROM;

    let a_vec = Vector3D::from([a, 0.0, 0.0]);
    let b_vec = Vector3D::from([b*f64::cos(gamma), b*f64::sin(gamma), 0.0]);
    let cx = c * f64::cos(beta);
    let cy = c*(f64::cos(alpha) - f64::cos(beta)*f64::cos(gamma))/f64::sin(gamma);
    let cz = c * f64::sqrt( 1.0 -cx*cx - cy*cy);
    let c_vec = Vector3D::from([cx, cy, cz]);


    return Ok(vec![a_vec,b_vec,c_vec]);
}
//------------------------------------------------------------------------------
fn parse_atoms(model: &pdbtbx::Model) -> Result<Vec::<Particle>,CluEError>
{
  let mut bath_particles = Vec::<Particle>::with_capacity(model.atom_count());
  let mut n_parse_failures = 0;

  for chain in model.chains(){
    for res in chain.residues(){   
      for conformer in res.conformers(){

        for atom in conformer.atoms(){

          let serial = atom.serial_number();
          let Some(elmt) = atom.element() else{
            n_parse_failures += 1;
            continue;
            //return Err(CluEError::AtomDoesNotSpecifyElement(serial))
          };

          let element = Element::from(&elmt.to_string())?;

          let coordinates = Vector3D::from(
              [atom.x()*ANGSTROM, atom.y()*ANGSTROM, atom.z()*ANGSTROM]);

          let residue = conformer.name().to_owned();

          let residue_sequence_number = res.serial_number() as u32;

          bath_particles.push(Particle{
            element,
            isotope: Isotope::most_common_for(&element),
            coordinates,
            serial: Some(serial as u32),
            residue: Some(residue),
            residue_sequence_number: Some(residue_sequence_number),
          });
        }
      }
    }
  }

  if n_parse_failures > 0 {
    println!("Warning, {} atoms could not be parsed.", n_parse_failures);
  }

  Ok(bath_particles)
}
//------------------------------------------------------------------------------
fn build_connections(pdb: &pdbtbx::PDB) -> Result<AdjacencyList, CluEError>
{
  let mut connections = AdjacencyList::with_capacity(pdb.atom_count());

  let mut n_connections_to_invalid_atoms = 0;
  
  let mut idx0 = 0;
  for atom0 in pdb.atoms(){
    if atom0.element() == None{
      continue;
    }
    let mut idx1 = 0;
    let n0 = atom0.serial_number();
    
    for atom1 in pdb.atoms(){
      if idx1 >= idx0 { break; }


      let are_connected: bool;
      match are_atoms_connected(atom0,atom1, pdb){
        Ok(connected_status) => are_connected = connected_status,
        Err(CluEError::AtomDoesNotSpecifyElement(serial)) => continue,
        Err(err) => return Err(err),
      }

      if are_connected{
        connections.connect(idx1,idx0);
      }

      idx1 += 1;
    }

    idx0 += 1;
  }

  Ok(connections)
}
//------------------------------------------------------------------------------
fn are_atoms_connected(atom0: &pdbtbx::Atom, atom1: &pdbtbx::Atom,
    pdb: &pdbtbx::PDB) -> Result<bool,CluEError> 
{
  
  let r: f64;
  if let Some(unit_cell) = &pdb.unit_cell{
    r = atom0.distance_wrapping(atom1,unit_cell);
  }else{
    r = atom0.distance(atom1);
  };

  let Some(elmt0) = atom0.element() else{
    let serial = atom0.serial_number();
    return Err(CluEError::AtomDoesNotSpecifyElement(serial))
  };
  let r0 = elmt0.atomic_radius().covalent_single;

  let Some(elmt1) = atom1.element() else{
    let serial = atom1.serial_number();
    return Err(CluEError::AtomDoesNotSpecifyElement(serial))
  };
  let r1 = elmt1.atomic_radius().covalent_single;

  const TOLERANCE: f64 = 0.1; 
  let are_connected = r <= (r0+r1)*(1.0 + TOLERANCE);

  Ok(are_connected)
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests {
  use super::*;
  #[test]
  fn test_parse_pdb(){
    let filename = "./assets/TEMPO.pdb";
    //let file = std::fs::read_to_string(filename).unwrap();
    let structures = parse_pdb(&filename).unwrap();

    assert_eq!(structures.len(),1);
    assert_eq!(structures[0].bath_particles.len(),29);
    assert_eq!(structures[0].bath_particles[27].element,Element::Nitrogen);
    assert_eq!(structures[0].bath_particles[28].element,Element::Oxygen);
    assert!(structures[0].connections.are_connected(27,28));

    assert_eq!(structures[0].bath_particles[27].serial, Some(28));
    assert_eq!(structures[0].bath_particles[27].residue,
        Some("TEM".to_string()));
    assert_eq!(structures[0].bath_particles[27].residue_sequence_number, 
        Some(1));
    assert_eq!(structures[0].bath_particles[27].coordinates.x(), 
        36.440*ANGSTROM);
    assert_eq!(structures[0].bath_particles[27].coordinates.y(), 
        36.900*ANGSTROM);
    assert_eq!(structures[0].bath_particles[27].coordinates.z(), 
        37.100*ANGSTROM);

    assert_eq!(structures[0].cell_offsets[0].x(),72.5676*ANGSTROM);
    assert_eq!(structures[0].cell_offsets[0].y(),0.0);
    assert_eq!(structures[0].cell_offsets[0].z(),0.0);

    assert!((structures[0].cell_offsets[1].x()-0.0).abs() < 1e12);
    assert_eq!(structures[0].cell_offsets[1].y(),72.5676*ANGSTROM);
    assert!((structures[0].cell_offsets[1].z()-0.0).abs() < 1e12);

    assert!((structures[0].cell_offsets[2].x()-0.0).abs() < 1e12);
    assert!((structures[0].cell_offsets[2].y()-0.0).abs() < 1e12);
    assert_eq!(structures[0].cell_offsets[2].z(),72.5676*ANGSTROM);


    //         O
    //   3HC   N  CH3
    //  3HC-C    C-CH3
    //    2HC    CH2 
    //        CH2

    assert!(structures[0].connections.are_connected(0,1)); 
    assert!(structures[0].connections.are_connected(0,5));
    assert!(structures[0].connections.are_connected(0,9));
    assert!(structures[0].connections.are_connected(0,27));
    assert!(structures[0].connections.are_connected(1,2));
    assert!(structures[0].connections.are_connected(1,3));
    assert!(structures[0].connections.are_connected(1,4));
    assert!(structures[0].connections.are_connected(5,6));
    assert!(structures[0].connections.are_connected(5,7));
    assert!(structures[0].connections.are_connected(5,8));
    assert!(structures[0].connections.are_connected(9,10));
    assert!(structures[0].connections.are_connected(9,11));
    assert!(structures[0].connections.are_connected(9,12));
    assert!(structures[0].connections.are_connected(12,13));
    assert!(structures[0].connections.are_connected(12,14));
    assert!(structures[0].connections.are_connected(12,15));
    assert!(structures[0].connections.are_connected(15,16));
    assert!(structures[0].connections.are_connected(15,17));
    assert!(structures[0].connections.are_connected(15,18));
    assert!(structures[0].connections.are_connected(18,19));
    assert!(structures[0].connections.are_connected(18,23));
    assert!(structures[0].connections.are_connected(18,27));
    assert!(structures[0].connections.are_connected(19,20));
    assert!(structures[0].connections.are_connected(19,21));
    assert!(structures[0].connections.are_connected(19,22));
    assert!(structures[0].connections.are_connected(23,24));
    assert!(structures[0].connections.are_connected(23,25));
    assert!(structures[0].connections.are_connected(23,26));
    assert!(structures[0].connections.are_connected(27,28));

    assert!(!structures[0].connections.are_connected(0,28));
  }

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
