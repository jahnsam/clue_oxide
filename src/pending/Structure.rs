pub struct Structure{ 
  number: usize,
  particles: Vec::<Particle>,
  pdb_indices: Vec::<usize>,
  pdb: PDB,
  number_unit_cells: usize,
}

pub struct Structure{
  number: usize,
  coordinates: mox::Mat,
  elements: Vec::<Element>,
  exchange_group_ids: Vec::<usize>,
  number_exchange_groups: usize,
  pdb_indices: Vec::<usize>,
  pdb: PDB,
  number_unit_cells: usize,
  particles: Vec::<Particle>,
}

impl Structure{
  pub fn with_capacity(n: usize) -> Structure {
  }

  //----------------------------------------------------------------------------
  pub fn push_element(&mut self, element: Element, idx_pdb: usize, pdb: &PDB )
  -> Result<(),CLUE_ERROR> {

    let n = structure.number;

    if pdb.x.len() <= n || pdb.y.len() <= n || pdb.z.len() <= n {
      return CLUE_ERROR;
    }
    let r = Mat::from(3,1,vec![pdb.x[ipdb], pdb.y, pdb.z]);

    if coordinates.n_cols() <= n || coordinates.n_rows() != 3 {
      return CLUE_ERROR;
    }

    self.coordinates.set_col(n, r);

    self.number += 1;

    Ok(())
  }

  //----------------------------------------------------------------------------
  pub fn trim(){

  }

  //----------------------------------------------------------------------------
  pub fn copy_particle(base_idx: usize, r_shift: &Mat){
  }

  //----------------------------------------------------------------------------
  pub fn set_tunnel_splittings(&mut self, tunnel_splittings: Vec::<f64>)
  -> Result<(),CLUE_ERROR> {
  }
  //----------------------------------------------------------------------------
  pub fn set_number_exchange_groups(n: usize)
  -> Result<(),CLUE_ERROR> {
  }
  //----------------------------------------------------------------------------
  pub fn set_exchange_groups(exchange_groups: Vec::<usize>)
  -> Result<(),CLUE_ERROR> {
  }
  //----------------------------------------------------------------------------
  pub fn set_exchange_group_ids(exchange_groups: Vec::<usize>,
      kept_atoms: Vec::<bool>) -> Result<(),CLUE_ERROR> {
  }
  //----------------------------------------------------------------------------

}
