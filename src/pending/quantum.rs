//extern crate matrix_oxide;
use matrix_oxide as mox;
use super::phys::*;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn get_propagators(
    hamiltonian: &mox::CxMat, 
    times: &Vec<mox::Z64>) -> Vec<mox::CxMat>{
  
  let mut propagators = Vec::with_capacity(times.len());

  let dim = hamiltonian.n_cols();
  let (eval, evec) = hamiltonian.eigh()
    .unwrap_or_else(|error| {
        hamiltonian.print("Could not diagonalize Hamiltonian:");
        panic!("Error: {:?}",error);
        });

  let evec_t = evec.trans();

  for t in times {
    let mut u = mox::CxMat::zeros(dim,dim);

    for ii in 0..dim {
      let value = (-I*2.0*PI/HBAR*t*eval.uget(ii,0) ).exp();
      u.set(ii,ii,value);
    }
    
    u = &(&evec * &u) * &evec_t;
    propagators.push(u);
  }

  propagators
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <m'|S-|m> = delta_{m'+1,m} * sqrt(S*(S+1) - m'*m).
pub fn spin_minus(spin_multiplicity: usize) -> mox::CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = mox::CxMat::zeros(spin_multiplicity,spin_multiplicity);

  for ii in 0..spin_multiplicity {
    let n = ii as f64;
    let ms = spin - n;
    op.set(ii+1,ii,
      ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt()
    );
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <m'|S+|m> = delta_{m',m+1} * sqrt(S*(S+1) - m'*m).
pub fn spin_plus(spin_multiplicity: usize) -> mox::CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = mox::CxMat::zeros(spin_multiplicity,spin_multiplicity);

  for ii in 0..spin_multiplicity {
    let n = ii as f64;
    let ms = spin - n;
    op.set(ii,ii+1,
      ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt()
    );
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <m'|S^2|m> = delta_{m',m}*S*(S+1).
pub fn spin_squared(spin_multiplicity: usize) -> mox::CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = mox::CxMat::zeros(spin_multiplicity,spin_multiplicity);

  for ii in 0..spin_multiplicity {
    let value = ONE*spin*(spin+1.0);
    op.set(ii,ii, value);
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <m'|Sx|m> = 1/2*(delta_{m',m+1} + delta_{m'+1,1})*sqrt(S*(S+1) - m'*m).
pub fn spin_x(spin_multiplicity: usize) -> mox::CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = mox::CxMat::zeros(spin_multiplicity,spin_multiplicity);

  for ii in 0..spin_multiplicity {
    let n = ii as f64;
    let ms = spin - n;
    let value = 0.5*ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
    op.set(ii,ii+1, value);
    op.set(ii+1,ii, value);
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <m'|Sy|m> = -i/2*(delta_{m',m+1} - delta_{m'+1,1})*sqrt(S*(S+1) - m'*m).
pub fn spin_y(spin_multiplicity: usize) -> mox::CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = mox::CxMat::zeros(spin_multiplicity,spin_multiplicity);

  for ii in 0..spin_multiplicity {
    let n = ii as f64;
    let ms = spin - n;
    let value = -0.5*I*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
    op.set(ii,ii+1, value);
    op.set(ii+1,ii, -value);
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <m'|Sz|m> = delta_{m',m}*m.
pub fn spin_z(spin_multiplicity: usize) -> mox::CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = mox::CxMat::zeros(spin_multiplicity,spin_multiplicity);

  for ii in 0..spin_multiplicity {
    let n = ii as f64;
    let ms = spin - n;
    op.set(ii,ii, ms*ONE);
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Copy, Clone)]
pub enum SpinOp{
  E,
  Sx,
  Sy,
  Sz,
  Sp,
  Sm,
  S2,
}



fn get_sop(spin_multiplicity: usize, sop: &SpinOp) -> mox::CxMat{
  match sop {
    SpinOp::E => mox::CxMat::eye(spin_multiplicity,spin_multiplicity),
    SpinOp::Sx => spin_x(spin_multiplicity),
    SpinOp::Sy => spin_y(spin_multiplicity),
    SpinOp::Sz => spin_z(spin_multiplicity),
    SpinOp::Sp => spin_plus(spin_multiplicity),
    SpinOp::Sm => spin_minus(spin_multiplicity),
    SpinOp::S2 => spin_squared(spin_multiplicity),
  }
}

pub fn kron_spin_op(spin_mults: &Vec<usize>, sops: &Vec<SpinOp>) -> 
  Result<mox::CxMat,String> {

  if spin_mults.len() !=  sops.len() {
    let err_str = format!("{} {} {} {} {}",
    "unequal dimensions: spin_multiplicitis have",
    spin_mults.len(), 
    "elements, but spin operators have",
    sops.len(),
    "elements.",
    );
    return Err(err_str);
  }

  let mut sop = mox::CxMat::eye(1,1); 
  for ii in 0..spin_mults.len() {
    let s = get_sop(spin_mults[ii],&sops[ii]);
    sop = sop.kron(&s);
  }

  Ok(sop)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct KronSpinOpList {
  sop_list: Vec::<mox::CxMat>,
}

impl<'a> KronSpinOpList {

  pub fn new(
      spin_multiplicity: usize, 
      spin_operator: SpinOp, 
      max_size: usize) -> Result<KronSpinOpList,String> {
  
    let n_ops = (max_size*( max_size+1) )/2;
    let mut sop_list = Vec::<mox::CxMat>::with_capacity(n_ops);

    for n_ops in 1..=max_size {
      let spin_mults = vec![spin_multiplicity; n_ops];

      for p_idx in 0..n_ops {

        let mut sops = vec![SpinOp::E; n_ops];
        sops[p_idx] = spin_operator;

        let sop = kron_spin_op(&spin_mults,&sops)?;
        sop_list.push(sop);
      }
    }

    Ok(KronSpinOpList {
      sop_list,
    })
  }

  pub fn get(&'a self, op_pos: usize, n_ops: usize) -> 
    Result<&'a mox::CxMat,String> {


    if op_pos >= n_ops {
      let err_str = format!("operator position = {} >= {} = n_positions",
          op_pos, n_ops);
      return Err(err_str);
    }

    let idx = KronSpinOpList::get_index(op_pos, n_ops);

    if idx >= self.sop_list.len() as i32 || idx < 0{
      let err_str = String::from("spin operator is unavailable");
      return Err(err_str);
    }

    Ok(&self.sop_list[idx as usize])
  }

  fn get_index(op_pos: usize, n_ops: usize) -> i32 {
    let op_pos = op_pos as i32;
    let n_ops = n_ops as i32;
    ( n_ops*(n_ops - 1) )/2 + op_pos
  }
}


pub struct KronSpinOpXYZ {
  sx_list: KronSpinOpList,
  sy_list: KronSpinOpList,
  sz_list: KronSpinOpList,
}

impl<'a> KronSpinOpXYZ {

  pub fn new(
      spin_multiplicity: usize, 
      max_size: usize) -> Result<KronSpinOpXYZ,String> {

    let sx_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sx, max_size)?;
    let sy_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sy, max_size)?;
    let sz_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sz, max_size)?;

    Ok(KronSpinOpXYZ{
      sx_list,  
      sy_list,  
      sz_list,  
    })
  }

  pub fn get(&'a self, sop: &SpinOp, op_pos: usize, n_ops: usize) -> 
   Result<&'a mox::CxMat,String> {

     match sop {
       SpinOp::Sx => self.sx_list.get(op_pos,n_ops),
       SpinOp::Sy => self.sy_list.get(op_pos,n_ops),
       SpinOp::Sz => self.sz_list.get(op_pos,n_ops),
       _ => {
         let err_str = String::from("spin operator unavailable");
         return Err(err_str);
       }
     }
  } 
}

pub struct ClusterSpinOperators {
  max_size: usize, 
  spin_multiplicities: Vec<usize>, 
  cluster_spin_ops: Vec<KronSpinOpXYZ>,
}

impl<'a> ClusterSpinOperators {
  pub fn new(
      spin_multiplicities: &Vec<usize>, max_size: usize) 
   -> Result<ClusterSpinOperators, String> {
    
    let n_mults = spin_multiplicities.len();

    let mut cluster_spin_ops = Vec::<KronSpinOpXYZ>::with_capacity(n_mults);

    for spin_multiplicity in spin_multiplicities.iter() {
      let sops = KronSpinOpXYZ::new(*spin_multiplicity, max_size)?;
      cluster_spin_ops.push(sops);
    }

    Ok(ClusterSpinOperators{
        max_size: max_size,
        spin_multiplicities: spin_multiplicities.clone(),
        cluster_spin_ops: cluster_spin_ops,
        })
  }

  pub fn get(&'a self, 
      sop: &SpinOp, spin_multiplicity: usize, cluster_size: usize,
      op_pos: usize) -> Result<&'a mox::CxMat,String> {
  
    if cluster_size > self.max_size {
      let err_str = format!("requested size {}, exceeds max size {} ", 
          cluster_size,self.max_size);
      return Err(err_str);
    }

    for (ii,ispin_mult) in self.spin_multiplicities.iter().enumerate() {
      let ispin_mult = *ispin_mult;
      if ispin_mult == spin_multiplicity {
        
       let sop = self.cluster_spin_ops[ii]
         .get(sop, op_pos,cluster_size)?;
       return Ok(sop);
      }
    }
    let err_str = format!(
     "could not find spin operator with multiplicity {} for a {}-cluster.",
     spin_multiplicity, cluster_size);
    Err(err_str)
  }


}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#[cfg(test)]
mod tests {
  use super::*;
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  fn test_get_propagators() {
    let sy = spin_y(2);

    let times = vec![ONE/4.0*HBAR];

    let propagators = get_propagators(&sy,&times);
    let u = &propagators[0];

    for irow in  0..2 {
      for icol in 1..2 {
        let err: f64 = mox::Z64::norm( 
            u.uget(irow,icol).conj()*u.uget(irow,icol) - 0.5
            );
        assert!(err < ERROR_THRESHOLD);
      }
    }

  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  #[allow(non_snake_case)]
  fn test_ClusterSpinOperators() {

    let spin_multiplicities = vec![2,3];
    let max_size = 3;
    let sops = ClusterSpinOperators::new(&spin_multiplicities,max_size).unwrap();


    for ispin_mult in &spin_multiplicities {
      let ispin_mult = *ispin_mult;
      for cluster_size in 1..=max_size{
        for op_pos in 0..cluster_size {
          let sx = sops.get(&SpinOp::Sx,ispin_mult,cluster_size,op_pos)
            .unwrap();

          let sy = sops.get(&SpinOp::Sy,ispin_mult,cluster_size,op_pos)
            .unwrap();

          let sz = sops.get(&SpinOp::Sz,ispin_mult,cluster_size,op_pos)
            .unwrap();

          let sp = sx + &sy.scalar_multiply(I);
          let sm = sx - &sy.scalar_multiply(I);
          let s2 = &(&(sx*sx) + &(sy*sy)) + &(sz*sz);

          assert_eq!(sx.n_cols(), ispin_mult.pow(cluster_size as u32));
          assert!(check_spin_ops(sx,sy,sz,&sp,&sm,&s2));
        }
      }
    }

  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  #[allow(non_snake_case)]
  fn test_KronSpinOpList() {

    let spin_multiplicity = 2;
    let max_size = 2;

    let sx_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sx, max_size).unwrap();

    let sy_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sy, max_size).unwrap();

    let sz_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sz, max_size).unwrap();

    let sp_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sp, max_size).unwrap();

    let sm_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sm, max_size).unwrap();

    let s2_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::S2, max_size).unwrap();

    for n_ops in 1..=max_size{
      for op_pos in 0..n_ops {
        let sx = sx_list.get(op_pos,n_ops).unwrap();
        let sy = sy_list.get(op_pos,n_ops).unwrap();
        let sz = sz_list.get(op_pos,n_ops).unwrap();
        let sp = sp_list.get(op_pos,n_ops).unwrap();
        let sm = sm_list.get(op_pos,n_ops).unwrap();
        let s2 = s2_list.get(op_pos,n_ops).unwrap();

        assert_eq!(sx.n_cols(), spin_multiplicity.pow(n_ops as u32));
        assert!(check_spin_ops(sx,sy,sz,sp,sm,s2));
      }
  }

  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  #[allow(non_snake_case)]
  fn test_KronSpinOpXYZ() {
    let spin_multiplicity = 2;
    let max_size = 2;
    let sops = KronSpinOpXYZ::new(spin_multiplicity, max_size).unwrap();

    for n_ops in 1..=max_size{
      for op_pos in 0..n_ops {
        let sx = sops.get(&SpinOp::Sx, op_pos,n_ops).unwrap(); 
        let sy = sops.get(&SpinOp::Sy, op_pos,n_ops).unwrap(); 
        let sz = sops.get(&SpinOp::Sz, op_pos,n_ops).unwrap(); 
        let sp = sx + &sy.scalar_multiply(I);
        let sm = sx - &sy.scalar_multiply(I);
        let s2 = &(&(sx*sx) + &(sy*sy)) + &(sz*sz);

        assert_eq!(sx.n_cols(), spin_multiplicity.pow(n_ops as u32));
        assert!(check_spin_ops(sx,sy,sz,&sp,&sm,&s2));
      }
    }

  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  fn test_spin_ops() {

    for spin_multiplicity in 0..7 {
      let sx = spin_x(spin_multiplicity);
      let sy = spin_y(spin_multiplicity);
      let sz = spin_z(spin_multiplicity);
      let sp = spin_plus(spin_multiplicity);
      let sm = spin_minus(spin_multiplicity);
      let s2 = spin_squared(spin_multiplicity);


      assert_eq!(sx.n_cols(), spin_multiplicity);

      assert!(check_spin_ops(&sx,&sy,&sz,&sp,&sm,&s2));
      

    }
  }

  fn check_spin_ops(
      sx: &mox::CxMat,
      sy: &mox::CxMat,
      sz: &mox::CxMat,
      sp: &mox::CxMat,
      sm: &mox::CxMat,
      s2: &mox::CxMat,
      ) -> bool {
  
      let sx2 = sx * sx;
      let sy2 = sy * sy;
      let sz2 = sz * sz;
    let spin2 = &(&sx2 + &sy2) + &sz2;

    let mut pass: bool = true;  

    pass &= mox::CxMat::approx_eq( 
          &mox::CxMat::commutator(sx,sy).unwrap(), 
          &sz.scalar_multiply(I) ) ; 
      
    pass &= mox::CxMat::approx_eq( 
          &mox::CxMat::commutator(sz,sx).unwrap(), 
          &sy.scalar_multiply(I) ) ; 

    pass &= mox::CxMat::approx_eq( 
          &mox::CxMat::commutator(sy,sz).unwrap(), 
          &sx.scalar_multiply(I) ) ; 

    pass &= mox::CxMat::approx_eq( 
          sp, 
          &(sx + &sy.scalar_multiply(I) ) 
          ); 

    pass &= mox::CxMat::approx_eq( 
          sm, 
          &(sx - &sy.scalar_multiply(I) ) 
          );
    pass &= mox::CxMat::approx_eq( 
            s2, 
            &spin2
            ); 
    pass 
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #[test]
  fn test_spin_op() {

    let spin_mults = Vec::<usize>::from([2,2]);

    let xx = kron_spin_op(&spin_mults,&vec![SpinOp::Sx,SpinOp::Sx]).unwrap();
    let yy = kron_spin_op(&spin_mults,&vec![SpinOp::Sy,SpinOp::Sy]).unwrap();
    let zz = kron_spin_op(&spin_mults,&vec![SpinOp::Sz,SpinOp::Sz]).unwrap();
    let pm = kron_spin_op(&spin_mults,&vec![SpinOp::Sp,SpinOp::Sm]).unwrap();
    let mp = kron_spin_op(&spin_mults,&vec![SpinOp::Sm,SpinOp::Sp]).unwrap();

    assert!(mox::CxMat::approx_eq( 
          &( &(&xx + &yy) + &zz), 
          &(&zz + &(&pm.scalar_multiply(0.5*ONE)+&mp.scalar_multiply(0.5*ONE)))
          )); 
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}
