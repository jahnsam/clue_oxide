use crate::physical_constants::*;
use crate::clue_errors::*;
use crate::config::Config;
use crate::HamiltonianTensors;

use std::fmt;
use ndarray::Array2;
use ndarray_linalg::{Eigh, UPLO};
use ndarray::linalg::kron;
use num_complex::Complex;

type CxMat = Array2::<Complex<f64>>;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct SpinHamiltonian{
  beta: CxMat,
  alpha: CxMat,
}
pub fn build_hamiltonian(spin_indices: &Vec::<usize>,
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config)
  -> Result<SpinHamiltonian,CluEError>
{


  let Some(central_spin) = config.detected_spin_identity else{
    return Err(CluEError::NoCentralSpin);
  };

  let central_spin_mult = central_spin.spin_multiplicity();
  let s = (central_spin_mult as f64 - 1.0)/2.0;
  let spin_ms: Vec::<f64> = (0..central_spin_mult).map(|n| (n as f64 - s))
    .collect();

  let Some(transition) = config.detected_spin_transition else {
    return Err(CluEError::NoCentralSpinTransition);
  };

  let ms_beta = spin_ms[transition[0]];
  let ms_alpha = spin_ms[transition[1]];


  let spin_multiplicities: Vec::<usize> = 
  spin_indices.iter().map(|idx| tensors.spin_multiplicities[*idx]).collect();

  let mut dim: usize = 1;
  spin_multiplicities.iter().for_each(|spin_mul| dim *= spin_mul);

  let mut ham0 = CxMat::zeros((dim,dim));
  let mut ham_ms = CxMat::zeros((dim,dim));
  
  // electron Zeeman
  if let Some(vec) = tensors.spin1_tensors.get(0){
    ham_ms = CxMat::eye(dim)*vec.z();
  }

  let cluster_size = spin_indices.len();


  for (sop_idx0, &ten_idx0) in spin_indices.iter().enumerate(){

    let spin_mult0 = tensors.spin_multiplicities[ten_idx0];
    let sx0 = spin_ops.get(&SpinOp::Sx,spin_mult0,cluster_size,sop_idx0)?;
    let sy0 = spin_ops.get(&SpinOp::Sy,spin_mult0,cluster_size,sop_idx0)?;
    let sz0 = spin_ops.get(&SpinOp::Sz,spin_mult0,cluster_size,sop_idx0)?;

    // nuclear Zeeman
    if let Some(vec) = tensors.spin1_tensors.get(ten_idx0){
      ham0 = ham0 + sx0*vec.x();
      ham0 = ham0 + sy0*vec.y();
      ham0 = ham0 + sz0*vec.z();
    }

    // nuclear hyperfine
    if let Some(ten) = tensors.spin2_tensors.get(0,ten_idx0){
      ham_ms = ham_ms + sx0*ten.xx();
      ham_ms = ham_ms + sx0*ten.xy();
      ham_ms = ham_ms + sx0*ten.xz();
      ham_ms = ham_ms + sy0*ten.yx();
      ham_ms = ham_ms + sy0*ten.yy();
      ham_ms = ham_ms + sy0*ten.yz();
      ham_ms = ham_ms + sz0*ten.zx();
      ham_ms = ham_ms + sz0*ten.zy();
      ham_ms = ham_ms + sz0*ten.zz();
    }

    for (sop_idx1, &ten_idx1) in spin_indices.iter().enumerate().skip(sop_idx0){

      let spin_mult1 = tensors.spin_multiplicities[ten_idx1];
      let sx1 = spin_ops.get(&SpinOp::Sx,spin_mult1,cluster_size,sop_idx1)?;
      let sy1 = spin_ops.get(&SpinOp::Sy,spin_mult1,cluster_size,sop_idx1)?;
      let sz1 = spin_ops.get(&SpinOp::Sz,spin_mult1,cluster_size,sop_idx1)?;
      
      // dipole-dipole, and electric quadrupole
      if let Some(ten) = tensors.spin2_tensors.get(ten_idx0,ten_idx1){
        ham0 = ham0 + sx0.dot(sx1)*ten.xx();
        ham0 = ham0 + sx0.dot(sy1)*ten.xy();
        ham0 = ham0 + sx0.dot(sz1)*ten.xz();
        ham0 = ham0 + sy0.dot(sx1)*ten.yx();
        ham0 = ham0 + sy0.dot(sy1)*ten.yy();
        ham0 = ham0 + sy0.dot(sz1)*ten.yz();
        ham0 = ham0 + sz0.dot(sx1)*ten.zx();
        ham0 = ham0 + sz0.dot(sy1)*ten.zy();
        ham0 = ham0 + sz0.dot(sz1)*ten.zz();
      }
     
      
    }
  }

  let beta = ham0.clone() + ham_ms.clone()*ms_beta;
  let alpha = ham0 + ham_ms*ms_alpha;
  
  Ok(SpinHamiltonian{beta,alpha})
}
//------------------------------------------------------------------------------
pub fn get_propagators(hamiltonian: &CxMat, times: &Vec::<f64>  )
  -> Result<Vec::<CxMat>,CluEError> 
{
  let Ok((eigvals, eigvecs)) = hamiltonian.eigh(UPLO::Lower) else{
    return Err(
        CluEError::CannotDiagonalizeHamiltonian(hamiltonian.to_string()));
  };

  let mut propagators = Vec::<CxMat>::with_capacity(times.len());

  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  for &t in times.iter(){
    let u_eig = CxMat::from_diag(&eigvals.map(|nu| 
          { let i_phase: Complex<f64> = (-I*2.0*PI*nu)*t;
            i_phase.exp() 
          }  
          )
        );

    let u = eigvecs.dot( &u_eig.dot( &inv_eigvecs) );

    propagators.push(u);

  }

  Ok(propagators)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct ClusterSpinOperators {
  max_size: usize, 
  spin_multiplicities: Vec<usize>, 
  cluster_spin_ops: Vec<KronSpinOpXYZ>,
}
//------------------------------------------------------------------------------
impl<'a> ClusterSpinOperators {
  pub fn new(spin_multiplicities: &Vec<usize>, max_size: usize) 
   -> Result<ClusterSpinOperators, CluEError> {
    
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
      op_pos: usize) -> Result<&'a CxMat,CluEError> {
  
    if cluster_size > self.max_size {
      return Err(CluEError::NoSpinOpForClusterSize(cluster_size,self.max_size));
    }

    for (ii,&ispin_mult) in self.spin_multiplicities.iter().enumerate() {
      if ispin_mult == spin_multiplicity {
        
       let sop = self.cluster_spin_ops[ii]
         .get(sop, op_pos,cluster_size)?;
       return Ok(sop);
      }
    }
    Err(CluEError::NoSpinOpWithMultiplicity(spin_multiplicity))
  }


}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct KronSpinOpXYZ {
  sx_list: KronSpinOpList,
  sy_list: KronSpinOpList,
  sz_list: KronSpinOpList,
}

impl<'a> KronSpinOpXYZ {

  pub fn new(spin_multiplicity: usize, max_size: usize) 
    -> Result<KronSpinOpXYZ,CluEError> 
  {

    let sx_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sx, max_size)?;
    let sy_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sy, max_size)?;
    let sz_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sz, max_size)?;

    Ok(KronSpinOpXYZ{
      sx_list,  
      sy_list,  
      sz_list,  
    })
  }
  //----------------------------------------------------------------------------
  pub fn get(&'a self, sop: &SpinOp, op_pos: usize, n_ops: usize) 
    -> Result<&'a CxMat,CluEError> 
  {

     match sop {
       SpinOp::Sx => self.sx_list.get(op_pos,n_ops),
       SpinOp::Sy => self.sy_list.get(op_pos,n_ops),
       SpinOp::Sz => self.sz_list.get(op_pos,n_ops),
       _ => Err(CluEError::CannotFindSpinOp(sop.to_string())),
     }
  } 
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct KronSpinOpList {
  sop_list: Vec::<CxMat>,
}
//------------------------------------------------------------------------------
impl<'a> KronSpinOpList {

  pub fn new(
      spin_multiplicity: usize,
      spin_operator: SpinOp,
      max_size: usize) -> Result<KronSpinOpList,CluEError> {

    let n_ops = (max_size*( max_size+1) )/2;
    let mut sop_list = Vec::<CxMat>::with_capacity(n_ops);

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
  //----------------------------------------------------------------------------
  pub fn get(&'a self, op_pos: usize, n_ops: usize) ->
    Result<&'a CxMat,CluEError> {


    if op_pos >= n_ops {
      return Err(CluEError::UnavailableSpinOp(op_pos,n_ops));
    }

    let idx = KronSpinOpList::get_index(op_pos, n_ops);

    if idx >= self.sop_list.len(){
      return Err(CluEError::UnavailableSpinOp(idx,self.sop_list.len()));
    }

    Ok(&self.sop_list[idx])
  }
//------------------------------------------------------------------------------
  fn get_index(op_pos: usize, n_ops: usize) -> usize {
    ( n_ops*(n_ops - 1) )/2 + op_pos
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//------------------------------------------------------------------------------
pub fn kron_spin_op(spin_mults: &Vec<usize>, sops: &Vec<SpinOp>) ->
  Result<CxMat,CluEError> 
{

  if spin_mults.len() !=  sops.len() {
    return Err(CluEError::UnequalLengths(
          "spin_multiplicities".to_string(),
          spin_mults.len(),
          "spin operators".to_string(),
          sops.len(),
    ));
  }

  let mut sop = CxMat::eye(1);
  for ii in 0..spin_mults.len() {
    let s = get_sop(spin_mults[ii],&sops[ii]);
    sop = kron(&sop,&s);
  }

  Ok(sop)
}
//------------------------------------------------------------------------------
#[derive(Copy,Debug,Clone,PartialEq)]
pub enum SpinOp{
  E,
  Sx,
  Sy,
  Sz,
  Sp,
  Sm,
  S2,
}
impl fmt::Display for SpinOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      match self{
        SpinOp::E => write!(f, "E"),
        SpinOp::Sx => write!(f, "Sx"),
        SpinOp::Sy => write!(f, "Sy"),
        SpinOp::Sz => write!(f, "Sz"),
        SpinOp::Sp => write!(f, "S+"),
        SpinOp::Sm => write!(f, "S-"),
        SpinOp::S2 => write!(f, "S^2"),
      }
    }
}
//------------------------------------------------------------------------------
fn get_sop(spin_multiplicity: usize, sop: &SpinOp) -> CxMat{
  match sop {
    SpinOp::E => CxMat::eye(spin_multiplicity),
    SpinOp::Sx => spin_x(spin_multiplicity),
    SpinOp::Sy => spin_y(spin_multiplicity),
    SpinOp::Sz => spin_z(spin_multiplicity),
    SpinOp::Sp => spin_plus(spin_multiplicity),
    SpinOp::Sm => spin_minus(spin_multiplicity),
    SpinOp::S2 => spin_squared(spin_multiplicity),
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//------------------------------------------------------------------------------
// <m'|Sx|m> = 1/2*(delta_{m',m+1} + delta_{m'+1,1})*sqrt(S*(S+1) - m'*m).
pub fn spin_x(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    let value = 0.5*ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
    op[[ii,ii+1]] = value;
    op[[ii+1,ii]] = value;
  }

  op
}
//------------------------------------------------------------------------------
// <m'|Sy|m> = -i/2*(delta_{m',m+1} - delta_{m'+1,1})*sqrt(S*(S+1) - m'*m).
pub fn spin_y(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    let value = -0.5*I*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
    op[[ii,ii+1]] = value;
    op[[ii+1,ii]] = -value;
  }

  op
}
//------------------------------------------------------------------------------
// <m'|Sz|m> = delta_{m',m}*m.
pub fn spin_z(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity  {
    let n = ii as f64;
    let ms = spin - n;
    op[[ii,ii]] =  ms*ONE;
  }

  op
}
//------------------------------------------------------------------------------
// <m'|S-|m> = delta_{m'+1,m} * sqrt(S*(S+1) - m'*m).
pub fn spin_minus(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    op[[ii+1,ii]] = ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
  }

  op
}
//------------------------------------------------------------------------------
// <m'|S+|m> = delta_{m',m+1} * sqrt(S*(S+1) - m'*m).
pub fn spin_plus(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    op[[ii,ii+1]] = ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
  }

  op
}
//------------------------------------------------------------------------------
// <m'|S^2|m> = delta_{m',m}*S*(S+1).
pub fn spin_squared(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity  {
    let value = ONE*spin*(spin+1.0);
    op[[ii,ii]] = value;
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#[cfg(test)]
mod tests {
  use super::*;


  //----------------------------------------------------------------------------
  #[test]
  fn test_get_propagators() {
    let sy = spin_y(2);
    let ham = sy*2.0;

    let times = vec![1.0/8.0];

    let propagators = get_propagators(&ham,&times).unwrap();
    let u = &propagators[0];

    for irow in  0..2 {
      for icol in 1..2 {
        let err: f64 = (u[[irow,icol]].conj()*u[[irow,icol]] - 0.5).norm();
        assert!(err < 1e-12);
      }
    }

  }  
  //----------------------------------------------------------------------------
  fn approx_eq(m0: &CxMat, m1: &CxMat, tol: f64) -> bool {

    assert!(tol > 0.0);

    if m0.ncols() != m1.ncols() || m0.nrows() != m1.nrows() {
      return false;
    }
  
    for irow in 0..m0.nrows() {
      for icol in 0..m0.ncols(){
        let err = ( m0[[irow,icol]] - m1[[irow,icol]] ).norm();
        if err >= tol {
          return false;
        }
      }
    }
    true
  }
  //----------------------------------------------------------------------------
  fn commutator(mat0: &CxMat, mat1: &CxMat) -> CxMat{
    mat0.dot(mat1) - mat1.dot(mat0)
  }
  //----------------------------------------------------------------------------
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

          let sp = sx + sy*I;
          let sm = sx - sy*I;
          let s2 = sx.dot(sx) + sy.dot(sy) + sz.dot(sz);

          assert_eq!(sx.ncols(), ispin_mult.pow(cluster_size as u32));
          assert!(check_spin_ops(sx,sy,sz,&sp,&sm,&s2));
        }
      }
    }

  }
  //----------------------------------------------------------------------------
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
        let sp = sx + sy*I;
        let sm = sx - sy*I;
        let s2 = sx.dot(sx) + sy.dot(sy) + sz.dot(sz);

        assert_eq!(sx.ncols(), spin_multiplicity.pow(n_ops as u32));
        assert!(check_spin_ops(sx,sy,sz,&sp,&sm,&s2));
      }
    }

  }
  //----------------------------------------------------------------------------
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

        assert_eq!(sx.ncols(), spin_multiplicity.pow(n_ops as u32));
        assert!(check_spin_ops(sx,sy,sz,sp,sm,s2));
      }
  }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kron_spin_op() {

    let spin_mults = Vec::<usize>::from([2,2]);

    let xx = kron_spin_op(&spin_mults,&vec![SpinOp::Sx,SpinOp::Sx]).unwrap();
    let yy = kron_spin_op(&spin_mults,&vec![SpinOp::Sy,SpinOp::Sy]).unwrap();
    let zz = kron_spin_op(&spin_mults,&vec![SpinOp::Sz,SpinOp::Sz]).unwrap();
    let pm = kron_spin_op(&spin_mults,&vec![SpinOp::Sp,SpinOp::Sm]).unwrap();
    let mp = kron_spin_op(&spin_mults,&vec![SpinOp::Sm,SpinOp::Sp]).unwrap();

    assert!(approx_eq(
          &( &(&xx + &yy) + &zz), 
          &(&zz + &(&(pm*(0.5*ONE))+&(mp*(0.5*ONE)))),1e-12)
        );
  }
//------------------------------------------------------------------------------
  #[test]
  fn test_spin_ops() {

    for spin_multiplicity in 0..7 {
      let sx = spin_x(spin_multiplicity);
      let sy = spin_y(spin_multiplicity);
      let sz = spin_z(spin_multiplicity);
      let sp = spin_plus(spin_multiplicity);
      let sm = spin_minus(spin_multiplicity);
      let s2 = spin_squared(spin_multiplicity);


      assert_eq!(sx.ncols(), sx.nrows());
      assert_eq!(sx.ncols(), spin_multiplicity);

      assert!(check_spin_ops(&sx,&sy,&sz,&sp,&sm,&s2));
      

    }
  }
  //----------------------------------------------------------------------------
  
  fn check_spin_ops(
      sx: &CxMat,
      sy: &CxMat,
      sz: &CxMat,
      sp: &CxMat,
      sm: &CxMat,
      s2: &CxMat,
      ) -> bool {

    let tol = 1e-12;

    let sx2 = sx.dot(sx);
    let sy2 = sy.dot(sy);
    let sz2 = sz.dot(sz);
    let spin2 = sx2 + sy2 + sz2;

    let mut pass: bool = true;  

    pass &= approx_eq( 
          &commutator(sx,sy), 
          &(I*sz), tol ) ; 
      
    pass &= approx_eq( 
          &commutator(sz,sx), 
          &(I*sy), tol ) ; 

    pass &= approx_eq( 
          &commutator(sy,sz), 
          &(I*sx), tol ) ; 
      
    pass &= approx_eq( 
          &sp, 
          &(sx + I*sy), tol ) ; 
      
    pass &= approx_eq( 
          &sm, 
          &(sx - I*sy), tol ) ; 
      
    pass &= approx_eq( 
          &s2, 
          &spin2, tol ) ; 
      
      
    pass 
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

}



