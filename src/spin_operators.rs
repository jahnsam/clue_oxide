use super::physical_constants::*;

use ndarray::Array2;
use ndarray_linalg;
use num_complex::Complex;

type CxMat = Array2::<Complex<f64>>;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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

fn commutator(mat0: &CxMat, mat1: &CxMat) -> CxMat{
  mat0.dot(mat1) - mat1.dot(mat0)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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



