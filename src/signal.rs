use std::ops::{Add,Sub,Mul,Div};
use num_complex::Complex;

pub struct Signal{
  data: Vec::<Complex<f64>>,
}
impl Signal{
  pub fn len(&self) -> usize{
    self.data.len()
  }
  //----------------------------------------------------------------------------
  pub fn ones(n: usize) -> Self {
    let mut data = Vec::<Complex<f64>>::with_capacity(n);

    for ii in 0..n{
      data.push(Complex::<f64>{re:1.0,im: 0.0});
    }
    Signal{data}
  }
  //----------------------------------------------------------------------------
}

/*
impl Add for Signal{
  fn add(&self,other: &Signal) -> Signal{
    
  }
}
impl Sub for Signal{}
impl Mul for Signal{}
impl Div for Signal{}
*/
