use std::ops::{Add,Sub,Mul,Div};
use num_complex::Complex;

pub struct Signal{
  data: Vec::<Complex<f64>>,
}
impl Signal{
  fn len(&self) -> usize{
    self.data.len()
  }
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
