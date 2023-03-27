pub mod calculate_cluster_signal;
pub mod calculate_analytic_restricted_2cluster_signals;
pub mod calculate_signal;
use std::ops::{Add,Sub,Mul,Div};
use num_complex::Complex;

#[derive(PartialEq,Debug,Clone,Default)]
pub struct Signal{
  pub data: Vec::<Complex<f64>>,
}
impl Signal{
  //----------------------------------------------------------------------------
  pub fn len(&self) -> usize{
    self.data.len()
  }
  //----------------------------------------------------------------------------
  pub fn new() -> Self {Signal::default()}
  //----------------------------------------------------------------------------
  pub fn ones(n: usize) -> Self {
    let mut data = Vec::<Complex<f64>>::with_capacity(n);

    for _ii in 0..n{
      data.push(Complex::<f64>{re:1.0,im: 0.0});
    }
    Signal{data}
  }
  //----------------------------------------------------------------------------
}

impl Add for Signal{
  type Output = Self;
 
  fn add(self,other: Signal) -> Self{
    assert_eq!(self.len(), other.len());
    Signal{
      data: self.data.iter().zip(other.data).map(|(x,y)| x+y).collect(),
    }
  }
}

impl Sub for Signal{
  type Output = Self;
 
  fn sub(self,other: Signal) -> Self{
    assert_eq!(self.len(), other.len());
    Signal{
      data: self.data.iter().zip(other.data).map(|(x,y)| x-y).collect(),
    }
  }
}


impl Mul for Signal{
  type Output = Self;
 
  fn mul(self,other: Signal) -> Self{
    assert_eq!(self.len(), other.len());
    Signal{
      data: self.data.iter().zip(other.data).map(|(x,y)| x*y).collect(),
    }
  }
}
impl Div for Signal{
  type Output = Self;
 
  fn div(self,other: Signal) -> Self{
    assert_eq!(self.len(), other.len());
    Signal{
      data: self.data.iter().zip(other.data).map(|(x,y)| x/y).collect(),
    }
  }
}

#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_signal_ops(){
    const ONE: Complex<f64> = Complex::<f64>{re:1.0,im:0.0};
    let signal0 = Signal{ data: vec![ONE,ONE],};
    let signal1 = Signal{ data: vec![ONE,0.5*ONE],};

    let signal2 = signal0.clone() + signal1.clone();
    assert_eq!(signal2.data,vec![2.0*ONE,1.5*ONE]);

    let signal3 = signal0.clone() - signal1.clone();
    assert_eq!(signal3.data,vec![0.0*ONE,0.5*ONE]);

    let signal4 = signal0.clone() * signal1.clone();
    assert_eq!(signal4.data,vec![ONE,0.5*ONE]);

    let signal5 = signal0.clone() / signal1.clone();
    assert_eq!(signal5.data,vec![ONE,2.0*ONE]);
  }
}


