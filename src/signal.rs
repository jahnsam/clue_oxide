pub mod calculate_analytic_restricted_2cluster_signals;
pub mod calculate_signal;
pub mod cluster_correlation_expansion;

use crate::cluster::find_clusters::ClusterSet;
use std::ops::{Add,Sub,Mul,Div};
use num_complex::Complex;

use crate::CluEError;
use std::error::Error;

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
  pub fn is_empty(&self) -> bool{
    self.data.is_empty()
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
  pub fn zeros(n: usize) -> Self {
    let mut data = Vec::<Complex<f64>>::with_capacity(n);

    for _ii in 0..n{
      data.push(Complex::<f64>{re:0.0,im: 0.0});
    }
    Signal{data}
  }
  //----------------------------------------------------------------------------
  pub fn scale(&mut self, scale_factor: Complex<f64>){

    for z in self.data.iter_mut(){
      *z *= scale_factor;
    }

  }
  //----------------------------------------------------------------------------
  pub fn read_from_csv(filename: &str) -> Result<Self,CluEError>
  {
    match Self::read_single_signal_from_csv(filename){
      Ok(data) => Ok(Signal{data}),
      Err(_) => Err(CluEError::CannotOpenFile(filename.to_string())),
    }

  }
  //----------------------------------------------------------------------------
  fn read_single_signal_from_csv(filename: &str)
    -> Result<Vec::<Complex<f64>>,Box<dyn Error>>
  {
    // Count data points.
    let mut rdr = csv::Reader::from_path(filename)?;
    let mut num_data = 0;
    for _result in rdr.records() {
        num_data += 1;
    }

    // Load signal.
    let mut signal = Vec::<Complex<f64>>::with_capacity(num_data);
    let mut rdr = csv::Reader::from_path(filename)?;
    for result in rdr.records() {
        let record = result?;
        if let Some(v) = record.get(0){
          let v = v.parse::<Complex<f64>>()?;
          signal.push(v);
        }
    }
    Ok(signal)
  }
  //----------------------------------------------------------------------------
  pub fn write_to_csv(&self, filename: &str) -> Result<(),CluEError>{
    match self.write_single_signal_to_csv(filename){
      Ok(()) => Ok(()),
      Err(_) => Err(CluEError::CannotWriteFile(filename.to_string())),
    }
  }
  //----------------------------------------------------------------------------
  fn write_single_signal_to_csv(&self, filename: &str) 
    -> Result<(),Box<dyn Error>>
  {
    let mut wtr = csv::Writer::from_path(filename)?;
    wtr.write_record(["signal"])?;
    for v in self.data.iter(){
      wtr.write_record(&[v.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
  }
  //----------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
pub fn write_vec_signals(signals: &Vec::<Signal>, filename: &str) 
  -> Result<(),CluEError>
{

   let headers = get_signal_headers(signals.len());
  
   if headers.len() != signals.len(){
     return Err(CluEError::MissingHeader(filename.to_string()));
   }

   for signal in signals.iter(){
     if signal.data.len() != signals[0].len(){
       return Err(CluEError::AllSignalsNotSameLength(filename.to_string()));
     }
   }

    match write_vec_signals_to_csv(signals,filename,headers){
      Ok(()) => Ok(()),
      Err(_) => Err(CluEError::CannotWriteFile(filename.to_string())),
    }
}
//------------------------------------------------------------------------------
fn get_signal_headers(n: usize) -> Vec::<String>
{
  let mut headers = Vec::<String>::with_capacity(n);
  for ii in 0..n{
    headers.push(format!("signal_{}",ii+1));
  }

  headers
}
//------------------------------------------------------------------------------
fn write_vec_signals_to_csv(signals: &[Signal],filename: &str,
    headers: Vec::<String>) -> Result<(),Box<dyn Error>>
{
  let n_data = signals[0].data.len();

  let mut wtr = csv::Writer::from_path(filename)?;
  
  wtr.write_record(&headers)?;

  for ii in 0..n_data{
    
    let mut rec = Vec::<String>::with_capacity(n_data);
   
    for signal in signals.iter(){
      rec.push(signal.data[ii].to_string());
    }

    wtr.write_record(&rec)?;
  }

  wtr.flush()?;
  
  Ok(())
}
//------------------------------------------------------------------------------
pub fn write_batch_signals(cluster_set: &ClusterSet, 
    cluster_size: usize, n_data: usize, idx: usize, batch_size: usize, 
    filename: &str) -> Result<(),CluEError>
{

  if cluster_size > cluster_set.clusters.len(){
    return Err(CluEError::NoClustersOfSize(cluster_size));
  }


  let headers = cluster_set.clusters[cluster_size - 1].iter()
    .skip(idx).take(batch_size).map(|cluster| cluster.to_header())
    .collect::<Vec::<String>>();

  let Ok(mut wtr) = csv::Writer::from_path(filename) else{
    return Err(CluEError::CannotWriteFile(filename.to_string()));
  };

  if wtr.write_record(&headers).is_err(){
    return Err(CluEError::CannotWriteFile(filename.to_string()));
  }

  for ii in 0..n_data{
    
    let mut rec = Vec::<String>::with_capacity(n_data);

    for cluster in cluster_set.clusters[cluster_size - 1].iter()
      .skip(idx).take(batch_size){

        match &cluster.signal{
          Err(err)  => return Err(err.clone()),
          _ => (),
        }

        let Ok(Some(signal)) = &cluster.signal else{
          return Err(CluEError::ClusterHasNoSignal(cluster.to_string()));
        };

        if signal.data.len() != n_data{
          return Err(CluEError::AllSignalsNotSameLength(filename.to_string()));
        }

        rec.push(signal.data[ii].to_string());
    }

    if wtr.write_record(&rec).is_err() {
      return Err(CluEError::CannotWriteFile(filename.to_string()));
    }
  }

  if wtr.flush().is_err(){
    return Err(CluEError::CannotWriteFile(filename.to_string()));
  }

  Ok(())
}
//------------------------------------------------------------------------------
impl Add for &Signal{
  type Output = Signal;
 
  fn add(self,rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v + rhs.data[ii]);
    }
    Signal{data}
  }
}

impl Sub for &Signal{
  type Output = Signal;
 
  fn sub(self,rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v - rhs.data[ii]);
    }
    Signal{data}
  }
}


impl Mul for &Signal{
  type Output = Signal;
 
  fn mul(self, rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v*rhs.data[ii]);
    }
    Signal{data}
  }
}
impl Div for &Signal{
  type Output = Signal;
 
  fn div(self,rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v/rhs.data[ii]);
    }
    Signal{data}
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

    let signal2 = &signal0 + &signal1;
    assert_eq!(signal2.data,vec![2.0*ONE,1.5*ONE]);

    let signal3 = &signal0 - &signal1;
    assert_eq!(signal3.data,vec![0.0*ONE,0.5*ONE]);

    let signal4 = &signal0 * &signal1;
    assert_eq!(signal4.data,vec![ONE,0.5*ONE]);

    let signal5 = &signal0 / &signal1;
    assert_eq!(signal5.data,vec![ONE,2.0*ONE]);
  }
}


