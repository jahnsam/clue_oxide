
#![allow(unused)]
fn main() {
use pyo3::prelude::*;
use std::iter::zip;

use std::collections::HashMap;

//use clue_oxide;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass]
struct CluEOptions{
  //base_config: String,
  config: HashMap<String,String>,
  filters: HashMap::<String,Filter>,
  //spin_properties: HashMap<(String,String),HashMap<String,String>>,
  //structure_properties: HashMap<String,HashMap<String,String>>,
}

#[pymethods] 
impl CluEOptions {

  #[new]
  #[args(config = "HashMap::<String,String>::new()")]
  fn new(config: HashMap::<String,String>) -> Self {
    CluEOptions{
      //base_config,
      config,
      filters: HashMap::<String,Filter>::new(),
      //spin_properties: HashMap<(String,String),HashMap<String,String>>::new(),
      //structure_properties: HashMap<String,HashMap<String,String>>::new(),
    }
  }
  //----------------------------------------------------------------------------
  pub fn print(&self) -> PyResult<()>{
    println!("#[config]");
    //println!("{}",self.base_config);

    for (key, value) in self.config.iter(){
      println!("{} = {};",key,value);
    } 

    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn set_config(&mut self, key: String, value: String) -> PyResult<()> {
    self.config.insert(key,value);
    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn set_filter(&mut self, label: String, key: String, value: String)
    -> PyResult<()>
  {
    if !self.filters.contains_key(){
      self.filters.insert(label.clone(), Filter::new(label.clone()));
    }

    let Some(filter) = self.filters.get_mut(&label) else{
      return Err(String::from("cannot find filter"))
    };

    filter

    Ok(())
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
struct Filter{
  label: String,
  fields: HashMap::<String,String>
}
impl Filter{
  pub fn new(label: String) -> Self{
    Filter{
      label,
      fields: HashMap::<String,String>::new(),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/// This function removes all points from a list that have a positive slope.
/*
#[pyfunction]
fn run(clue_options: &CluEOptions) 
  -> PyResult<Vec::<f64>,Vec::<Complex::<f64>>> {

  clue_oxide::run(config)
}
*/
/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn clue_oxide(_py: Python, m: &PyModule) -> PyResult<()> {
    //m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add_class::<CluEOptions>()?;
    Ok(())
}
}

