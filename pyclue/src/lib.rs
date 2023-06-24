
#![allow(unused)]
fn main() {
use pyo3::prelude::*;
use pyo3::exceptions::PyTypeError;
use std::iter::zip;

use std::collections::HashMap;

use num_complex::Complex;

use clue_oxide as clue;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass]
struct CluEOptions{
  //base_config: String,
  config: HashMap<String,String>,
  filters: HashMap::<String,Filter>,
  structure_properties: HashMap<String,StructureProperties>,
  spin_properties: HashMap<(String,String),SpinProperties>,
}
//------------------------------------------------------------------------------
impl ToString for CluEOptions{
  fn to_string(&self) -> String{
    let mut out = String::from("#[config]\n");

    for (key, value) in self.config.iter(){
      out = format!("{}{} = {};\n",out, key,value);
    } 
    out = format!("{}\n",out);
    
    for (_key, filter) in self.filters.iter(){
      out = format!("{}{}", out, filter.to_string());
    }

    out
  }
}
//------------------------------------------------------------------------------
#[pymethods] 
impl CluEOptions {

  #[new]
  #[args(config = "HashMap::<String,String>::new()",
     filters = "HashMap::<String,Filter>::new()",
     structure_properties = "HashMap::<String,StructureProperties>::new()",
     spin_properties = "HashMap::<(String,String),SpinProperties>::new()",
     )]
  fn new(config: HashMap::<String,String>,
    filters: HashMap::<String,Filter>,
    structure_properties: HashMap::<String,StructureProperties>,
    spin_properties: HashMap<(String,String),SpinProperties>,
    ) -> Self 
  {
    CluEOptions{
      config,
      filters,
      structure_properties,
      spin_properties,
    }
  }
  //----------------------------------------------------------------------------
  pub fn run(&self) -> PyResult<(Vec::<f64>, Vec::<Complex::<f64>>)>{
    let expressions 
      = match clue::config::lexer::get_tokens_from_line(&self.to_string()){
        Ok(exps) => exps,
        Err(err) => return Err(PyErr::new::<PyTypeError, _>(err.to_string())),
      };

    let mut config = clue_oxide::config::Config::new();

    for expression in expressions.iter(){
      match config.parse_config_line(expression){
        Ok(_) => (),
        Err(err) => return Err(PyErr::new::<PyTypeError, _>(err.to_string())),
      }
    }
    match clue::run(config){
      Ok((time, signal)) => Ok((time,signal)),
      Err(err) => return Err(PyErr::new::<PyTypeError, _>(err.to_string())),
    }
  }
  //----------------------------------------------------------------------------
  pub fn print(&self) -> PyResult<()>{

    println!("{}",self.to_string());

    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn set_config(&mut self, key: String, value: String) -> PyResult<()> {
    self.config.insert(key,value);
    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn add_filter(&mut self, filter: Filter) -> PyResult<()>{

    self.filters.insert(filter.label.clone(), filter);
    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn set_filter(&mut self, label: String, key: String, value: String)
    -> PyResult<()>
  {
    if !self.filters.contains_key(&key){
      self.filters.insert(label.clone(), 
          Filter::new(label.clone(), HashMap::<String,String>::new() ));
    }

    let Some(filter) = self.filters.get_mut(&label) else{
      return Err(PyErr::new::<PyTypeError, _>("cannot find filter"));
    };

    filter.criteria.insert(key,value);

    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn add_structure_properties(&mut self, 
      structure_properties: StructureProperties) -> PyResult<()>{

    self.structure_properties.insert(structure_properties.label.clone(), 
        structure_properties);
    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn set_structure_properties(&mut self, label: String, key: String, value: String)
    -> PyResult<()>
  {
    if !self.structure_properties.contains_key(&key){
      self.structure_properties.insert(label.clone(), 
          StructureProperties::new(label.clone(), 
            HashMap::<String,String>::new() ));
    }

    let Some(structure_properties) = self.structure_properties.get_mut(&label) 
    else{
      return Err(PyErr::new::<PyTypeError, _>(
            "cannot find structure_properties"));
    };

    structure_properties.properties.insert(key,value);

    Ok(())
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass]
#[derive(Clone)]
struct Filter{
  pub label: String,
  criteria: HashMap::<String,String>
}
//------------------------------------------------------------------------------
impl ToString for Filter{
  fn to_string(&self) -> String {
    let mut out = format!("#[filter(label = {})]\n",self.label);

    for (key,val) in &self.criteria{
      out = format!("{}  {} {}\n",out,key,val);
    }
    out
  }
}
//------------------------------------------------------------------------------
#[pymethods] 
impl Filter{
  #[new]
  #[args(criteria = "HashMap::<String,String>::new()")]
  pub fn new(label: String, criteria: HashMap::<String,String>) -> Self{
    Filter{
      label,
      criteria,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass]
#[derive(Clone)]
struct StructureProperties{
  pub label: String,
  properties: HashMap::<String,String>
}
//------------------------------------------------------------------------------
impl ToString for StructureProperties{
  fn to_string(&self) -> String {
    let mut out = format!("#[structure_properties(label = {})]\n",self.label);

    for (key,val) in &self.properties{
      out = format!("{}  {} {}\n",out,key,val);
    }
    out
  }
}
//------------------------------------------------------------------------------
#[pymethods] 
impl StructureProperties{
  #[new]
  #[args(properties = "HashMap::<String,String>::new()")]
  pub fn new(label: String, properties: HashMap::<String,String>) -> Self{
    StructureProperties{
      label,
      properties,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass]
#[derive(Clone)]
struct SpinProperties{
  pub label: String,
  pub isotope: String,
  properties: HashMap::<String,String>
}
//------------------------------------------------------------------------------
impl ToString for SpinProperties{
  fn to_string(&self) -> String {
    let mut out = format!("#[spin_properties(label = {}, isotope = {})]\n",
        self.label, self.isotope);

    for (key,val) in &self.properties{
      out = format!("{}  {} {}\n",out,key,val);
    }
    out
  }
}
//------------------------------------------------------------------------------
#[pymethods] 
impl SpinProperties{
  #[new]
  #[args(properties = "HashMap::<String,String>::new()")]
  pub fn new(label: String, isotope: String,
      properties: HashMap::<String,String>) -> Self
  {
    SpinProperties{
      label,
      isotope,
      properties,
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn clue_oxide(_py: Python, m: &PyModule) -> PyResult<()> {
    //m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add_class::<CluEOptions>()?;
    m.add_class::<Filter>()?;
    m.add_class::<SpinProperties>()?;
    m.add_class::<StructureProperties>()?;
    //m.add_class::<SpinProperties>()?;
    Ok(())
}
}

