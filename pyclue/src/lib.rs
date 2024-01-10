
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
  config: HashMap<String,String>,
  groups: HashMap::<String,Group>,
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
    
    for (_key, group) in self.groups.iter(){
      out = format!("{}{}\n", out, group.to_string());
    }

    for (_key, properties) in self.structure_properties.iter(){
      out = format!("{}{}\n", out, properties.to_string());
    }

    for (_key, properties) in self.spin_properties.iter(){
      out = format!("{}{}\n", out, properties.to_string());
    }

    out
  }
}
//------------------------------------------------------------------------------
#[pymethods] 
impl CluEOptions {

  #[new]
  #[args(config = "HashMap::<String,String>::new()",
     group_list = "Vec::<Group>::new()",
     structure_properties_list = "Vec::<StructureProperties>::new()",
     spin_properties_list = "Vec::<SpinProperties>::new()",
     )]
  fn new(config: HashMap::<String,String>,
    group_list: Vec::<Group>,
    structure_properties_list: Vec::<StructureProperties>,
    spin_properties_list: Vec::<SpinProperties>,
    ) -> Self 
  {

    let mut groups = HashMap::<String,Group>::with_capacity(
        group_list.len());
    for group in group_list{
      groups.insert(group.label.clone(),group);
    }

    let mut structure_properties = HashMap::<String,StructureProperties>::
      with_capacity(structure_properties_list.len());
    for structure_property in structure_properties_list{
      structure_properties.insert(structure_property.label.clone(),
          structure_property);
    }

    let mut spin_properties = HashMap::<(String,String),SpinProperties>::
      with_capacity(spin_properties_list.len());
    for spin_property in spin_properties_list{
      spin_properties.insert(
          (spin_property.label.clone(), spin_property.isotope.clone()),
          spin_property);
    }

    CluEOptions{
      config,
      groups,
      structure_properties,
      spin_properties,
    }
  }
  //----------------------------------------------------------------------------
  pub fn run(&self) -> PyResult<(Vec::<f64>, Vec::<Complex::<f64>>)>{
    let config = match clue::config::Config::from(&self.to_string()){
      Ok(cfg) => cfg,
      Err(err) => return Err(PyErr::new::<PyTypeError, _>(err.to_string())),
    };

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
  pub fn add_group(&mut self, group: Group) -> PyResult<()>{

    self.groups.insert(group.label.clone(), group);
    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn set_group(&mut self, label: String, key: String, value: String)
    -> PyResult<()>
  {
    if !self.groups.contains_key(&key){
      self.groups.insert(label.clone(), 
          Group::new(label.clone(), HashMap::<String,String>::new() ));
    }

    let Some(group) = self.groups.get_mut(&label) else{
      return Err(PyErr::new::<PyTypeError, _>("cannot find group"));
    };

    group.criteria.insert(key,value);

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
  pub fn set_structure_properties(&mut self, label: String,
      key: String, value: String)
    -> PyResult<()>
  {
    if !self.structure_properties.contains_key(&label){
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
  pub fn add_spin_properties(&mut self, 
      spin_properties: SpinProperties) -> PyResult<()>{

    self.spin_properties.insert(
        (spin_properties.label.clone(), spin_properties.isotope.clone()),
        spin_properties);
    Ok(())
  }
  //----------------------------------------------------------------------------
  pub fn set_spin_properties(&mut self, label: String, isotope: String,
      key: String, value: String)
    -> PyResult<()>
  {
    let spin_key = (label.clone(),isotope.clone());

    if !self.spin_properties.contains_key(&spin_key){
      self.spin_properties.insert((label.clone(), isotope.clone()), 
          SpinProperties::new(label.clone(), isotope.clone(),
            HashMap::<String,String>::new() ));
    }

    let Some(spin_properties) = self.spin_properties.get_mut( &spin_key ) 
    else{
      return Err(PyErr::new::<PyTypeError, _>(
            "cannot find spine_properties"));
    };

    spin_properties.properties.insert(key,value);

    Ok(())
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass]
#[derive(Clone)]
struct Group{
  pub label: String,
  criteria: HashMap::<String,String>
}
//------------------------------------------------------------------------------
impl ToString for Group{
  fn to_string(&self) -> String {
    let mut out = format!("#[group(label = {})]\n",self.label);

    for (key,val) in &self.criteria{
      out = format!("{}  {} {};\n",out,key,val);
    }
    out
  }
}
//------------------------------------------------------------------------------
#[pymethods] 
impl Group{
  #[new]
  #[args(criteria = "HashMap::<String,String>::new()")]
  pub fn new(label: String, criteria: HashMap::<String,String>) -> Self{
    Group{
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
      out = format!("{}  {} = {};\n",out,key,val);
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
      out = format!("{}  {} = {};\n",out,key,val);
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
    m.add_class::<Group>()?;
    m.add_class::<StructureProperties>()?;
    m.add_class::<SpinProperties>()?;
    Ok(())
}
}

