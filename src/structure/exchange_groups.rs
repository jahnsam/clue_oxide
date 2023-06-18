use crate::CluEError;
use crate::physical_constants::ANGSTROM;
use crate::space_3d::Vector3D;
//use crate::structure::Structure;

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub trait Translate{ fn translate(&mut self, r: &Vector3D); }

pub trait GetCentroid{fn centroid(&self) -> &Vector3D; }

pub trait GetNormal{fn normal(&self) -> &Vector3D; }

pub trait GetIndices{fn indices(&self) -> Vec::<usize>; }

pub trait GetExchangeCoupling{fn exchange_coupling(&self) -> f64; }

//pub trait ToStringResult{fn to_string_result(&self, structure: &Structure) 
//  -> Result<String, CluEError>;}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct ExchangeGroupManager{
  pub exchange_groups: Vec::<ExchangeGroup>,
  pub exchange_group_ids: Vec::<Option<usize>>,
  pub exchange_couplings: Vec::<f64>, // one entry per exchange_group
}

impl ExchangeGroupManager{
  pub fn to_csv(&self, filename: &str) -> Result<(),CluEError>
  {
    let Ok(file) = File::create(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };

    let n_char_f64 = 16;
    let bytes_per_char = 32;
    let n_bytes = 8*bytes_per_char*n_char_f64*(self.exchange_groups.len()+1);

    let mut stream = BufWriter::with_capacity(n_bytes,file);

    let line = "exchange_group,\
center_x,center_y,center_z,\
normal_x,normal_y,normal_z,\
exchange_coupling\n".to_string();
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    for (ii,exchange_group) in self.exchange_groups.iter().enumerate(){

      let methyl_str = exchange_group.to_string();
      let r = exchange_group.centroid().scale(1.0/ANGSTROM);
      let n = exchange_group.normal();
      let j = self.exchange_couplings[ii];

      let line = format!("{},{},{},{},{},{},{},{}\n",
          methyl_str, r.x(), r.y(), r.z(), n.x(),n.y(), n.z(), j);
    
      if stream.write(line.as_bytes()).is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }
    }


    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug, Clone)]
pub enum ExchangeGroup{
  Methyl(C3Rotor),
  PrimaryAmonium(C3Rotor),
}

impl ToString for ExchangeGroup{
  fn to_string(&self) -> String{
    match self{
      ExchangeGroup::Methyl(rotor) => format!("methyl_{}",rotor.to_string()),
      ExchangeGroup::PrimaryAmonium(rotor) 
        => format!("primary_amonium_{}",rotor.to_string()),
    }
  }
}
//------------------------------------------------------------------------------
impl GetCentroid for ExchangeGroup{
  fn centroid(&self) -> &Vector3D{
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.centroid(),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.centroid(),
    }
  }
}

impl GetIndices for ExchangeGroup{
  fn indices(&self) -> Vec::<usize>{
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.indices(),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.indices(),
    }
  }
}

impl GetNormal for ExchangeGroup{
  fn normal(&self) -> &Vector3D{
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.normal(),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.normal(),
    }
  }
}

impl Translate for ExchangeGroup{
  fn translate(&mut self, r: &Vector3D){
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.translate(r),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.translate(r),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
#[derive(Debug, Clone)]
pub struct C3Rotor{
  pub center: Vector3D,
  pub normal: Vector3D,
  pub indices: [usize;3],
}

//------------------------------------------------------------------------------
impl ToString for C3Rotor{
  // This function returns the reference indices of the rotor as a string.
  fn to_string(&self) -> String{
    // The output should be the reference indices, but the rotor stores the
    // structure indices, so a "+1" is used to convert from structure indices
    // to reference indices.
    format!("{}_{}_{}", self.indices[0]+1,self.indices[1]+1,self.indices[2]+1)
  }
}
//------------------------------------------------------------------------------
impl GetCentroid for C3Rotor{
  fn centroid(&self) -> &Vector3D {                                                
    &self.center                                                           
  }
}
impl C3Rotor{                                                                     
                                                                                 
pub fn from(r_carbon: Vector3D,                                                   
    h0: Vector3D,                                                                 
    h1: Vector3D,                                                                 
    h2: Vector3D,
    indices: [usize;3],    
    ) -> Self{                                                                 
                                                                                 
  let center = (&( &h0 + &h1) + &h2).scale(1.0/3.0);                             
                                                                                 
  let x = (&h0 - &center).normalize();                                           
                                                                                 
  let y = (&h1 - &center).normalize();                                           
  let y = &y - &x.scale(x.dot(&y));                                              
                                                                                 
                                                                                 
  let mut normal = x.cross(&y).normalize();                                                  
                                                                                 
  let axis = &center -  &r_carbon;                                               
  if normal.dot(&axis) < 0.0 {                                                   
    normal = normal.scale(-1.0);                                                 
  }                                                                              
                                                                                 
  C3Rotor{ center, normal, indices}
}                                                                                


pub fn normal(&self) -> &Vector3D {                                                
  &self.normal                                                      
}
}

impl GetIndices for C3Rotor{
  fn indices(&self) -> Vec::<usize>{
    
    let mut out = Vec::<usize>::with_capacity(3);

    for h in self.indices{
      out.push(h);
    }
    out
  }
}

impl Translate for C3Rotor{
  fn translate(&mut self, r: &Vector3D){
   self.center = &self.center + r;
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


