use crate::space_3d::Vector3D;

pub trait Translate{ fn translate(&mut self, r: &Vector3D); }

pub trait GetCentroid{fn centroid(&self) -> &Vector3D; }

pub trait GetIndices{fn indices(&self) -> Vec::<usize>; }

pub trait GetExchangeCoupling{fn exchange_coupling(&self) -> f64; }

#[derive(Debug,Clone)]
pub struct ExchangeGroupManager{
  pub exchange_groups: Vec::<ExchangeGroup>,
  pub exchange_group_ids: Vec::<Option<usize>>,
  pub exchange_couplings: Vec::<f64>, // one entry per exchange_group
}

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

impl Translate for ExchangeGroup{
  fn translate(&mut self, r: &Vector3D){
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.translate(r),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.translate(r),
    }
  }
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
#[derive(Debug, Clone)]
pub struct C3Rotor{
  pub center: Vector3D,
  pub normal: Vector3D,
  pub indices: [usize;3],
}

//------------------------------------------------------------------------------
impl ToString for C3Rotor{
  fn to_string(&self) -> String{
    format!("{}_{}_{}", self.indices[0],self.indices[1],self.indices[2])
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
                                                                                 
                                                                                 
  let mut normal = x.cross(&y);                                                  
                                                                                 
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


